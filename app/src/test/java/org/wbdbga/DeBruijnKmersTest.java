/*
DeBruijnKmersTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
Copyright (C) 2025 William Bruns

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

package org.wbdbga;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import java.util.*;

class DeBruijnKmersTest {
    @Test
    void DeBruijnDefaultConstructor() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        assertNotNull(dbk);
    }

    @Test
    void setupFromReadStringsAndSpecifiedKThrowsAnExceptionIfTheReadsAreNull() {
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnKmers dbk = new DeBruijnKmers();
                dbk.setupFromReadStringsAndSpecifiedK(30, null);
            });
        assertEquals("reads list must not be null", exception.getMessage());
    }

    @Test
    void setupFromReadStringsAndSpecifiedKThrowsAnExceptionIfTheReadsAreEmpty() {
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnKmers dbk = new DeBruijnKmers();
                dbk.setupFromReadStringsAndSpecifiedK(30, new ArrayList<byte[]>());
            });
        assertEquals("reads list must not be empty", exception.getMessage());
    }

    @Test
    void setupFromReadStringsAndSpecifiedKIfFirstReadLengthIsLessThanKThrowsAnException() {
        int k = 30;
        String firstRead = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 29 As when K = 30
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnKmers dbk = new DeBruijnKmers();
                List<byte[]> reads = new ArrayList<byte[]>();
                reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead)); 
                dbk.setupFromReadStringsAndSpecifiedK(k, reads);
            });
        assertEquals("k (==" + k + ") must be <= read length (== " + firstRead.length() + ")", exception.getMessage());
    }

    @Test
    void setupFromReadStringsAndSpecifiedK() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        assertNull(dbk.kmerCoverageMap);
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new LinkedList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        assertNull(dbk.reads);
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        assertEquals(reads, dbk.reads);
        assertTrue(reads instanceof LinkedList);
        // should have been turned into an arraylist for constant time by index access
        assertTrue(dbk.reads instanceof ArrayList);
        assertEquals(k, dbk.k);
        assertNotNull(dbk.kmerCoverageMap);
        assertEquals(2, dbk.kmerCoverageMap.keySet().size());
        assertEquals(Byte.valueOf((byte)1), dbk.kmerCoverageMap.get(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}))); // "AAAAA" (readIdx 0 in first 3 bytes, offset 0 in last byte), covered by 1 read
        assertEquals(Byte.valueOf((byte)1), dbk.kmerCoverageMap.get(ValueComparedByteArray.from(new byte[]{0, 0, 0, 1}))); // "AAAAT" (readIdx 0 in first 3 bytes, offset 1 in last byte), covered by 1 read
    }

    @Test
    void setupFromReadStringsAndSpecifiedKWithUseTrieCorrectionsForSmallGenomesEnabled() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        dbk.useTrieCorrectionsForSmallGenomes = true;
        assertNull(dbk.kmerCoverageMap);
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        assertNull(dbk.reads);
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        assertEquals(reads, dbk.reads);
        assertEquals(k, dbk.k);
        assertNotNull(dbk.kmerCoverageMap);
        // should NOT make read corrections, same result as previous unit test
        assertEquals(2, dbk.kmerCoverageMap.keySet().size());
        assertEquals(Byte.valueOf((byte)1), dbk.kmerCoverageMap.get(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}))); // "AAAAA" (readIdx 0 in first 3 bytes, offset 0 in last byte), covered by 1 read
        assertEquals(Byte.valueOf((byte)1), dbk.kmerCoverageMap.get(ValueComparedByteArray.from(new byte[]{0, 0, 0, 1}))); // "AAAAT" (readIdx 0 in first 3 bytes, offset 1 in last byte), covered by 1 read
    }

    // TODO add integration test where setupFromReadStringsAndSpecifiedK triggers KmerTrieFilterReads.recommendReadCorrections to return read corrections

    @Test
    void loadRead() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        assertNull(dbk.kmerCoverageMap);
        int k = 5;
        String baseRead = "TTTTTT"; // 0th read - present during `dbk.setupFromReadStringsAndSpecifiedK(k, reads);` call
        String firstRead = "AAAAAA"; // 1st read - NOT present during `dbk.setupFromReadStringsAndSpecifiedK(k, reads);` call during this test specifically
        String expectedKmer = "AAAAA";
        byte[] expectedKmerReadAndIndexPairRaw = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(1, 0);
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(baseRead));
        // in this test we don't add firstRead here before calling setupFromReadStringsAndSpecifiedK
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        assertNull(dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        int addlCoverage = dbk.loadRead(1);
        assertNotNull(dbk.kmerCoverageMap);
        assertEquals(2, addlCoverage); // AAAAA (Ax5) is covered 2x in AAAAAA (Ax6), at positions 0 and 1, though they have the same key deduplicated to (read 1, offset 0)
        assertEquals(2, dbk.kmerCoverageMap.keySet().size()); // 1 key for TTTTT, 1 key for AAAAA
        assertNotNull(dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
        assertArrayEquals(expectedKmerReadAndIndexPairRaw, dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
    }

    @Test
    void loadReadReturnsZeroForInvalidReadsThatCauseExceptionsInGetStringFromBPBytes() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        List<byte[]> reads = new ArrayList<byte[]>();
        int k = 30;
        reads.add(new byte[]{5,0}); // invalid
        dbk.reads = reads;
        int result = dbk.loadRead(0);
        assertEquals(0, result); // no exception expected
    }

    @Test
    void unloadRead() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        assertNull(dbk.kmerCoverageMap);
        int k = 5;
        String baseRead = "TTTTTT"; // 0th read - present during `dbk.setupFromReadStringsAndSpecifiedK(k, reads);` call
        String firstRead = "AAAAAA"; // 1st read - present during `dbk.setupFromReadStringsAndSpecifiedK(k, reads);` call
        String expectedKmer = "AAAAA";
        byte[] expectedKmerReadAndIndexPairRaw = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(1, 0);
        List reads = new ArrayList();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(baseRead));
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        dbk.freeKmerHashCodeToReadIndexAndSubstringIndexMapAtEndOfSetup = false;
        dbk.disableTrieCorrections = true;
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        // note this would be freed up at the end of setupFromReadStringsAndSpecifiedK, but we disabled that via freeKmerHashCodeToReadIndexAndSubstringIndexMapAtEndOfSetup
        assertNotNull(dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
        assertEquals(2, dbk.kmerCoverageMap.keySet().size()); // 1 key for TTTTT, 1 key for AAAAA
        int coverageRemoved = dbk.unloadRead(1);
        assertEquals(2, coverageRemoved); // 2x AAAAA (Ax5) in AAAAAA (Ax6)
        // note this would be freed up at the end of setupFromReadStringsAndSpecifiedK, but we disabled that via freeKmerHashCodeToReadIndexAndSubstringIndexMapAtEndOfSetup
        //assertNull(dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
        //System.out.println("dbk.kmerCoverageMap.keySet(): " + dbk.kmerCoverageMap.keySet());
        // TODO: not sure I want there to be a fake read added during the unloading process in this case, review this more closely
        assertEquals(1, dbk.kmerCoverageMap.keySet().size()); // 1 key for TTTTT
        assertNull(dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
    }

    @Test
    void unloadReadWhenUnloadedReadIsTheDefiningReadButNoEntryInKmerHashCodeToReadIndexAndSubstringIndexMapProceedsWithoutFailing() {
        // same setup as the last unit test
        DeBruijnKmers dbk = new DeBruijnKmers();
        assertNull(dbk.kmerCoverageMap);
        int k = 5;
        String baseRead = "TTTTTT"; // 0th read - present during `dbk.setupFromReadStringsAndSpecifiedK(k, reads);` call
        String firstRead = "AAAAAA"; // 1st read - present during `dbk.setupFromReadStringsAndSpecifiedK(k, reads);` call
        String expectedKmer = "AAAAA";
        byte[] expectedKmerReadAndIndexPairRaw = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(1, 0);
        List reads = new ArrayList();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(baseRead));
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        dbk.freeKmerHashCodeToReadIndexAndSubstringIndexMapAtEndOfSetup = false;
        dbk.disableTrieCorrections = true;
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        // note this would be freed up at the end of setupFromReadStringsAndSpecifiedK, but we disabled that via freeKmerHashCodeToReadIndexAndSubstringIndexMapAtEndOfSetup
        assertNotNull(dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
        assertEquals(2, dbk.kmerCoverageMap.keySet().size()); // 1 key for TTTTT, 1 key for AAAAA
        // note that collisions are expected and allowed in the hashCode keyed map and handled manually by having the value be an array (we cannot have a keySet of kmer strings due to excess memory requirement but need constant time lookup by string to dedupe across reads the kmer strings), when lookups occur, the actual kmer key is always compared, see the implementation of getKmerReadAndIndexPairForKmerString , or the existing unit tests like getKmerReadAndIndexPairForKmerStringExistingEntriesForHashcode
        assertEquals(2, dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.keySet().size());
        // then we remove this key from kmerHashCodeToReadIndexAndSubstringIndexMap (there is not a collision in this case so we can do it this way in this unit test)
        // this may seem artificial but we are checking that this method does not fail if something else drops the low coverage kmer directly,
        dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.remove("AAAAA".hashCode()); // we can remove "AAAAA" as if the low coverage kmer had already been dropped and also unload that read
        int coverageRemoved = dbk.unloadRead(1);
        assertEquals(2, coverageRemoved); // 2x AAAAA (Ax5) in AAAAAA (Ax6) 
        // So again, more coverage can be reported as removed in the current implementation than is registered as existing because we allow dropping of low coverage kmers AND unloading of them again with the reads. If kmers are not dropped by e.g. removeLowCoverageKmers, and unloadReads is only called 1x per read, then coverage removed by unloadReads will balance with updates to read coverage in kmerCoverageMap, otherwise it will not necessarily balance as no such balance is assumed (making sure readers do not assume it).
        assertEquals(1, dbk.kmerCoverageMap.keySet().size()); // 1 key for TTTTT
        assertNull(dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
    }

    @Test
    void unloadReadCanUpdateCoverageForAReadWhenNotRemovingTheDefiningReadOrOnlyReadForThatKmer() {
        // same setup as the last unit test
        DeBruijnKmers dbk = new DeBruijnKmers();
        assertNull(dbk.kmerCoverageMap);
        int k = 5;
        String baseRead = "TTTTTT"; // 0th read
        String firstRead = "AAAAAA"; // 1st read
        // 2nd read, means "AAAAA" is covered by multiple reads,
        // we will in this test unload the 2nd read which is not where "AAAAA" was originally defined
        String secondRead = "AAAAAA";
        String expectedKmer = "AAAAA";
        byte[] expectedKmerReadAndIndexPairRaw = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(1, 0);
        List reads = new ArrayList();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(baseRead));
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(secondRead));
        dbk.freeKmerHashCodeToReadIndexAndSubstringIndexMapAtEndOfSetup = false;
        dbk.disableTrieCorrections = true;
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);        
        assertEquals(2, dbk.kmerCoverageMap.keySet().size()); // 1 key for TTTTT, 1 key for AAAAA
        int coverageRemoved = dbk.unloadRead(2);
        // adding assertions here shortly
        assertEquals(2, coverageRemoved); // 2x AAAAA (Ax5) in AAAAAA (Ax6)
        // however in this case this kmer should still have read coverage 2, and still be present
        assertEquals(2, dbk.kmerCoverageMap.keySet().size()); // 1 key for TTTTT, 1 key for AAAAA, compare to unloadRead() test where this was reduced to 1 key
        assertNotNull(dbk.getKmerReadAndIndexPairForKmerString(expectedKmer));
    }

    @Test
    void getKmerReadAndIndexPairForKmerStringNoExistingEntries() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        byte[] kmerReadAndIndexPair = dbk.getKmerReadAndIndexPairForKmerString("AAA");
        assertNull(kmerReadAndIndexPair);
    }

    @Test
    void getKmerReadAndIndexPairForKmerStringThrowsNotInitializedExceptionIfNotInitializedWhenCalled() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        // cannot use putKmerReadAndIndexPairForKmerString because its error handling will stop us from adding it
        // manually add
        String kmer = "AAA";
        byte[] kmerReadAndIndexPairRaw = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(0, 1);
        List<ValueComparedByteArray> entries = new LinkedList<ValueComparedByteArray>();
        entries.add(ValueComparedByteArray.from(kmerReadAndIndexPairRaw));
        // note hashCode is used since the memory savings of avoiding having kmer strings as keys is enormous, but it means that the method under test does it own collision handling (tested in a separate unit test)
        dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.put(kmer.hashCode(), entries);
        // finished adding manually
        DeBruijnKmers.NotInitializedException exception = assertThrows(DeBruijnKmers.NotInitializedException.class, () -> {
                byte[] kmerReadAndIndexPairRetrieved = dbk.getKmerReadAndIndexPairForKmerString("AAA");
            });
        assertEquals("getKmerReadAndIndexPairForKmerString has been called on a DeBruijnKmers instance that has not been properly initialized; please ensure setupFromReadStringsAndSpecifiedK is called first and that DeBruijnKmers.reads is not being set to null", exception.getMessage());
    }

    @Test
    void putKmerReadAndIndexPairForKmerStringThrowsNotInitializedExceptionIfNotInitialized() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        DeBruijnKmers.NotInitializedException exception = assertThrows(DeBruijnKmers.NotInitializedException.class, () -> {
                byte[] kmerReadAndIndexPair = dbk.putKmerReadAndIndexPairForKmerString("AAA", 0, 1);        
            });
        assertEquals("putKmerReadAndIndexPairForKmerString has been called on a DeBruijnKmers instance that has not been properly initialized; please ensure setupFromReadStringsAndSpecifiedK is called first and that DeBruijnKmers.reads is not being set to null", exception.getMessage());
    }

    @Test
    void getKmerReadAndIndexPairForKmerStringThrowsCorruptedReadIndexAndSubstringIndexMapExceptionIfThereIsAnOutOfBoundsReadIndex() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        dbk.reads = new ArrayList<byte[]>();
        String firstRead = "TAAA";
        // the read needs to actually exist to avoid an exception during putKmerReadAndIndexPairForKmerString now
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        DeBruijnKmers.CorruptedReadIndexAndSubstringIndexMapException exception = assertThrows(DeBruijnKmers.CorruptedReadIndexAndSubstringIndexMapException.class, () -> {
                // entry is there, but no reads list exists which is needed to verify the retrieved item to avoid being confused by hash collisions
                byte[] kmerReadAndIndexPair = dbk.putKmerReadAndIndexPairForKmerString("AAA", 0, 1);
                // cannot put an invalid read index above due to new error checking on putKmerReadAndIndexPairForKmerString
                List<ValueComparedByteArray> entries = dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.get("AAA".hashCode());
                assertEquals(1, entries.size()); // in this test case we only put one entry in, we are not testing the collision handling tested in a separate case
                byte[] corrupted = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(1, 1);
                entries.get(0).value[0] = corrupted[0]; // corrupt the read index directly
                entries.get(0).value[1] = corrupted[1]; // corrupt the read index directly
                entries.get(0).value[2] = corrupted[2]; // corrupt the read index directly
                entries.get(0).value[3] = corrupted[3]; // corrupt the read index directly
                byte[] kmerReadAndIndexPairRetrieved = dbk.getKmerReadAndIndexPairForKmerString("AAA");
        });
        assertEquals("getKmerReadAndIndexPairForKmerString found an invalid entry referring to a read index that is outside the range of this.reads, this should not ever happen.", exception.getMessage());
    }

    @Test
    void getKmerReadAndIndexPairForKmerStringThrowsCorruptedReadIndexAndSubstringIndexMapExceptionIfTheOffsetIsOutOfBounds() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        dbk.reads = new ArrayList<byte[]>();
        String firstRead = "TAAA";
        // the read needs to actually exist so getKmerReadAndIndexPairForKmerString can doublecheck what it retrieved since it uses hashCode to avoid maintenance of a large keyset in this memory constrained setting
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        DeBruijnKmers.CorruptedReadIndexAndSubstringIndexMapException exception = assertThrows(DeBruijnKmers.CorruptedReadIndexAndSubstringIndexMapException.class, () -> {
                byte[] kmerReadAndIndexPair = dbk.putKmerReadAndIndexPairForKmerString("AAA", 0, 1); // ok offset on insertion to avoid triggering exception tested separately in put
                // cannot put an invalid offst above due to new error checking on putKmerReadAndIndexPairForKmerString
                List<ValueComparedByteArray> entries = dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.get("AAA".hashCode());
                assertEquals(1, entries.size()); // in this test case we only put one entry in, we are not testing the collision handling tested in a separate case
                byte[] corrupted = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(0, 2);
                entries.get(0).value[0] = corrupted[0]; // corrupt the offset directly
                entries.get(0).value[1] = corrupted[1]; // corrupt the offset directly
                entries.get(0).value[2] = corrupted[2]; // corrupt the offset directly
                entries.get(0).value[3] = corrupted[3]; // corrupt the offset directly
                byte[] retrievedKmerReadAndIndexPair = dbk.getKmerReadAndIndexPairForKmerString("AAA");
        });
        assertEquals("getKmerReadAndIndexPairForKmerString found an invalid entry, with an offset which overflows the read length, this should not ever happen.", exception.getMessage());
    }

    @Test
    void getKmerReadAndIndexPairForKmerStringHasEntryManuallyAdded() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        dbk.reads = new ArrayList<byte[]>();
        String firstRead = "TAAA";
        // the read needs to actually exist so getKmerReadAndIndexPairForKmerString can doublecheck what it retrieved since it uses hashCode to avoid maintenance of a large keyset in this memory constrained setting
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        // manually add
        String kmer = "AAA";
        byte[] kmerReadAndIndexPairRaw = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(0, 1);
        List<ValueComparedByteArray> entries = new LinkedList<ValueComparedByteArray>();
        entries.add(ValueComparedByteArray.from(kmerReadAndIndexPairRaw));
        // note hashCode is used since the memory savings of avoiding having kmer strings as keys is enormous, but it means that the method under test does it own collision handling (tested in a separate unit test)
        dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.put(kmer.hashCode(), entries);
        // finished adding manually
        byte[] kmerReadAndIndexPair = dbk.getKmerReadAndIndexPairForKmerString(kmer);
        assertArrayEquals(kmerReadAndIndexPairRaw, kmerReadAndIndexPair);
    }

    @Test
    void getKmerReadAndIndexPairForKmerStringHasEntryAddedViaPutKmerReadAndIndexPairForKmerString() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        String firstRead = "TAAA";
        dbk.reads = new ArrayList<byte[]>();
        // the read needs to actually exist so getKmerReadAndIndexPairForKmerString can double check what it retrieved since it uses hashCode to avoid maintenance of a large keyset in this memory constrained setting
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        assertNull(dbk.getKmerReadAndIndexPairForKmerString("AAA"));
        byte[] kmerReadAndIndexPair = dbk.putKmerReadAndIndexPairForKmerString("AAA", 0, 1);
        byte[] retrievedKmerReadAndIndexPair = dbk.getKmerReadAndIndexPairForKmerString("AAA");
        assertNotNull(retrievedKmerReadAndIndexPair);
        assertArrayEquals(kmerReadAndIndexPair, retrievedKmerReadAndIndexPair);
    }

    @Test
    void getKmerReadAndIndexPairForKmerStringExistingEntriesForHashcode() {
        // e.g. "TTAATGCTTATTTTTTTTAAAGTAACTAAATTATTTACAA" 
        // and "AAAAATAATTTTTTAAAGTTTAACTTTATTAAAAATAAAA" for example both have hashCode 1456205320
        String kmer1 = "TTAATGCTTATTTTTTTTAAAGTAACTAAATTATTTACAA";
        String kmer2 = "AAAAATAATTTTTTAAAGTTTAACTTTATTAAAAATAAAA";
        // make sure that our test data is still valid , that these share a hashcode
        assertEquals(kmer1.length(), kmer2.length());
        assertEquals(kmer1.hashCode(), kmer2.hashCode());
        assertNotEquals(kmer1, kmer2);
        // ok, then we have two kmers with the same hashcode which are not equal

        // construct a read featuring these:
        // 100bp read that includes both kmers with kmer1 at offset 0 and kmer2 at offset 50
        String firstRead = kmer1 + "ACTGACTGAC" + kmer2 + "TGACTGACTG";

        DeBruijnKmers dbk = new DeBruijnKmers();
        dbk.k = kmer1.length(); // not actually used in this test since we don't test setupFromReadStringsAndSpecifiedK here but to avoid possible concerns of readers
        dbk.reads = new ArrayList<byte[]>();
        // the read needs to actually exist so getKmerReadAndIndexPairForKmerString can double check what it retrieved since it uses hashCode to avoid maintenance of a large keyset in this memory constrained setting
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        // we did not initialize dbk with `dbk.setupFromReadStringsAndSpecifiedK(k, reads);`, so neither kmer should be loaded yet:
        assertNull(dbk.getKmerReadAndIndexPairForKmerString(kmer1));
        assertNull(dbk.getKmerReadAndIndexPairForKmerString(kmer2));
        // now we load them
        byte[] kmerReadAndIndexPair1 = dbk.putKmerReadAndIndexPairForKmerString(kmer1, 0, 0);
        byte[] kmerReadAndIndexPair2 = dbk.putKmerReadAndIndexPairForKmerString(kmer2, 0, 50);
        byte[] retrievedKmerReadAndIndexPair1 = dbk.getKmerReadAndIndexPairForKmerString(kmer1);
        byte[] retrievedKmerReadAndIndexPair2 = dbk.getKmerReadAndIndexPairForKmerString(kmer2);
        assertArrayEquals(kmerReadAndIndexPair1, retrievedKmerReadAndIndexPair1);
        assertArrayEquals(kmerReadAndIndexPair2, retrievedKmerReadAndIndexPair2);
        assertEquals(ValueComparedByteArray.from(kmerReadAndIndexPair1), ValueComparedByteArray.from(retrievedKmerReadAndIndexPair1));
        assertEquals(ValueComparedByteArray.from(kmerReadAndIndexPair2), ValueComparedByteArray.from(retrievedKmerReadAndIndexPair2));
        assertNotEquals(ValueComparedByteArray.from(kmerReadAndIndexPair1), ValueComparedByteArray.from(kmerReadAndIndexPair2));
        assertNotEquals(ValueComparedByteArray.from(retrievedKmerReadAndIndexPair1), ValueComparedByteArray.from(retrievedKmerReadAndIndexPair2));
        // confirm they indeed were stored under the same key, though this is an implementation detail
        // maybe change in implementation could make this test accidentally invalid due to other k value overriding
        assertNotNull(dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmer1.hashCode()));
        assertNotNull(dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmer2.hashCode()));
        assertEquals(2, dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmer1.hashCode()).size());
        assertEquals(2, dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmer2.hashCode()).size());
    }

    // putKmerReadAndIndexPairForKmerString
    // already partially tested above in tests that check both put and get
    // - no entry in putKmerReadAndIndexPairForKmerString for the kmer hashCode -> return new entry
    // - already an entry in putKmerReadAndIndexPairForKmerString for the kmer hashCode but no entry for this kmer -> return new entry
    // - already an entry for this specific kmer -> return existing entry

    // we can test input validation here
    // also add validation that we do not enter read indices that are out of bounds or offsets that are out of bounds
    @Test
    void putKmerReadAndIndexPairForKmerStringThrowsExceptionIfReadIndexExceedsAvailableReads() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        String firstRead = "TAAA"; // length 4 read
        dbk.reads = new ArrayList<byte[]>();
        // the read needs to actually exist so getKmerReadAndIndexPairForKmerString can double check what it retrieved since it uses hashCode to avoid maintenance of a large keyset in this memory constrained setting
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        DeBruijnKmers.OutOfRangeReadIndexException exception = assertThrows(DeBruijnKmers.OutOfRangeReadIndexException.class, () -> {
                byte[] kmerReadAndIndexPair = dbk.putKmerReadAndIndexPairForKmerString("AAA", 1, 1); // intentional mistake here tries to insert at read index 1
            });
        assertEquals("putKmerReadAndIndexPairForKmerString has been called with an invalid entry referring to a read index that is outside the range of this.reads, this should not ever happen.", exception.getMessage());
    }

    @Test
    void putKmerReadAndIndexPairForKmerStringThrowsExceptionIfOffsetPlusKmerLengthExceedsReadLength() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        String firstRead = "TAAA"; // length 4 read
        dbk.reads = new ArrayList<byte[]>();
        // the read needs to actually exist so getKmerReadAndIndexPairForKmerString can double check what it retrieved since it uses hashCode to avoid maintenance of a large keyset in this memory constrained setting
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        DeBruijnKmers.OutOfRangeReadOffsetException exception = assertThrows(DeBruijnKmers.OutOfRangeReadOffsetException.class, () -> {
                byte[] kmerReadAndIndexPair = dbk.putKmerReadAndIndexPairForKmerString("AAA", 0, 2); // intentional mistake here tries to insert at offset 2
            });
        assertEquals("putKmerReadAndIndexPairForKmerString has been called with an invalid entry referring to a offset that denotes a kmer extending beyond the range of the read, this should not ever happen.", exception.getMessage());
    }

    @Test
    void getDeBruijnMapThrowsExceptionIfCalledWhenDeBruijnMapIsNotYetSet() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        DeBruijnKmers.NotInitializedException exception = assertThrows(DeBruijnKmers.NotInitializedException.class, () -> {
                Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = dbk.getDeBruijnMap();
            });
        assertEquals("debruijnMap not yet created. call setupFromReadStringsAndSpecifiedK and processKmers before calling getDebruijnMap()", exception.getMessage());
    }

    @Test
    void getDeBruijnMapReturnsDeBruijnMapIfAvailable() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        String firstRead = "TAAA";
        dbk.reads = new ArrayList<byte[]>();
        dbk.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        assertNull(dbk.getKmerReadAndIndexPairForKmerString("AAA"));
        byte[] kmerReadAndIndexPair = dbk.putKmerReadAndIndexPairForKmerString("AAA", 0, 1);
        dbk.debruijnMap.put(ValueComparedByteArray.from(kmerReadAndIndexPair), new HashMap<ValueComparedByteArray, Byte>());
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = dbk.getDeBruijnMap();
        assertEquals(dbk.debruijnMap, debruijnMap);
    }

    @Test
    void processKmersThrowsNotInitializedExceptionIfCalledBeforeKmerCoverageMapIsReady() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        DeBruijnKmers.NotInitializedException exception = assertThrows(DeBruijnKmers.NotInitializedException.class, () -> {
                dbk.processKmers();
            });
        assertEquals("kmerCoverageMap is not yet loaded. setupFromReadStringsAndSpecifiedK must be called before this method.", exception.getMessage());
    }

    /*
    // TODO add validation that processKmers is not called more than once
    */
    @Test
    void processKmers() {
        DeBruijnKmers dbk = new DeBruijnKmers();
        assertNull(dbk.kmerCoverageMap);
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        // should replace kmerHashCodeToReadIndexAndSubstringIndexMap with (k-1)-mers prefix/suffix of the original kmers 
        // and populate debruijnMap
        assertTrue(dbk.debruijnMap == null || dbk.debruijnMap.keySet().size() == 0);
        assertNotNull(dbk.kmerCoverageMap);
        assertEquals(2, dbk.kmerCoverageMap.keySet().size()); // "AAAAA", "AAAAT" // 2 kmers
        dbk.processKmers();
        // debruijnMap should now be populated
        assertNotNull(dbk.debruijnMap);
        assertTrue(dbk.debruijnMap.keySet().size() > 0);
        assertEquals(1, dbk.debruijnMap.keySet().size()); // only 1 (k-1)mer in this case // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        assertEquals(2, dbk.debruijnMap.get(dbk.debruijnMap.keySet().iterator().next()).size());
        // should have cleared kmerHashCodeToReadIndexAndSubstringIndexMap to free memory
        assertTrue(dbk.kmerHashCodeToReadIndexAndSubstringIndexMap == null || dbk.kmerHashCodeToReadIndexAndSubstringIndexMap.keySet().size() == 0);
    }

    @Test
    void applyReadCorrectionSimpleCase() {
DeBruijnKmers dbk = new DeBruijnKmers();
        assertNull(dbk.kmerCoverageMap);
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);

        KmerTrieFilterReads.ReadCorrection rc = new KmerTrieFilterReads.ReadCorrection(0, 5, 'A');
        assertEquals((byte)1, dbk.kmerCoverageMap.getOrDefault(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), (byte)0)); // AAAAA , covered 1x
        int netReadCoverageRemoved = dbk.applyReadCorrection(rc);
        assertEquals(0, netReadCoverageRemoved); // removed 1x coverage from AAAAT, added 1x coverage to AAAAA, balances out
        assertEquals((byte)2, dbk.kmerCoverageMap.getOrDefault(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), (byte)0)); // AAAAA , now covered 2x
        /*
          started loadKmerCoverageMapFromReads
          Updating count for AAAAA ([0, 0, 0, 0]) to 1
          Updating count for AAAAT ([0, 0, 0, 1]) to 1
          Apply read correction:
          Updating count for AAAAA ([0, 0, 0, 0]) to 0
          Updating count for AAAAT ([0, 0, 0, 1]) to 0
          Read coverage removed: 2
          Updating count for AAAAA ([0, 0, 0, 0]) to 1
          Updating count for AAAAA ([0, 0, 0, 0]) to 2
          Read coverage added: 2
          Added and removed coverage, should balance
        */
    }
}
