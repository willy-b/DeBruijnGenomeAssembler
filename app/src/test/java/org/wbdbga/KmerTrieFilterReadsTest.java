/*
KmerTrieFilterReadsTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;

class KmerTrieFilterReadsTest {
    @Test
    void readCorrectionConstructor() {
        int readIdx = 0;
        int withinReadIdx = 1;
        char newCharacter = 'C';
        KmerTrieFilterReads.ReadCorrection rc = new KmerTrieFilterReads.ReadCorrection(readIdx, withinReadIdx, newCharacter);
        assertEquals(readIdx, rc.readIdx);
        assertEquals(withinReadIdx, rc.withinReadIdx);
        assertEquals(newCharacter, rc.newCharacter);
    }

    @Test
    void readCorrectionConstructorThrowsInvalidCharacterExceptionForInvalidCharacters() {
        int readIdx = 0;
        int withinReadIdx = 1;
        char newCharacter = 'F';
        
        KmerTrieFilterReads.ReadCorrection.InvalidCharacterException exception = assertThrows(KmerTrieFilterReads.ReadCorrection.InvalidCharacterException.class, () -> {
                KmerTrieFilterReads.ReadCorrection rc = new KmerTrieFilterReads.ReadCorrection(readIdx, withinReadIdx, newCharacter);
            });
        assertEquals("ReadCorrections must specify a new character of 'A','T','C', or 'G', found 'F'", exception.getMessage());
    }

    @Test
    public void KmerTrieFilterReadsConstructor() {
        List<byte[]> bpEncodedReadsList = new ArrayList<byte[]>();
        Set<ValueComparedByteArray> kmers = new HashSet<ValueComparedByteArray>();
        int k = 35;
        int numFakeReadsAddedForOrphanedKmers = 1;

        KmerTrieFilterReads ktfr = new KmerTrieFilterReads(bpEncodedReadsList, kmers, k, numFakeReadsAddedForOrphanedKmers);
        assertEquals(bpEncodedReadsList, ktfr.reads);
        assertEquals(kmers, ktfr.kmers);
        assertEquals(k, ktfr.k);
        assertEquals(numFakeReadsAddedForOrphanedKmers, ktfr.numFakeReadsAddedForOrphanedKmers);
        assertNotNull(ktfr.trie);
    }

    @Test
    public void KmerTrieFilterReadsConstructorThrowsExceptionIfReadsArgumentIsNull() {
        Set<ValueComparedByteArray> kmers = new HashSet<ValueComparedByteArray>();
        int k = 35;
        int numFakeReadsAddedForOrphanedKmers = 1;

        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                KmerTrieFilterReads ktfr = new KmerTrieFilterReads(null, kmers, k, numFakeReadsAddedForOrphanedKmers);
            });
        assertEquals("reads argument must not be null", exception.getMessage());
    }

    @Test
    public void KmerTrieFilterReadsConstructorThrowsExceptionIfKIsInvalid() {
        List<byte[]> bpEncodedReadsList = new ArrayList<byte[]>();
        Set<ValueComparedByteArray> kmers = new HashSet<ValueComparedByteArray>();
        int numFakeReadsAddedForOrphanedKmers = 1;

        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                KmerTrieFilterReads ktfr = new KmerTrieFilterReads(bpEncodedReadsList, kmers, 0, numFakeReadsAddedForOrphanedKmers);
            });
        assertEquals("Invalid value of k", exception.getMessage());
    }

    @Test
    public void KmerTrieFilterReadsConstructorThrowsExceptionIfKmersArgumentIsNull() {
        List<byte[]> bpEncodedReadsList = new ArrayList<byte[]>();
        Set<ValueComparedByteArray> kmers = new HashSet<ValueComparedByteArray>();
        int k = 35;
        int numFakeReadsAddedForOrphanedKmers = 1;

        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                KmerTrieFilterReads ktfr = new KmerTrieFilterReads(bpEncodedReadsList, null, k, numFakeReadsAddedForOrphanedKmers);
            });
        assertEquals("kmers argument must not be null", exception.getMessage());
    }

    public KmerTrieFilterReads setupKmerTrieFilterReadsWithExampleDataForTests() {
        List<byte[]> bpEncodedReadsList = new ArrayList<byte[]>();
        int k = 5;
        int readLen = 10;

        // here is a 22bp short sequence from the 1st 100bp of the Candidatus Carsonella Ruddii genome for example
        // ATGAATACTATATTTTCAAGAA
        // let us sample 10bp reads that we will get 5-mers from

        // TODO replace by loop here as well
        // error free reads, we include 1 copy here for kmers to reference but only feature "solid" kmers which are covered by many other reads (not shown in the test)
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("ATGAATACTA")); // 0 = readIdx
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("AATACTATAT")); // 1
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("ACTATATTTT")); // 2
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("ATATTTTCAA")); // 3
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("TTTTCAAGAA")); // 4
        // reads with errors, these will have errors reported
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("ATGAAGACTA")); // 5 = readIdx
        /*
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("CAGACTATAT")); // 6
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("ACTATAGTTT")); // 7
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("ATATGTTCAA")); // 8
        bpEncodedReadsList.add(GenomeSequenceEncodingUtil.getBPBytesFromString("CTTTCAAGAA")); // 9
        */
        
        // Kmers (k=5) , presumed filtered for "solid" kmers with threshold coverage in caller method (not responsibility of this class to do thresholding for solid kmers, happens in DeBruijnKmers)
        // Can generate with a loop but also good to illustrate by writing them out
        // "ATGAA"*, "TGAAT"*, "GAATA"*, "AATAC"*, "ATACT"*, "TACTA"*, 0   readIdx (we don't reference the reads from 5 onwards in defining the solid kmers as in this example we have put the reads with all solid kmers in the first 5 indices)
        // "AATAC", "ATACT", "TACTA", "ACTAT"*, "CTATA"*, "TATAT"*,    1      |
        // "ACTAT", "CTATA", "TATAT", "ATATT"*, "TATTT"*, "ATTTT"*,    2      v
        // "ATATT", "TATTT", "ATTTT", "TTTTC"*, "TTTCA"*, "TTCAA"*,    3
        // "TTTTC", "TTTCA", "TTCAA", "TCAAG"*, "CAAGA"*, "AAGAA"*     4
        //      0        1        2        3         4         5       <- offsetIdx

        Set<ValueComparedByteArray> kmers = new HashSet<ValueComparedByteArray>();

        // the reads are spaced by 3 positions, there are 5 reads of length 10 (the subsequent reads feature kmers that are not solid, are reads with errors; in real data these are not separated but there are more overlapping/shifted copies of the reads 0-4 and few copies of what is in 5-9, so upstream method filters out the error driven kmers from 5-9 as not solid)
        // kmers are deduplicated upstream to refer to their earliest read only, so no duplicate kmers are entered
        for (int offsetIdx = 0; offsetIdx < readLen - (k - 1); offsetIdx++) {
            for (int readIdx = 0; readIdx < 5; readIdx++) {
                if (readIdx == 0 || offsetIdx > 2) {
                    //System.out.println("In unit test, adding kmer " + GenomeSequenceEncodingUtil.getStringFromBPBytes(bpEncodedReadsList.get(readIdx)).substring(offsetIdx, offsetIdx+k));
                    kmers.add(ValueComparedByteArray.from(GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(readIdx, offsetIdx)));
                }
            }
        }

        int numFakeReadsAddedForOrphanedKmers = 0;
        KmerTrieFilterReads ktfr = new KmerTrieFilterReads(bpEncodedReadsList, kmers, k, numFakeReadsAddedForOrphanedKmers);
        return ktfr;
    }

    @Test
    public void buildTrie() {
        KmerTrieFilterReads ktfr = setupKmerTrieFilterReadsWithExampleDataForTests();
        ktfr.buildTrie();
        assertNotNull(ktfr.trie);

        // root -> A, T, C, G
        assertTrue(ktfr.trie.get(0).containsKey('A'));
        assertTrue(ktfr.trie.get(0).containsKey('T'));
        assertTrue(ktfr.trie.get(0).containsKey('C'));
        assertTrue(ktfr.trie.get(0).containsKey('G'));

        // root -> A -> A only two paths AATAC or AAGAA
        assertTrue(ktfr.trie.get(ktfr.trie.get(0).get('A')).containsKey('A'));
        int aaIdx = ktfr.trie.get(ktfr.trie.get(0).get('A')).get('A');
        assertFalse(ktfr.trie.get(aaIdx).containsKey('A'));
        assertTrue(ktfr.trie.get(aaIdx).containsKey('T')); // AA*T*AC
        assertFalse(ktfr.trie.get(aaIdx).containsKey('C'));
        assertTrue(ktfr.trie.get(aaIdx).containsKey('G')); // AA*G*AA

        // AA -> GAA  AAGAA
        int aagIdx = ktfr.trie.get(aaIdx).get('G');
        assertTrue(ktfr.trie.get(aagIdx).containsKey('A'));
        assertFalse(ktfr.trie.get(aagIdx).containsKey('T'));
        assertFalse(ktfr.trie.get(aagIdx).containsKey('C'));
        assertFalse(ktfr.trie.get(aagIdx).containsKey('G'));
        int aagaIdx = ktfr.trie.get(aagIdx).get('A');
        assertTrue(ktfr.trie.get(aagaIdx).containsKey('A'));
        assertFalse(ktfr.trie.get(aagaIdx).containsKey('T'));
        assertFalse(ktfr.trie.get(aagaIdx).containsKey('C'));
        assertFalse(ktfr.trie.get(aagaIdx).containsKey('G'));

        // AA -> TAC AATAC
        int aatIdx = ktfr.trie.get(aaIdx).get('T');
        assertTrue(ktfr.trie.get(aatIdx).containsKey('A')); // AAT*A*C
        assertFalse(ktfr.trie.get(aatIdx).containsKey('T'));
        assertFalse(ktfr.trie.get(aatIdx).containsKey('C'));
        assertFalse(ktfr.trie.get(aatIdx).containsKey('G'));
        int aataIdx = ktfr.trie.get(aatIdx).get('A');
        assertFalse(ktfr.trie.get(aataIdx).containsKey('A'));
        assertFalse(ktfr.trie.get(aataIdx).containsKey('T'));
        assertTrue(ktfr.trie.get(aataIdx).containsKey('C')); // AATA*C*
        assertFalse(ktfr.trie.get(aataIdx).containsKey('G'));
    }

    @Test
    public void readSpectralAlignmentResultConstructor() {
        int readIdx = 20;
        Map<Integer, Integer> brokenKmersPerWithinReadIdx = new HashMap<Integer, Integer>();
        Map<Integer, Set<Character>> resolvingCharactersByPosition = new HashMap<Integer, Set<Character>>();
        KmerTrieFilterReads.ReadSpectralAlignmentResult rsar = new KmerTrieFilterReads.ReadSpectralAlignmentResult(readIdx, brokenKmersPerWithinReadIdx, resolvingCharactersByPosition);
        assertEquals(readIdx, rsar.readIdx);
        assertEquals(brokenKmersPerWithinReadIdx, rsar.brokenKmersPerWithinReadIdx);
        assertEquals(resolvingCharactersByPosition, rsar.resolvingCharactersByPosition);
    }

    @Test
    public void readSpectralAlignmentResultConstructorThrowsAnExceptionIfArgumentsAreInvalid() {
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                KmerTrieFilterReads.ReadSpectralAlignmentResult rsar = new KmerTrieFilterReads.ReadSpectralAlignmentResult(-1, new HashMap<Integer, Integer>(), new HashMap<Integer, Set<Character>>());
            });
        assertEquals("readIdx argument must be non-negative", exception.getMessage());

        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                KmerTrieFilterReads.ReadSpectralAlignmentResult rsar = new KmerTrieFilterReads.ReadSpectralAlignmentResult(1, null, new HashMap<Integer, Set<Character>>());
            });
        assertEquals("brokenKmersPerWithinReadIdx argument must be non-null", exception2.getMessage());

        IllegalArgumentException exception3 = assertThrows(IllegalArgumentException.class, () -> {
                KmerTrieFilterReads.ReadSpectralAlignmentResult rsar = new KmerTrieFilterReads.ReadSpectralAlignmentResult(1, new HashMap<Integer, Integer>(), null);
            });
        assertEquals("resolvingCharactersByPosition argument must be non-null", exception3.getMessage());
    }

    @Test
    public void checkReadBeforeInitializedThrowsException() {
        KmerTrieFilterReads ktfr = setupKmerTrieFilterReadsWithExampleDataForTests();
        // but intentionally do NOT call buildTrie()
        KmerTrieFilterReads.NotInitializedException exception = assertThrows(KmerTrieFilterReads.NotInitializedException.class, () -> {
                KmerTrieFilterReads.ReadSpectralAlignmentResult result = ktfr.checkRead(0);
            });
        assertEquals("ensure buildTrie is called before attempting to check or repair reads", exception.getMessage());
    }

    @Test
    public void checkReadWithNegativeReadIdxThrowsException() {
        KmerTrieFilterReads ktfr = setupKmerTrieFilterReadsWithExampleDataForTests();
        ktfr.buildTrie();
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                KmerTrieFilterReads.ReadSpectralAlignmentResult result = ktfr.checkRead(-1);
            });
        assertEquals("readIdx must be non-negative", exception.getMessage());
    }

    @Test
    public void checkReadWithOutOfBoundsIndexThrowsException() {
        KmerTrieFilterReads ktfr = setupKmerTrieFilterReadsWithExampleDataForTests();
        ktfr.buildTrie();
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                KmerTrieFilterReads.ReadSpectralAlignmentResult result = ktfr.checkRead(ktfr.reads.size());
            });
        assertEquals("readIdx must be within the bounds of this.reads", exception.getMessage());
    }

    @Test
    public void checkRead() {
        KmerTrieFilterReads ktfr = setupKmerTrieFilterReadsWithExampleDataForTests();
        ktfr.buildTrie();
        assertNotNull(ktfr.trie);

        // expect no errors for the error-free reads
        KmerTrieFilterReads.ReadSpectralAlignmentResult result1 = ktfr.checkRead(1); // skip index 0 as all of its indices are represented in kmer definitions so it is trivially error free
        assertEquals(1, result1.readIdx);
        assertEquals(0, result1.numSitesCausingMisalignedKmers);
        assertEquals(0, result1.numMisalignedKmers);
        assertEquals(0, result1.brokenKmersPerWithinReadIdx.keySet().size());

        KmerTrieFilterReads.ReadSpectralAlignmentResult result2 = ktfr.checkRead(2);
        assertEquals(2, result2.readIdx);
        assertEquals(0, result2.numSitesCausingMisalignedKmers);
        assertEquals(0, result2.numMisalignedKmers);
        assertEquals(0, result2.brokenKmersPerWithinReadIdx.keySet().size());

        KmerTrieFilterReads.ReadSpectralAlignmentResult result3 = ktfr.checkRead(3);
        assertEquals(3, result3.readIdx);
        assertEquals(0, result3.numSitesCausingMisalignedKmers);
        assertEquals(0, result3.numMisalignedKmers);
        assertEquals(0, result3.brokenKmersPerWithinReadIdx.keySet().size());

        KmerTrieFilterReads.ReadSpectralAlignmentResult result4 = ktfr.checkRead(4);
        assertEquals(4, result4.readIdx);
        assertEquals(0, result4.numSitesCausingMisalignedKmers);
        assertEquals(0, result4.numMisalignedKmers);
        assertEquals(0, result4.brokenKmersPerWithinReadIdx.keySet().size());

        // next the reads with errors, these should report issues        
        KmerTrieFilterReads.ReadSpectralAlignmentResult result5 = ktfr.checkRead(5);
        // read 5 has "ATGAAGACTA" where it should be "ATGAATACTA" (G->T in 0-based index 5)
        assertEquals(5, result5.readIdx);
        /*
        Starting trie check of "ATGAA":
        ATGAA

        Starting trie check of "TGAAG":
        TGAAG (no matching solid kmer) (last G not found, index 5 in 0-based index)

        Starting trie check of "GAAGA":
        GAAG (no matching solid kmer) (last G not found, index 5 in 0-based index)

        Starting trie check of "AAGAC":
        AAGAC (no matching solid kmer) (last C  not found, index 7 in 0-based index)

        Starting trie check of "AGACT":
        AG (no matching solid kmer)  (G not found, index 5 in 0-based index)

        Starting trie check of "GACTA":
        GAC (no matching solid kmer) (C not found, index 7 in 0-based index)
        {5=3, 7=2} // num breaks by 0-based index into this read
        {5=[T], 7=[A]} // characters that could first order resolve the break by 0-based index into this read; correcting index 5 to 'T' instead of 'G' would in fact fix the error
        // why correct reads and recompute the kmer spectrum iteratively rather than just using the initial solid spectrum? because the coverage ratios are used in later error correction and the iterative process can make more rare kmers solid
        */
        assertEquals(5, result5.numMisalignedKmers); // ideally each point error introduces k broken kmers (k=5 here), ideally all at the same position, but here 2 are shifted, we take the worst position in terms of kmer breaks to avoid correcting the wrong position
        assertEquals(2, result5.numSitesCausingMisalignedKmers);
        assertEquals(2, result5.brokenKmersPerWithinReadIdx.keySet().size());
        assertEquals(2, result5.resolvingCharactersByPosition.keySet().size());
        assertNotNull(result5.resolvingCharactersByPosition.get(5));
        assertEquals(1, result5.resolvingCharactersByPosition.get(5).size());
        assertEquals('T', result5.resolvingCharactersByPosition.get(5).iterator().next());
        assertNotNull(result5.resolvingCharactersByPosition.get(7));
        assertEquals(1, result5.resolvingCharactersByPosition.get(7).size());
        assertEquals('A', result5.resolvingCharactersByPosition.get(7).iterator().next());
    }

    @Test
    public void recommendReadCorrections() {
        KmerTrieFilterReads ktfr = setupKmerTrieFilterReadsWithExampleDataForTests();
        ktfr.buildTrie();
        assertNotNull(ktfr.trie);

        Set<KmerTrieFilterReads.ReadCorrection> rcs = ktfr.recommendReadCorrections();
        assertEquals(1, rcs.size());
        KmerTrieFilterReads.ReadCorrection rc = rcs.iterator().next();
        assertEquals(5, rc.readIdx);
        assertEquals(5, rc.withinReadIdx);
        assertEquals('T', rc.newCharacter);
    }

    @Test
    void cleanup() {
        KmerTrieFilterReads ktfr = setupKmerTrieFilterReadsWithExampleDataForTests();
        ktfr.buildTrie();
        assertNotNull(ktfr.trie);
        assertNotNull(ktfr.kmers);
        assertNotNull(ktfr.reads);
        assertNotNull(ktfr.spectrallyMisalignedReadsToIgnore);
        assertNotNull(ktfr.readCorrectionsToMake);

        ktfr.cleanup();

        assertNull(ktfr.trie);
        assertNull(ktfr.kmers);
        assertNull(ktfr.reads);
        assertNull(ktfr.spectrallyMisalignedReadsToIgnore);
        assertNull(ktfr.readCorrectionsToMake);
    }
}
