/*
DeBruijnGraphTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
import org.junit.jupiter.api.AfterEach;
import static org.junit.jupiter.api.Assertions.*;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.PrintStream;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.io.ByteArrayOutputStream;

class DeBruijnGraphTest {
    @Test 
    void DeBruijnGraphConstructor() {
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 2}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);
        assertNotNull(dbg);
        assertEquals(k, dbg.k);
        assertEquals(reads, dbg.reads);
        // test the effects of setupFromMap() called by the constructor
        assertEquals(firstRead.length(), dbg.readLen);
        assertNotNull(dbg.nodeValues);
        // 2 kmers (AAAAA,AAAAT) split into 1 k-1mer with two suffixes (AAAA-AAAA,AAAA-AAAT), 2 unique k-1mers overall
        assertEquals(2, dbg.numNodes);
        assertEquals(2, dbg.nodeValues.size());
        assertNotNull(dbg.nodeValueIndexMap);
        assertEquals(dbg.nodeValues.size(), dbg.nodeValueIndexMap.keySet().size());
        assertNotNull(dbg.adjacency);
        assertEquals(2, dbg.adjacency.length);
        assertEquals(2, dbg.adjacency[0].size());
        assertNotNull(dbg.adjacency[0]);
        assertTrue(dbg.adjacency[0].contains(0));
        assertTrue(dbg.adjacency[0].contains(1));
        assertTrue(dbg.adjacency[1] == null || dbg.adjacency[1].size() == 0); // TODO consider whether we want these outdegree==0 nodes to have null entries or empty list entries instead
        assertNotNull(dbg.adjacencyCoverageCount);
        assertEquals(2, dbg.adjacencyCoverageCount.length);
        assertNotNull(dbg.backAdjacency);
        assertEquals(2, dbg.backAdjacency.length);
        assertNotNull(dbg.backAdjacency[0]);
        assertTrue(dbg.backAdjacency[0].contains(0));
        assertNotNull(dbg.backAdjacency[1]);
        assertTrue(dbg.backAdjacency[1].contains(0));
        assertNotNull(dbg.indegree);
        assertEquals(2, dbg.indegree.length);
        assertEquals(1, dbg.indegree[0]);
        assertEquals(1, dbg.indegree[1]);
        assertNotNull(dbg.outdegree);
        assertEquals(2, dbg.outdegree.length);
        assertEquals(2, dbg.outdegree[0]);
        assertEquals(0, dbg.outdegree[1]);
        // todo check that debruijnMap has had memory freed
    }

    // TODO assertions on much more complicated artificial input data getting translated to the appropriate DeBruijnGraph

    @Test 
    void DeBruijnGraphConstructorThrowsIllegalArgumentExceptionIfDeBruijnMapIsNullOrEmpty() {
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));

        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = new DeBruijnGraph(null, reads, k);
            });
        assertEquals("DeBruijnGraph cannot be constructed from a null or empty debruijnMap", exception.getMessage());

        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = new DeBruijnGraph(new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>(), reads, k);
            });
        assertEquals("DeBruijnGraph cannot be constructed from a null or empty debruijnMap", exception2.getMessage());
    }

    @Test 
    void DeBruijnGraphConstructorThrowsIllegalArgumentExceptionIfReadsListIsNullOrEmpty() {
        int k = 5;
        // don't add any reads in this test
        List<byte[]> reads = new ArrayList<byte[]>();
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 2}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, null, k);
            });
        assertEquals("DeBruijnGraph cannot be constructed from a null or empty reads list", exception.getMessage());

        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);
            });
        assertEquals("DeBruijnGraph cannot be constructed from a null or empty reads list", exception2.getMessage());
    }

    @Test
    void DeBruijnGraphConstructorThrowsIllegalArgumentExceptionIfKIsLessThanZeroOrGreaterThanOrEqualToReadLength() {
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 2}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, 0);
            });
        assertEquals("DeBruijnGraph cannot be constructed with k values that are <= 0 or equal to the length of the read (should be substantially shorter than the read length)", exception.getMessage());

        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, 6);
            });
        assertEquals("DeBruijnGraph cannot be constructed with k values that are <= 0 or equal to the length of the read (should be substantially shorter than the read length)", exception2.getMessage());
    }

    @Test
    void getKMinus1MerString() {
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 2}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);
        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);

        // DeBruijnGraph.getKMinus1MerString should return the node kmer (k-1) string starting in the read specified by the first 3 bytes and the offset specified by the last byte
        assertEquals("AAAA", dbg.getKMinus1MerString(new byte[]{0, 0, 0, 0}));
        assertEquals("AAAT", dbg.getKMinus1MerString(new byte[]{0, 0, 0, 2}));
    }

    @Test
    void getKMinus1MerStringThrowsReadsNullExceptionIfReadsIsNull() {
        int k = 5;

        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));

        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 2}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);
        dbg.reads = null;

        DeBruijnGraph.ReadsNullException exception = assertThrows(DeBruijnGraph.ReadsNullException.class, () -> {
                dbg.getKMinus1MerString(new byte[]{0, 0, 0, 0});
        });
        assertEquals("This DeBruijnGraph reads instance is missing, cannot decode pointers into reads into kmer strings!", exception.getMessage());
    }

    @Test
    void getKMinus1MerStringThrowsReadIdxOutOfRangeExceptionIfReadIdxIsOutOfRange() {
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 2}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);

        DeBruijnGraph.ReadIdxOutOfRangeException exception = assertThrows(DeBruijnGraph.ReadIdxOutOfRangeException.class, () -> {
                dbg.getKMinus1MerString(new byte[]{1, 0, 0, 0});
        });
        assertEquals("DeBruijnGraph getKMinus1MerString called with invalid zero-based read index (1) (# of reads: 1)", exception.getMessage());
    }

    @Test
    void getKMinus1MerStringThrowsOffsetOutOfRangeExceptionIfOffsetIsOutOfRange() {
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes 
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 2}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);

        DeBruijnGraph.OffsetOutOfRangeException exception = assertThrows(DeBruijnGraph.OffsetOutOfRangeException.class, () -> {
                dbg.getKMinus1MerString(new byte[]{0, 0, 0, 6});
        });
        assertEquals("DeBruijnGraph getKMinus1MerString called with invalid zero-based offset (6) (read length: 6)", exception.getMessage());
    }

    @Test
    void saveToDiskThrowsIllegalArgumentExceptionIfAttemptingToSaveWithoutExpectedExtension() {
        DeBruijnGraph dbg = setupDeBruijnGraphForTests();
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph.saveToDisk(dbg, "deBruijnGraph");
        });
        assertEquals("DeBruijnGraph.saveToDisk will only attempt to save '" + DeBruijnGraph.DBG_FILE_EXTENSION + "' extension files", exception.getMessage());
    }

    @Test
    void saveToDiskReturnsFalseIfCannotSave() throws Exception {
        String invalidDbgFilename = "directory_not_existing/invalidDbg.dbg";
        DeBruijnGraph dbg = setupDeBruijnGraphForTests();
        assertFalse(DeBruijnGraph.saveToDisk(dbg, invalidDbgFilename));
    }

    private final PrintStream standardErr = System.err;
    @AfterEach
    public void tearDown() {
        System.setErr(standardErr);
    }

    @Test
    void clearFromDiskLogsAStacktraceButDoesNotThrowAnExceptionIfRequestedToDeleteAnInvalidPath() {
        final ByteArrayOutputStream toCaptureOutputStream = new ByteArrayOutputStream();
        System.setErr(new PrintStream(toCaptureOutputStream));
        DeBruijnGraph.clearFromDisk("not_existing_directory/ADeBruijnGraphThatDoesNotExist.dbg");
        String output = toCaptureOutputStream.toString();
        assertNotNull(output);
        assertEquals("java.nio.file.NoSuchFileException: not_existing_directory/ADeBruijnGraphThatDoesNotExist.dbg", output.split("\n")[0]);
    }

    @Test
    void clearFromDiskThrowsIllegalArgumentExceptionIfAttemptingToClearFileWithoutExpectedExtension() {
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph.clearFromDisk("deBruijnGraph");
        });
        assertEquals("DeBruijnGraph.clearFromDisk will only remove '" + DeBruijnGraph.DBG_FILE_EXTENSION + "' extension files", exception.getMessage());
    }

    @Test
    void loadFromDiskThrowsIllegalArgumentExceptionIfAttemptingToLoadFileWithoutExpectedExtension() {
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph.loadFromDisk("deBruijnGraph");
        });
        assertEquals("DeBruijnGraph.loadFromDisk will only attempt to load '" + DeBruijnGraph.DBG_FILE_EXTENSION + "' extension files", exception.getMessage());
    }

    @Test
    void loadFromDiskReturnsNullIfFileCannotBeDeserialized() {
        String invalidDbgFilename = "invalidDbg.dbg";
        try {
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(invalidDbgFilename)));
            out.println("Not a serialized dbg.");
            out.flush();
            DeBruijnGraph dbg = DeBruijnGraph.loadFromDisk(invalidDbgFilename);
            assertNull(dbg);
        } catch (Exception e) {
            fail("Threw an exception when expected it to return null");
        } finally {
            try {
                Files.deleteIfExists(FileSystems.getDefault().getPath(invalidDbgFilename));
            } catch (java.io.IOException e) {
                e.printStackTrace();
            }
        }
    }

    // TODO replace setup in preceding unit tests with call to here when possible
    // right now this is used by other classes that need a pre-initialized DeBruijnGraph with expected fake data for their own unit/integration testing purposes
    public static DeBruijnGraph setupDeBruijnGraphForTests() {
        int k = 5;
        String firstRead = "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 2}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);
        assertNotNull(dbg);
        assertEquals(k, dbg.k);
        assertEquals(reads, dbg.reads);
        // test the effects of setupFromMap() called by the constructor
        assertEquals(firstRead.length(), dbg.readLen);
        assertNotNull(dbg.nodeValues);
        // 2 kmers (AAAAA,AAAAT) split into 1 k-1mer with two suffixes (AAAA-AAAA,AAAA-AAAT), 2 unique k-1mers overall
        assertEquals(2, dbg.numNodes);

        return dbg;
    }

    public static DeBruijnGraph setupDeBruijnGraphForTestsSimpleSingleNodeCycle() {
        int k = 5;
        String firstRead = "AAAAAA"; 
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);
        assertNotNull(dbg);
        assertEquals(k, dbg.k);
        assertEquals(reads, dbg.reads);
        // test the effects of setupFromMap() called by the constructor
        assertEquals(firstRead.length(), dbg.readLen);
        assertNotNull(dbg.nodeValues);
        assertEquals(1, dbg.numNodes);

        return dbg;
    }

    public static DeBruijnGraph getSimpleBubbleGraphBasedOnPhiX174First70bpWith1Error() {
        // first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
        // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
        //                                     A <- on one read we introduce a point error at this position, which causes a low coverage bubble in the de Bruijn graph
        //                                     *
        int k = 15;
        int readLen = 30;
        List<byte[]> reads = new ArrayList<byte[]>();

        String phix174sampleSequence = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        int coverageCount = 10;
        int readStep = 5;
        for (int coverageRepeatIdx = 0; coverageRepeatIdx < coverageCount; coverageRepeatIdx++) {
            for (int i = 0; i < phix174sampleSequence.length() - readLen; i+=readStep) {
                String read = phix174sampleSequence.substring(i, i+readLen);
                if (read.equals("CGCAGAAGTTAACACTTTCGGATATTTCTG") && coverageRepeatIdx == 4) {
                    read = "CGCAGAAGTTAACACTATCGGATATTTCTG";
                    //                      *
                }
                reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(read));
            }
        }

        DeBruijnKmers dbk = new DeBruijnKmers();
        dbk.disableTrieCorrections = true; // for these tests we want to leave the errors in
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        dbk.processKmers();
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = dbk.getDeBruijnMap();
        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);

        return dbg;
    }

    public static DeBruijnGraph getSimpleIncomingBubbleGraphBasedOnPhiX174First70bpWith1ErrorNearStart() {
        // first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
        // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
        //        T <- on one read we introduce a point error at this position within k bases of the start, which causes a low coverage incoming bubble in the de Bruijn graph
        //        *
        int k = 15;
        int readLen = 30;
        List<byte[]> reads = new ArrayList<byte[]>();

        String phix174sampleSequence = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        int coverageCount = 10;
        int readStep = 5;
        for (int coverageRepeatIdx = 0; coverageRepeatIdx < coverageCount; coverageRepeatIdx++) {
            for (int i = 0; i < phix174sampleSequence.length() - readLen; i+=readStep) {
                String read = phix174sampleSequence.substring(i, i+readLen);
                if (read.equals("GAGTTTTATCGCTTCCATGACGCAGAAGTT") && coverageRepeatIdx == 4) {
                    read = "GAGTTTTTTCGCTTCCATGACGCAGAAGTT";
                    //             *
                }
                reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(read));
            }
        }

        DeBruijnKmers dbk = new DeBruijnKmers();
        dbk.disableTrieCorrections = true; // for these tests we want to leave the errors in
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        dbk.processKmers();
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = dbk.getDeBruijnMap();
        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);

        return dbg;
    }

    public static DeBruijnGraph getEccDNA80DeBruijnGraphWithOrWithoutCircularOverlapInReads(boolean withCircularOverlap) {
        int k = 15;
        int readLen = 30;
        int coverageCount = 5;
        int readStep = 5;
        return getEccDNA80DeBruijnGraphWithOrWithoutCircularOverlapInReads(withCircularOverlap, k, readLen, coverageCount, readStep);
    }

    public static DeBruijnGraph getEccDNA80DeBruijnGraphWithOrWithoutCircularOverlapInReads(boolean withCircularOverlap, int k, int readLen, int coverageCount, int readStep) {
        List<String> readStrings = DeBruijnGenomeAssemblerTest.getEccDNA80Reads(withCircularOverlap, k, readLen, coverageCount, readStep);
        List<byte[]> reads = new ArrayList<byte[]>();
        for (String readString : readStrings) {
            reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(readString));
        }

        DeBruijnKmers dbk = new DeBruijnKmers();
        dbk.disableTrieCorrections = true; // for these tests we want to leave the errors in
        dbk.setupFromReadStringsAndSpecifiedK(k, reads);
        dbk.processKmers();
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = dbk.getDeBruijnMap();
        DeBruijnGraph dbg = new DeBruijnGraph(debruijnMap, reads, k);

        return dbg;
    }
}
