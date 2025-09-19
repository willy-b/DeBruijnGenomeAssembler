/*
KmerSizeSelectorTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.io.ByteArrayInputStream;

class KmerSizeSelectorTest {
    @Test
    public void KmerSizeSelectorConstructor() {
        KmerSizeSelector kss = new KmerSizeSelector();
        assertNotNull(kss);
    }

    @Test
    public void getDeBruijnGraphThrowsSolutionNotReadyExceptionIfCalledBeforeSolve() {
        KmerSizeSelector kss = new KmerSizeSelector();
        KmerSizeSelector.SolutionNotReadyException exception = assertThrows(KmerSizeSelector.SolutionNotReadyException.class, () -> {
                DeBruijnGraph dbg = kss.getDeBruijnGraph();
            });
        assertEquals("KmerSizeSelector.solve() should be called before getDeBruijnGraph can be used to retrieve the solution", exception.getMessage());
    }

    // TODO add unit test for getDeBruijnGraph when solution is ready (after solve has been called), will add after adding unit tests for solve

    @Test
    public void setReadsThrowsVariableReadLengthNotYetSupportedException() {
        KmerSizeSelector kss = new KmerSizeSelector();
        KmerSizeSelector.VariableReadLengthNotYetSupportedException exception = assertThrows(KmerSizeSelector.VariableReadLengthNotYetSupportedException.class, () -> {
                List<String> reads = new ArrayList<String>();
                reads.add("AAAAA");
                reads.add("AAAA");
                kss.setReads(reads);
            });
        assertEquals("Variable read lengths are disabled pending adequate support. Please preprocess data to uniform length (e.g. exclude reads below a minimum length and trim longer than typical reads to a common supported length.)", exception.getMessage());
    }

    @Test
    public void setReads() {
        KmerSizeSelector kss = new KmerSizeSelector();
        List<String> reads = new ArrayList<String>();
        reads.add("AAAAA");
        reads.add("AAAAT");
        assertEquals(2, reads.size());
        kss.setReads(reads);
        assertNotNull(kss.reads);
        assertEquals(reads.size(), kss.reads.size());
        assertArrayEquals(GenomeSequenceEncodingUtil.getBPBytesFromString(reads.get(0)), kss.reads.get(0));
        assertArrayEquals(GenomeSequenceEncodingUtil.getBPBytesFromString(reads.get(1)), kss.reads.get(1));
    }

    public static DeBruijnGraph getSimpleBalancedGraph() {
        int k = 5;
        String firstRead = "AAAAAA";
        // AAAA-AAAA
        List<byte[]> reads = new ArrayList<byte[]>();
        reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(firstRead));
        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

        Map<ValueComparedByteArray, Byte> kminus1merEdgeMap = new HashMap<ValueComparedByteArray, Byte>();
        kminus1merEdgeMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), Byte.valueOf((byte)1));
        debruijnMap.put(ValueComparedByteArray.from(new byte[]{0, 0, 0, 0}), kminus1merEdgeMap);

        return new DeBruijnGraph(debruijnMap, reads, k);
    }

    // just making sure behavior on this edge case does not change unexpectedly by asserting it
    // it is not a meaningful case
    @Test
    public void getN50ExtremelySimpleCase1EdgeButNoEdgesInContigNotCoveringCase() {
        DeBruijnGraph dbg = getSimpleBalancedGraph();
        // only 1 edge in graph, and 1 contig with 1 (k-1)mer node in it since we don't follow self-edges
        KmerSizeSelector kss = new KmerSizeSelector();
        int n50estimate = kss.getN50(dbg);
        // it doesn't cover the sequence since there is no edge retained in the contig (1 node in graph with self-loop is not followed)
        assertEquals(-1, n50estimate);
    }

    @Test
    public void getN50ExtremelySimple2EdgeCase() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        KmerSizeSelector kss = new KmerSizeSelector();
        int n50estimate = kss.getN50(dbg);
        // contigs AAAA (self-loop does not get followed), AAAAT (1 edge) ; also AAAAT is length 5 and covers 50% of the edges (1 out of 2), also >=50% of the sequence "AAAAAT" length 6
        assertEquals(1 + 5 - 1, n50estimate);
    }

    @Test
    void solveSimpleCase() {
        // first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
        // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
        //                                     A <- on one read we introduce a point error at this position, which causes a low coverage bubble in the de Bruijn graph
        //                                     *
        int k = 15;
        int readLen = 30;
        List<String> reads = new ArrayList<String>();

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
                reads.add(read);
            }
        }

        KmerSizeSelector kmerSizeSelector = new KmerSizeSelector();
        kmerSizeSelector.startingKmerSizeGuess = k;
        kmerSizeSelector.setReads(reads);
        int resultK = -1;
        try {
            resultK = kmerSizeSelector.solve();
        } catch (KmerSizeSelector.NoSuitableKmerLengthFoundException e) {
            fail("No suitable kmer found", e);
        }
        DeBruijnGraph dbg = kmerSizeSelector.getDeBruijnGraph();
        assertNotNull(dbg);
        assertTrue(dbg.k >= k);
        assertEquals(dbg.k, resultK);
    }

    @Test
    void solveThrowsNoSuitableKmerLengthFoundException() {
        int k = 15;

        List<String> reads = new ArrayList<String>();
        reads.add("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        reads.add("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

        KmerSizeSelector kmerSizeSelector = new KmerSizeSelector();
        kmerSizeSelector.startingKmerSizeGuess = k;
        kmerSizeSelector.setReads(reads);

        KmerSizeSelector.NoSuitableKmerLengthFoundException exception = assertThrows(KmerSizeSelector.NoSuitableKmerLengthFoundException.class, () -> {
                int resultK = kmerSizeSelector.solve();
        });
        assertEquals("no suitable kmer length found, currently we will not attempt to assemble if we do not think we can even assemble half of the underlying sequence, there may be an option to relax this in the future", exception.getMessage());
    }

    // credit to https://www.baeldung.com/java-junit-testing-system-in
    // see also https://docs.oracle.com/en/java/javase/21/docs/api/java.base/java/lang/System.html#setIn(java.io.InputStream)
    void overrideSystemInWithString(String simulatedStdinString) {
        ByteArrayInputStream simulatedStdinIS = new ByteArrayInputStream(simulatedStdinString.getBytes());
        System.setIn(simulatedStdinIS);
    }

    @Test
    public void readInParametersCanAcceptReadsViaStdin() {
        KmerSizeSelector kmerSizeSelector = new KmerSizeSelector();
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead);
        int k = 15;
        kmerSizeSelector.readInParameters(k);
        assertNotNull(kmerSizeSelector.reads);
        assertEquals(1, kmerSizeSelector.reads.size());
        assertEquals(firstRead, GenomeSequenceEncodingUtil.getStringFromBPBytes(kmerSizeSelector.reads.get(0)));
    }
    
    @Test
    public void readInParametersCanAcceptANumberIndicatingHowManyLinesToReadFirst() {
        KmerSizeSelector kmerSizeSelector = new KmerSizeSelector();
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        String secondReadShouldBeIgnored = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        overrideSystemInWithString("1\n" + firstRead + "\n" + secondReadShouldBeIgnored);
        int k = 15;
        kmerSizeSelector.readInParameters(k);
        assertNotNull(kmerSizeSelector.reads);
        assertEquals(1, kmerSizeSelector.reads.size());
        assertEquals(firstRead, GenomeSequenceEncodingUtil.getStringFromBPBytes(kmerSizeSelector.reads.get(0)));
    }

    @Test
    public void readInParametersThrowsExceptionWhenReadLinesAreNonUniformLength() {
        KmerSizeSelector kmerSizeSelector = new KmerSizeSelector();
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        String secondReadTooLong = firstRead + "T";
        overrideSystemInWithString(firstRead + "\n" + secondReadTooLong);
        int k = 15;
        KmerSizeSelector.VariableReadLengthNotYetSupportedException exception = assertThrows(KmerSizeSelector.VariableReadLengthNotYetSupportedException.class, () -> {
                kmerSizeSelector.readInParameters(k);
            });
        assertEquals("Variable read lengths are disabled pending adequate support. Please preprocess data to uniform length (e.g. exclude reads below a minimum length and trim longer than typical reads to a common supported length.)", exception.getMessage());
    }
}
