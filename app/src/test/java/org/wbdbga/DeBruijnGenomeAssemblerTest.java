/*
DeBruijnGenomeAssemblerTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.FileSystems;
import java.nio.file.Files;

class DeBruijnGenomeAssemblerTest {
    @Test 
    void deBruijnAssemblyOfDeBruijnGenomeAssemblerTestConstructor() {
        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler();
        assertNotNull(assembler);
        assertNotNull(assembler.kSelector);
    }

    @Test
    void solveSimpleCase() {
        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler();

        // first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
        // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
        //                                     A <- on one read we introduce a point error at this position, which causes a low coverage bubble in the de Bruijn graph
        //                                     *
        int k = 15;
        int readLen = 30;
        List<String> reads = new ArrayList<String>();

        String phix174sampleSequence = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        int coverageCount = 5;
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

        // bypass readInParameters call
        assembler.kSelector.startingKmerSizeGuess = k;
        assembler.kSelector.setReads(reads);

        // run solve and check the output
        List<String> contigStrings = assembler.solve();
        assertNotNull(assembler.kSelector.getDeBruijnGraph());
        assertNotNull(contigStrings);
        assertTrue(contigStrings.size() > 0);
        assertEquals(phix174sampleSequence.substring(0, phix174sampleSequence.length() - readStep), contigStrings.get(0));
    }

    // TODO I may instead propagate the exception KmerSizeSelector.NoSuitableKmerLengthFoundException up
    @Test
    void solveReturnsNullIfNoSuitableKmerLengthFound() {
        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler();

        int k = 15;

        List<String> reads = new ArrayList<String>();
        reads.add("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        reads.add("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

        // bypass readInParameters call
        assembler.kSelector.startingKmerSizeGuess = k;
        assembler.kSelector.setReads(reads);

        // run solve and check the output
        List<String> contigStrings = assembler.solve();
        assertNull(contigStrings);
    }

    // credit to https://www.baeldung.com/java-junit-testing-system-in
    // see also https://docs.oracle.com/en/java/javase/21/docs/api/java.base/java/lang/System.html#setIn(java.io.InputStream)
    void overrideSystemInWithString(String simulatedStdinString) {
        ByteArrayInputStream simulatedStdinIS = new ByteArrayInputStream(simulatedStdinString.getBytes());
        System.setIn(simulatedStdinIS);
    }

    @Test
    void readInParametersCanAcceptReadsViaStdin() {
        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler();
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead);
        int k = 15;
        assembler.readInParameters(k);
        assertNotNull(assembler.kSelector.reads);
        assertEquals(1, assembler.kSelector.reads.size());
        assertEquals(firstRead, GenomeSequenceEncodingUtil.getStringFromBPBytes(assembler.kSelector.reads.get(0)));
    }

    private final PrintStream standardOut = System.out;
    @AfterEach
    public void tearDown() {
        System.setOut(standardOut);
    }

    @Test
    void simulateMainWithoutArgs() {
        // credit to https://www.baeldung.com/java-testing-system-out-println
        // for this approach to check system.out output
        final ByteArrayOutputStream toCaptureOutputStream = new ByteArrayOutputStream();
        System.setOut(new PrintStream(toCaptureOutputStream));
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead);
        DeBruijnGenomeAssembler.main(new String[]{});
        String output = toCaptureOutputStream.toString();
        // should end with outputting ">CONTIG1 (<LEN>bp)" followed by the sequence
        assertEquals(">CONTIG1 (" + firstRead.length() + "bp)\n" + firstRead + "\n", output.substring(output.length() - firstRead.length() - 1 - 16));
        try {
            Files.deleteIfExists(FileSystems.getDefault().getPath("output.fasta"));
        } catch (Exception e) {}
    }

    @Test
    void simulateMainWithKArg() {
        // credit to https://www.baeldung.com/java-testing-system-out-println
        // for this approach to check system.out output
        final ByteArrayOutputStream toCaptureOutputStream = new ByteArrayOutputStream();
        System.setOut(new PrintStream(toCaptureOutputStream));
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead);
        DeBruijnGenomeAssembler.main(new String[]{"--k", "20"});
        String output = toCaptureOutputStream.toString();
        // should end with outputting ">CONTIG1 (<LEN>bp)" followed by the sequence
        assertEquals(">CONTIG1 (" + firstRead.length() + "bp)\n" + firstRead + "\n", output.substring(output.length() - firstRead.length() - 1 - 16));
        try {
            Files.deleteIfExists(FileSystems.getDefault().getPath("output.fasta"));
        } catch (Exception e) {}
    }

    @Test
    void simulateMainWithOutputFastaArg() {
        // trivial, not actually checking assembly here, just that we can control where it outputs the file
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead);
        String outputFasta = "output_assembly_attempt_in_unit_test.fasta";

        // delete the output file if it already exists before we run the code under test
        try {
            Files.deleteIfExists(FileSystems.getDefault().getPath(outputFasta));
        } catch (Exception e) {}

        // run main with an output fasta argument, which should create that file
        DeBruijnGenomeAssembler.main(new String[]{"--output_fasta", outputFasta});

        // check that the file was created with the expected contents
        File f = null;
        Scanner sc = null;
        try {
            f = new File(outputFasta);
            sc = new Scanner(f);
        } catch (FileNotFoundException e) {
            fail("Expected file was not created");
        }

        String line1 = sc.nextLine();
        String line2 = sc.nextLine();
        try {
            Files.deleteIfExists(FileSystems.getDefault().getPath(outputFasta));
        } catch (Exception e) {}        
        assertEquals(">CONTIG1 (" + firstRead.length() + "bp)", line1);
        assertEquals(firstRead, line2);
    }

    @Test
    void simulateMainWithInvalidKArgThrowsInvalidArgumentException() {
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead);
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGenomeAssembler.main(new String[]{"--k", "101"});
            });
        assertEquals("Overriding the kmer_size_starting_guess must be with an integer between 1 and 100 inclusive", exception.getMessage());
    }

    @Test
    void simulateMainWithMissingKArgNumberThrowsInvalidArgumentException() {
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead);
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGenomeAssembler.main(new String[]{"--k"});
            });
        assertEquals("'--kmer_size_starting_guess' must be followed by an integer specifying the kmer size to use, e.g. '--kmer_size_starting_guess 40'", exception.getMessage());
    }

    @Test
    void simulateMainWithInvalidArgThrowsInvalidArgumentException() {
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead);
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                    DeBruijnGenomeAssembler.main(new String[]{"--invalidarg", "1"});
            });
        assertEquals("unrecognized argument '--invalidarg', only '--output_fasta <output filename>' and '--k <positive integer>'/'--kmer_size_starting_guess <positive integer>' are supported right now", exception.getMessage());
    }

    @Test
    void simulateMainWithOutputFastaEmptyArgThrowsInvalidArgumentException() {
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead);
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                    DeBruijnGenomeAssembler.main(new String[]{"--output_fasta"});
            });
        assertEquals("'--output_fasta' must be followed by a filename, e.g. '--output_fasta carsonella_ruddii_assembly_attempt.fasta'", exception.getMessage());
    }

    @Test
    void simulateMainWithBadOutputFastaArgThrowsInvalidArgumentException() {
        String firstRead = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
        overrideSystemInWithString(firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead + "\n" + firstRead);
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGenomeAssembler.main(new String[]{"--output_fasta", "not_ending_in_fasta.txt"});
            });
        assertEquals("Invalid output FASTA filename 'not_ending_in_fasta.txt' doesn't end in '.fasta'", exception.getMessage());
    }

    // DNA sequence example for circularity corrections, might as well use a real eccDNA sequence so took the shortest one I could find in NIH NUCCORE:
    // {eccDNA 80, extrachromosomal circular DNA} [human, HeLa S3 cells, Genomic, 394 nt, but includes overlap, I use 388nt sequence from paper Fig 1 and Results section 1] (https://www.ncbi.nlm.nih.gov/nuccore/S71553.1)
    // A real extrachromosomal circular human DNA example from van Loon, N., Miller, D., & Murnane, J. P. (1994). Formation of extrachromosomal circular DNA in HeLa cells by nonhomologous recombination. Nucleic acids research, 22(13), 2447–2452. https://doi.org/10.1093/nar/22.13.2447 , not from a healthy cell, except with the one "N" (unknown) replaced with "T":
    final static String eccDNA80FromHeLaS3 = "AATTCAAAGGTGTCTCTAATCCTCTTCCACTGAACCTCTTGCTGTAACAGGCAAGGATCCTCTGCCAGGCAGCCTTGAGCCCAGACGAGGGGAAACCTGACTATACTCCGTGCATGATTTTCCCCAATAGCTTTCTTTATTGAGATAATTTCATTGTTGTTGAAGGAGAGGTAGAGGTGTTCTCAAGCCCACTGAGGGAAAGGGTTTGGCTTATTTGGGGGAAAACATGTAGTGGAGAATGTAGGGGGAGGATGTGGGAAGATGGTCCTCGAAATGCAGGCAGGCAGACCAACCAAGAAGGCCTTACTATATGCCTTGTTGAGTCTGGAGTTTATGCATGTGGGAAGTTAGGATCCAAACTTACACAATTTTCACATGATAGAAGTGG";
    // (But small eccDNAs do also occur naturally in healthy human cells see Møller, H.D., Mohiyuddin, M., Prada-Luengo, I. et al. Circular DNA elements of chromosomal origin are common in healthy human somatic tissue. Nat Commun 9, 1069 (2018). https://doi.org/10.1038/s41467-018-03369-8 .)

    // By the way, you will notice that the above sequence contains at least one repeated 10mer, e.g. ATGTGGGAAG, meaning single contig assembly gets tricky when k is <= 11 (de Bruijn graph edges are between nodes of (k-1), but for k > 11 should be easily single contig.)
    // (One can confirm this using, from the util folder, `java GenerateNoisyReadsFromWholeGenomeFasta --input_fasta eccDNA80S3_full_sequence.fasta --read_length 25 --error_probability 0 --read_coverage 10 > ../examples/eccDNA80_no_errors_25bp_reads_10x_coverage.txt` 
    // where eccDNA80S3_full_sequence.fasta is created from the sequence above, and then `java -jar ./app/build/libs/app.jar --k 11 < examples/eccDNA80_no_errors_25bp_reads_10x_coverage.txt`.)

    public static List<String> getEccDNA80Reads(boolean withCircularOverlap) {
        return getEccDNA80Reads(withCircularOverlap, 15, 30, 5, 5);
    }

    public static List<String> getEccDNA80Reads(boolean withCircularOverlap, int k, int readLen, int coverageCount, int readStep) {
        List<String> reads = new ArrayList<String>();
        for (int coverageRepeatIdx = 0; coverageRepeatIdx < coverageCount; coverageRepeatIdx++) {
            for (int i = 0; i < eccDNA80FromHeLaS3.length() - (withCircularOverlap ? 0 : readLen); i+=readStep) {
                StringBuilder read = new StringBuilder("");
                for (int j = 0; j < readLen; j++) {
                    read.append(eccDNA80FromHeLaS3.charAt((i + j) % eccDNA80FromHeLaS3.length()));
                }
                String readString = read.toString();
                reads.add(readString);
            }
        }
        return reads;
    }

    @Test
    void removeCircularEndOverlapOnContigsAndConvertToStringsCanMakeCircularCorrectionsWhenNeeded() {
        int k = 15;
        int readLen = 30;
        int coverageCount = 5;
        int readStep = 5;
        boolean withCircularOverlap = true;
        DeBruijnGraph dbg = DeBruijnGraphTest.getEccDNA80DeBruijnGraphWithOrWithoutCircularOverlapInReads(withCircularOverlap, k, readLen, coverageCount, readStep);

        DeBruijnGraphContigSolver dbgcs = new DeBruijnGraphContigSolver(dbg);
        dbgcs.solve();
        List<Contig> contigs = dbgcs.getContigs();
        Contig mainContig = contigs.get(0);
        String contigSequenceBeforeCircularityCorrection = mainContig.getContigSequence();
        assertTrue(contigSequenceBeforeCircularityCorrection.length() > eccDNA80FromHeLaS3.length());
        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler();
        List<String> contigStrings = assembler.removeCircularEndOverlapOnContigsAndConvertToStrings(contigs);
        String contigSequenceAfterCircularityCorrection = contigStrings.get(0);
        assertEquals(eccDNA80FromHeLaS3.length(), contigSequenceAfterCircularityCorrection.length());
        assertTrue(contigSequenceAfterCircularityCorrection.length() < contigSequenceBeforeCircularityCorrection.length());
        assertEquals(dbg.k - 1, contigSequenceBeforeCircularityCorrection.length() - contigSequenceAfterCircularityCorrection.length());
    }

    @Test
    void removeCircularEndOverlapOnContigsAndConvertToStringsDoesNotMakeCircularCorrectionsWhenNotNeeded() {
        int k = 15;
        int readLen = 30;
        int coverageCount = 5;
        int readStep = 5;
        boolean withCircularOverlap = false;
        DeBruijnGraph dbg = DeBruijnGraphTest.getEccDNA80DeBruijnGraphWithOrWithoutCircularOverlapInReads(withCircularOverlap, k, readLen, coverageCount, readStep);
        DeBruijnGraphContigSolver dbgcs = new DeBruijnGraphContigSolver(dbg);
        dbgcs.solve();
        List<Contig> contigs = dbgcs.getContigs();
        Contig mainContig = contigs.get(0);
        String contigSequenceBeforeCallingMethod = mainContig.getContigSequence();
        assertEquals(eccDNA80FromHeLaS3.length() - (eccDNA80FromHeLaS3.length() % readStep), contigSequenceBeforeCallingMethod.length());
        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler();
        List<String> contigStrings = assembler.removeCircularEndOverlapOnContigsAndConvertToStrings(contigs);
        String contigSequenceAfterCallingMethod = contigStrings.get(0);
        assertEquals(eccDNA80FromHeLaS3.length() - (eccDNA80FromHeLaS3.length() % readStep), contigSequenceAfterCallingMethod.length());
        assertEquals(0, contigSequenceBeforeCallingMethod.length() - contigSequenceAfterCallingMethod.length());
    }

    @Test
    void solveCanMakeCircularCorrectionsWhenNeeded() {
        int k = 15;
        int readLen = 30;
        int coverageCount = 5;
        int readStep = 5;
        boolean withCircularOverlap = true;
        List<String> reads = getEccDNA80Reads(withCircularOverlap, k, readLen, coverageCount, readStep);

        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler();

        // bypass readInParameters call
        assembler.kSelector.startingKmerSizeGuess = k;
        assembler.kSelector.setReads(reads);
        assertEquals(0, assembler.numBasesRemovedInCircularityCorrection);

        // run solve and check the output
        List<String> contigStrings = assembler.solve();
        assertNotNull(assembler.kSelector.getDeBruijnGraph());
        assertNotNull(contigStrings);
        assertTrue(contigStrings.size() > 0);
        assertEquals(eccDNA80FromHeLaS3.length(), contigStrings.get(0).length());
        assertNotEquals(eccDNA80FromHeLaS3, contigStrings.get(0)); // it is a rotation of the original string, not the original string
        assertTrue(assembler.numBasesRemovedInCircularityCorrection > 0);
        assertEquals(assembler.kSelector.k - 1, assembler.numBasesRemovedInCircularityCorrection);
    }

    @Test
    void solveDoesNotMakeCircularCorrectionsWhenNotNeeded() {
        int k = 15;
        int readLen = 30;
        int coverageCount = 5;
        int readStep = 5;
        boolean withCircularOverlap = false;
        List<String> reads = getEccDNA80Reads(withCircularOverlap, k, readLen, coverageCount, readStep);

        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler();

        // bypass readInParameters call
        assembler.kSelector.startingKmerSizeGuess = k;
        assembler.kSelector.setReads(reads);
        assertEquals(0, assembler.numBasesRemovedInCircularityCorrection);

        // run solve and check the output
        List<String> contigStrings = assembler.solve();
        assertNotNull(assembler.kSelector.getDeBruijnGraph());
        assertNotNull(contigStrings);
        assertTrue(contigStrings.size() > 0);
        assertEquals(eccDNA80FromHeLaS3.length()  - (eccDNA80FromHeLaS3.length() % readStep), contigStrings.get(0).length());
        // should not be circular correction here
        assertEquals(0, assembler.numBasesRemovedInCircularityCorrection);
    }
}
