/*
GenerateNoisyReadsFromWholeGenomeFasta - generate sample noisy unpaired short reads to practice the associated DeBruijnGenomeAssembler unpaired short read (~100bp) de Bruijn graph based genome assembler on
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

import java.lang.Math;
import java.io.File;
import java.lang.StringBuffer;
import java.util.Scanner;
import java.util.regex.Pattern;
import java.io.FileNotFoundException;
import java.util.Random;

// Simple utility to give us synthetic genome sequencing data as input to test our assembler.
// Here we sample reads (e.g. by default of 100 base pairs with 1% errors at 30x coverage) to give us something to assemble.
// (We can also use real genome sequencing reads from NCBI.)
public class GenerateNoisyReadsFromWholeGenomeFasta {
    Random random;
    public GenerateNoisyReadsFromWholeGenomeFasta(long randomSeed) {
        this.random = new Random(randomSeed);
    }
    void generateSyntheticReads(String sequence, int readLength, int targetCoverage, double errorProbability) {
        int sequenceLength = sequence.length();
        int numReads = (sequenceLength*targetCoverage)/readLength;
        for (int readIdx = 0; readIdx < numReads; readIdx++) {
            // extract a read
            int readStart = (int)Math.round(this.random.nextDouble()*(sequenceLength-1));
            // Note this generates synthetic sequencer reads with circular overlap on purpose
            // as sequencer reads from our example datasets should have this property
            // (e.g. definitely phi x174, carsonella ruddii, and E. coli, and I think B. brevis as well)
            // and to show off that our assembler can correct it.
            // A flag may be added to allow generating data without circular overlap in the near future.
            int readEnd = readStart + readLength; // can go beyond sequence length, we will sample modulo for circular sample
            StringBuffer readSb = new StringBuffer();
            for (int i = readStart; i < readEnd; i++) {
                int idx = i % sequenceLength;
                char c = sequence.charAt(idx);
                if (this.random.nextDouble() <= errorProbability) {
                    char c2 = c;
                    while (c2 == c) {
                        c2 = "ATCG".charAt((int)Math.round(3*this.random.nextDouble()));
                    }
                    c = c2;
                }
                readSb.append(c);
            }
            // print the synthetic reads out to stdout
            System.out.println(readSb.toString());
        }
    }

    static String readSequenceFromWholeGenomeSingleEntryFASTA(String inputFasta) throws FileNotFoundException {
        // read in a whole genome sequence from the provided FASTA file
        File f = new File(inputFasta);
        Scanner sc = new Scanner(f);
        StringBuffer sb = new StringBuffer();

        // crude, take 2nd line of single entry full genome FASTA
        // we do check that it is properly formatted below
        int l = 0;
        while(sc.hasNextLine()) {
            String line = sc.nextLine();
            sb.append(line);
            if (l == 0) { 
                if (line.charAt(0) != '>') {
                    throw new IllegalArgumentException("Please ensure the first line begins with '>' and that this is the only heading in the FASTA file (whole genome single fasta entry expected)");
                }
                sb.append("\n");
            } else {
                if (line.contains(">")) {
                    throw new IllegalArgumentException("Please ensure only the first line begins with '>'; that this is the only heading in the FASTA file (whole genome single fasta entry expected)");
                }
            }
            l++;
        }
        String s = sb.toString();
        String[] parts = s.split("\n");
        if (parts.length != 2) {
            throw new IllegalArgumentException("Unexpected format; expect whole genome fasta with single entry");
        }
        String sequence = s.split("\n")[1];
        if (Pattern.matches(".*[^ATCG].*", sequence)) {
            throw new IllegalArgumentException("Input sequence should only contain A,T,C,G");
        }
        return sequence;
    }

    public static void main(String args[]) throws Exception {
        int targetCoverage = 30;
        int readLength = 100;
        double errorProbability = 0.01;
        String inputFasta = null; //"carsonella_ruddii_complete_genome_AP009180.1.fasta";
        long randomSeed = 1L;
        final String HELP_STRING = """
HELP
Usage:
```
javac GenerateNoisyReadsFromWholeGenomeFasta.java
java GenerateNoisyReadsFromWholeGenomeFasta --input_fasta carsonella_ruddii_complete_genome_AP009180.1.fasta > ../examples/noisyCarsonellaRuddiiReads.txt
```
The noisy reads text file can then be used as input to DeBruijnGenomeAssembler in the parent folder to check ability to reassemble the genome from the noisy reads (fragments).

Supported arguments:
  --input_fasta: specify an input FASTA file containing a whole genome to generate synthetic reads from
 (e.g. carsonella_ruddii_complete_genome_AP009180.1.fasta from https://www.ncbi.nlm.nih.gov/nuccore/AP009180.1?report=fasta&log$=seqview&format=text , 
  or run 
`wget 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=116235183&extrafeat=null&conwithfeat=on&hide-cdd=on' -O carsonella_ruddii_complete_genome_AP009180.1.fasta` 
  to create it)
  --error_probability: probability of single base pair error in synthetic reads generated (double; 0 <= error_probability <= 1; defaults to 0.01)
    e.g. "--error_probability 0.01" for 1% errors
  --read_length: number of base pairs long each synthetic read should be (integer 25 <= read length <= 200; defaults to 100 bp)
  --read_coverage: how many synthetic reads to generate based on target coverage == (number of reads)*(read length)/(genome length), average number of times we expect each base to be sampled (integer: 1 <= coverage; defaults to 30x coverage)
  --random_seed: some long, using the same values for the same input fasta gives a reproducible random sampling
""";

        // read in arguments from command line
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            switch (arg) {
            case "--help":
                System.out.println(HELP_STRING);
                return;
            case "--error_probability":
                ++i;
                if (i >= args.length) {
                    System.err.println("--error_probability should be followed by a double specifying the error probability, e.g. '--error_probability 0.01");
                    return;
                }
                errorProbability = Double.valueOf(args[i]);
                if (errorProbability < 0) {
                    throw new IllegalArgumentException("Error probability must be >= 0");
                }
                if (errorProbability > 1) {
                    throw new IllegalArgumentException("Error probability must be <= 1");
                }
                break;
            case "--read_length":
                ++i;
                if (i >= args.length) {
                    System.err.println("--read_length should be followed by a int specifying the read length, e.g. '--read_length 100");
                    return;
                }
                readLength = Integer.valueOf(args[i]);
                if (readLength < 25) {
                    throw new IllegalArgumentException("Read length must be >= 25");
                }
                if (readLength > 200) {
                    throw new IllegalArgumentException("Read length must be <= 200");
                }
                break;
            case "--read_coverage":
                ++i;
                if (i >= args.length) {
                    System.err.println("--read_coverage should be followed by a int specifying the read coverage, e.g. '--read_coverage 30");
                    return;
                }
                targetCoverage = Integer.valueOf(args[i]);
                if (targetCoverage < 1) {
                    throw new IllegalArgumentException("Read coverage must be >= 1");
                }
                break;
            case "--input_fasta":
                ++i;
                if (i >= args.length) {
                    System.err.println("--input_fasta should be followed by a string specifying a fasta filename, e.g. '--input_fasta carsonella_ruddii_complete_genome_AP009180.1.fasta', where that file is obtainable from https://www.ncbi.nlm.nih.gov/nuccore/AP009180.1?report=fasta&log$=seqview&format=text");
                    return;
                }
                inputFasta = args[i];
                // TODO add validation of path here
                break;
            case "--random_seed":
                ++i;
                if (i >= args.length) {
                    System.err.println("--random_seed should be followed by a number (long) specifying a random seed, e.g. '--random_seed 1'");
                    return;
                }
                randomSeed = Long.parseLong(args[i]);
                break;
            }
        }

        if (inputFasta == null) {
            System.out.println(HELP_STRING);
            return;
        }

        String sequence = readSequenceFromWholeGenomeSingleEntryFASTA(inputFasta);

        // generate synthetic reads
        GenerateNoisyReadsFromWholeGenomeFasta generator = new GenerateNoisyReadsFromWholeGenomeFasta(randomSeed);
        generator.generateSyntheticReads(sequence, readLength, targetCoverage, errorProbability);
    }
}
