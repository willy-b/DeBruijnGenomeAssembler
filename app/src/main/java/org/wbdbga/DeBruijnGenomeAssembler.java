/*
DeBruijnGenomeAssembler for DeBruijnGenomeAssembler - unpaired short read (~100bp) de Bruijn graph based genome assembler
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

import java.util.*;
import java.io.*;

/**
Takes error-prone (unpaired) reads from a genome and generates contigs in de novo assembly (no reference is used) of a genome.
<p>
Uses iterative kmer spectral alignment, filtering by coverage, de Bruijn graph spur (tip) removal / simple bubble resolution, and finally de Bruijn graph readout of non-branching paths to get the contigs.
<p>
Has the ability to attempt to dynamically pick the k value but may be best to manually select the best performing k for your problem (using the `--k` argument), as due to the wide variety of genome sizes and coverage levels the iterative k selection can be very slow.
*/
public class DeBruijnGenomeAssembler {
    static final boolean VERBOSE = true;

    KmerSizeSelector kSelector;
    int k;
    int numBasesRemovedInCircularityCorrection = 0;


    final static int MAX_CHARACTERS_PER_CONTIG_LINE = 100;
    final static boolean LINE_BREAK_CONTIGS = true;
    final static boolean WRITE_OUTPUT_FASTA = true;

    /**
       Construct the DeBruijnGenomeAssembler, after which one can call readInParameters(k) with a starting guess for k (e.g. 30) and then solve().
       Sets kSelector.
       @param takeFirstSatisfyingSolution If in doubt, use false. If false, the KmerSizeSelector will try increasing the kmer size incrementally and estimating (reference-free) whether the N50 statistic is increasing, and then take the kmer size with the best estimated N50 and use the associated assembly. If true, do not try different kmer sizes to try to determine the best one, use the first value where we are fairly confident we could assemble (across many contigs) at least 50% of the sequence (advised to only set to true if you have hand chosen the kmer size as if reproducing a known good value). True is much much faster.
    */
    public DeBruijnGenomeAssembler(boolean takeFirstSatisfyingSolution) {
        kSelector = new KmerSizeSelector(takeFirstSatisfyingSolution);
    }

    /**
       Construct the DeBruijnGenomeAssembler, after which one can call readInParameters(k) with a starting guess for k (e.g. 30) and then solve().
       Sets kSelector.
    */
    public DeBruijnGenomeAssembler() {
        this(false);
    }

    /**
       Accepts a starting guess for the k value, and then uses stdin to take in reads as ATCG strings, 1 per line, after which solve() can be called.
       @param startingKmerSizeGuess Starting guess for length k of the kmers (a kmer is a substring of length k shorter than the read length for which we expect to have coverage in our samples of all unique values present in the genome being sequenced/assembled). In the de Bruijn graph the edges will represent the kmers of the reads and the nodes represent (k-1)mers.
     */
    public void readInParameters(int startingKmerSizeGuess) {
        kSelector.readInParameters(startingKmerSizeGuess);
    }

    /**
       Method that triggers assembly, just make sure readInParameters() is called first. Finds a good value of k and performs assembly if possible, returning a list of contig strings in descending order of length.
       @return List of assembled contig strings in descending order of contig length. For smaller genomes, the entire assembly may be in the first contig. Usually short contigs less than 500bp at the end can be ignored but are included for completeness.
     */
    public List<String> solve() {
        try {
            this.k = kSelector.solve();
        } catch (KmerSizeSelector.NoSuitableKmerLengthFoundException e) {
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.out.println("Could not find a suitable k, exiting.");
            return null;
        }
        DeBruijnGraph dbg = kSelector.getDeBruijnGraph();

        DeBruijnGraphContigSolver dbgcs = new DeBruijnGraphContigSolver(dbg);
        dbgcs.solve();

        List<Contig> contigs = dbgcs.getContigs();

        List<String> contigStrings = removeCircularEndOverlapOnContigsAndConvertToStrings(contigs);

        return contigStrings;
    }

    /**
       Trim less than k length circular genome overlap from the endings of the contigs and convert to strings.
       @param contigs list of contigs ("contiguous consensus sequences") from de Bruijn assembly
       @return string ATCG sequences for contigs obtained using Contig#getContigSequence but with circular genome overlap removed using overlap alignment
       @see org.wbdbga.OverlapAlignment
     */
    public List<String> removeCircularEndOverlapOnContigsAndConvertToStrings(List<Contig> contigs) {
        int k = contigs.get(0).dbg.k;

        Contig templateContig = contigs.remove(0); // they are already ordered by length descending

        StringBuilder sequenceBuilder = new StringBuilder(templateContig.getContigSequence());

        // Handle circularity (wrap around of genome e.g. as in phi x174) by clipping off any redundant bit (applies to largest contig only - does not require that any other contigs were combined into it).
        // Theoretically we only expect up to k - 1 bases (de Bruijn node size) of circularity correction to be needed, 
        // so this only allows up to k - 1 bases to be trimmed instead of k or more for that reason to avoid a small chance of truncating further by chance. 
        // However, one will observe in the unit tests (and the examples referenced in the repo) that if you update this to allow up to k bases or more to be clipped here, only up to the theoretically expected k - 1 will be clipped (not more),
        // showing the alignment is working.
        AlignmentResult alignResult = OverlapAlignment.align(1, 10, 10, sequenceBuilder.substring(sequenceBuilder.length() - (k - 1)), sequenceBuilder.substring(0, (k - 1)));
        if (alignResult.score >= Math.min(5, (k - 1)) && alignResult.originateI >= 0) {
            if (VERBOSE) { System.out.println("Detected circularity in longest contig, removed " + ((k - 1) - alignResult.originateI) + " bases"); }
            numBasesRemovedInCircularityCorrection = (k - 1) - alignResult.originateI;
            sequenceBuilder.setLength(sequenceBuilder.length() - numBasesRemovedInCircularityCorrection);
        }

        List<String> result = new ArrayList<String>();
        result.add(sequenceBuilder.toString());

        // add in the strings for contigs not used
        for (int ci = 0; ci < contigs.size(); ci++) {
            if (contigs.size() >= 1000) {
                if (contigs.get(ci).getEdgeLength() <= 1) { continue; } // only include the trivial contigs (1 node) when it is toy example (not if 1000 contigs or more)
            }
            result.add(contigs.get(ci).getContigSequence());
        }

        return result;
    }

    /**
       Run the de Bruijn assembler!
       See the README.md for instructions, but "java -jar ./app/build/libs/app.jar &lt; examples/noisyCarsonellaRuddiiReads.txt" is a working example (after following instructions in ./util/README.md to create examples/noisyCarsonellaRuddiiReads.txt and running "./gradlew build").

       Outputs the assembled contigs to standard output after showing the assembly process and writes the latest assembly to either "output.fasta" in the current directory or whatever file is specified using the "--output_fasta" argument.

       @param args command line arguments, currently only "--output_fasta","--k", and "--kmer_size_starting_guess" are supported
     */
    public static void main(String[] args) {
        System.out.println("DeBruijnGenomeAssembler Copyright (C) 2025 William Bruns\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\nYou should have received a copy of the GNU General Public License\nalong with this program.  If not, see https://www.gnu.org/licenses/.\nThis is free software, and you are welcome to redistribute it\nunder certain conditions, see the GNU General Public License for details.\n");
        int startingKmerSizeGuess = 30;
        boolean takeFirstSatisfyingSolution = false;
        String outputFastaFilename = "output.fasta";
        for (int i = 0; i < args.length; i++) {
            switch(args[i]) {
            case "--output_fasta":
                ++i;
                if (i >= args.length) {
                    throw new IllegalArgumentException("'--output_fasta' must be followed by a filename, e.g. '--output_fasta carsonella_ruddii_assembly_attempt.fasta'");
                }
                outputFastaFilename = args[i].trim();
                break;
            case "--k":
                takeFirstSatisfyingSolution = true;
                // do NOT break here
            case "--kmer_size_starting_guess":
                ++i;
                if (i >= args.length) {
                    throw new IllegalArgumentException("'--kmer_size_starting_guess' must be followed by an integer specifying the kmer size to use, e.g. '--kmer_size_starting_guess 40'");
                }
                startingKmerSizeGuess = Integer.parseInt(args[i]);
                if (startingKmerSizeGuess < 1 || startingKmerSizeGuess > 100) {
                    throw new IllegalArgumentException("Overriding the kmer_size_starting_guess must be with an integer between 1 and 100 inclusive");
                }
                break;
            default:
                throw new IllegalArgumentException("unrecognized argument '" + args[i] + "', only '--output_fasta <output filename>' and '--k <positive integer>'/'--kmer_size_starting_guess <positive integer>' are supported right now");
            }
        }

        if (WRITE_OUTPUT_FASTA) {
            validateOutputFastaFilename(outputFastaFilename);
        }

        DeBruijnGenomeAssembler assembler = new DeBruijnGenomeAssembler(takeFirstSatisfyingSolution);
        assembler.readInParameters(startingKmerSizeGuess);
        List<String> contigs = assembler.solve();
        writeContigsToOutput(contigs, WRITE_OUTPUT_FASTA ? outputFastaFilename : null);
    }

    private static void validateOutputFastaFilename(String outputFastaFilename) {
        if (outputFastaFilename.length() <= ".fasta".length() ||
            !outputFastaFilename.substring(outputFastaFilename.length() - ".fasta".length()).toLowerCase().equals(".fasta")) {
            throw new IllegalArgumentException("Invalid output FASTA filename '" + outputFastaFilename + "' doesn't end in '.fasta'");
        }
    }

    private static void writeContigsToOutput(List<String> contigStrings, String outputFastaFilename) {
        PrintWriter out = null;
        try {
            if (outputFastaFilename != null) {
                validateOutputFastaFilename(outputFastaFilename);
                out = new PrintWriter(new BufferedWriter(new FileWriter(outputFastaFilename)));
            }
            for (int contigIdx = 0; contigIdx < contigStrings.size(); contigIdx++) {
                String contigDescriptor = ">CONTIG" + (contigIdx + 1) + (VERBOSE ? " (" + contigStrings.get(contigIdx).length() + "bp)" : ""); // fasta format (start at contig1 instead of contig0, this is arbitrary, these labels can be anything, in fact and are typically not just indices but have metadata about the machine, etc)
                String s = contigStrings.get(contigIdx);

                // optionally line breaking to make it easier to embed in PDFs / read
                if (LINE_BREAK_CONTIGS) {
                    StringBuilder sb = new StringBuilder("");
                    for (int i = 0; i <= s.length() / MAX_CHARACTERS_PER_CONTIG_LINE; i++) {
                        if (i > 0) { sb.append("\n"); }
                        sb.append(s.substring(i * MAX_CHARACTERS_PER_CONTIG_LINE, Math.min((i + 1) * MAX_CHARACTERS_PER_CONTIG_LINE, s.length())));
                    }

                    s = sb.toString();
                }

                String s2 = contigDescriptor + "\n" + s;
                System.out.println(s2);
                if (out != null) {
                    out.println(s2);
                }
            }
            if (outputFastaFilename != null) {
                out.flush();
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (out != null) {
                out.close();
            }
        }
     }
}
