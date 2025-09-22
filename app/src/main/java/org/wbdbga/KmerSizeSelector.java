/*
KmerSizeSelector for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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

/**
Class that coordinates/triggers most of the assembly (DeBruijnKmers, DeBruijnGraph, DeBruijnGraphContigSolver actually solve for the sequence): Call readInParameters(k) with your initial guess at k to load the genome reads from standard input, then solve(), then getDeBruijnGraph() to get the resulting de Bruijn graph, ready for extracting contigs using DeBruijnGraphContigSolver.

@see org.wbdbga.KmerSizeSelector#readInParameters
@see org.wbdbga.KmerSizeSelector#solve
@see org.wbdbga.KmerSizeSelector#getDeBruijnGraph
@see org.wbdbga.DeBruijnGraphContigSolver
 */
public class KmerSizeSelector {
    int readLen = -1;
    int startingKmerSizeGuess = 30;
    List<byte[]> reads = new ArrayList<byte[]>();

    boolean foundKValue = false;
    boolean takeFirstSatisfyingSolution = false;
    DeBruijnGraph bestDbgSoFar;
    
    int k = -1;
    static final boolean VERBOSE = true;
    //final boolean PRINT_DEBRUIJN_IN_GRAPHVIZ = false;
    Integer numInputLines = null;
    static final String BEST_SOLUTION_SO_FAR_TEMPORARY_FILENAME_PREFIX = "deBruijnGraphForBestSolutionSoFarTmp";

    /**
       Construct a KmerSizeSelector, which can be followed by calling readInParameters(k) with your initial guess at k to load the genome reads from standard input, then solve(), then getDeBruijnGraph() to get the resulting de Bruijn graph, ready for extracting contigs using DeBruijnGraphContigSolver.
     */
    public KmerSizeSelector(boolean takeFirstSatisfyingSolution) {
        this.takeFirstSatisfyingSolution = takeFirstSatisfyingSolution;
    }

    public KmerSizeSelector() {
        this(false);
    }
    /**
       After calling .solve() on KmerSizeSelector , one can call getDeBruijnGraph() to retrieve the de Bruijn graph which can be used with DeBruijnGraphContigSolver to extract contigs. Do not call this before solve() or it will throw a SolutionNotReadyException.
       @return de Bruijn graph for a good choice of k found by iterative search
       @throws SolutionNotReadyException if called before solve()
     */
    public DeBruijnGraph getDeBruijnGraph() throws SolutionNotReadyException {
        if (!foundKValue) {
            throw new SolutionNotReadyException("KmerSizeSelector.solve() should be called before getDeBruijnGraph can be used to retrieve the solution");
        }
        return bestDbgSoFar;
    }

    static class SolutionNotReadyException extends RuntimeException {
        public SolutionNotReadyException(String message) {
            super(message);
        }
    }

    /**
This method will read from stdin reads as 1 read per line with consistent read length expected, an example being examples/noisyCarsonellaRuddiiReads.txt (see examples/README.md for a worked example piping that into stdin). Optionally, the first line can specify the number of lines to read out of what follows. The result of calling this method, if the reads are provided via stdin, is that the reads field will be populated with GenomeSequenceEncodingUtil.getBPBytesFromString encoded reads, setting one up for the next step, which is calling the solve method, followed by getDeBruijnGraph().

@param startingKmerSizeGuess the starting kmer size k, an iterative process will try to improve on this but it helps to be in the right ballpark, go higher if you have a smaller genome or higher expected read coverage, use a lower k if you have a larger genome or low expected read coverage, e.g. for 50K reads of a 160K bp genome with 1% error, k==40 is a good start, the default is k==30.

@see org.wbdbga.KmerSizeSelector#solve
@see org.wbdbga.KmerSizeSelector#getDeBruijnGraph
     */
    public void readInParameters(int startingKmerSizeGuess) {
        this.startingKmerSizeGuess = startingKmerSizeGuess;
        Scanner sc = new Scanner(System.in);
        int numReadsSkipped = 0;
        int numLines = 0;

        while(sc.hasNextLine()) {
            if (numLines > 0 && numInputLines != null && numLines >= numInputLines) {
                break;
            }

            String line = sc.nextLine().trim(); // ignore whitespace at end
            if ("".equals(line)) { continue; } // skip blank lines

            // first line can be the number of lines
            if (readLen == -1 && numInputLines == null && line.length() > 0 && line.length() < 9 && line.matches("^[0-9]+$")) {
                numInputLines = Integer.valueOf(line);
                continue;
            }

            validateInputLineLength(line);

            numLines++;

            // process unpaired data ("READ")
            byte[] lineBytes = GenomeSequenceEncodingUtil.getBPBytesFromString(line);
            reads.add(lineBytes);
        }

        if (VERBOSE) { System.out.println("numReadsSkipped: " + numReadsSkipped); }
    }

    private void validateInputLineLength(String line) {
        // set the read length from the input length if it is not set yet 
        // (in the way appropriate for whether data is paired or not)
        if (readLen == -1) {
            readLen = line.length();
        } else if (readLen != line.length()) {
            throw new VariableReadLengthNotYetSupportedException("Variable read lengths are disabled pending adequate support. Please preprocess data to uniform length (e.g. exclude reads below a minimum length and trim longer than typical reads to a common supported length.)");
        }
    }

    static class VariableReadLengthNotYetSupportedException extends RuntimeException {
        public VariableReadLengthNotYetSupportedException(String message) {
            super(message);
        }
    }

    /**
       Set the reads (raw ATCG gene sequencing output strings, uniform length) and validate that they have the same length (variable length reads are not supported in this version).
       Populates KmerSizeSelector.reads with GenomeSequenceEncodingUtil.getBPBytesFromString encoded reads.
       @param reads list of read Strings
     */
    public void setReads(List<String> reads) {
        for (String line : reads) {
            validateInputLineLength(line);
        }

        List<byte[]> readsBytes = new ArrayList<byte[]>(reads.size());
        for (String read : reads) {
            readsBytes.add(GenomeSequenceEncodingUtil.getBPBytesFromString(read));
        }

        this.reads = readsBytes;
    }

    private DeBruijnGraph setupDeBruijnGraphForReadsAndK(List<byte[]> reads, int k) {
        DeBruijnKmers dbk = new DeBruijnKmers();

        dbk.setupFromReadStringsAndSpecifiedK(k, reads);

        dbk.processKmers();
        /*
          if (dbk.kmersNotUniqueInIndividualReads) {
          if (VERBOSE) { System.out.println("For k == " + k + ", kmers are not unique WITHIN READS. This is not a hard requirement but it suggests k should be increased. Paired reads (not supported in this version) can help with the problems associated with non-unique kmers within a read, which occur in highly repetitive genomes."); }
          //continue;
          }
        */

        Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> deBruijnMap = dbk.getDeBruijnMap();

        //System.gc();
        DeBruijnGraph dbg;
        dbk = null;
        dbg = new DeBruijnGraph(deBruijnMap, reads, k);

        deBruijnMap = null;
        //System.gc();

        System.out.println("Built DeBruijn graph (of k-1 length strings for nodes, kmer strings as edges) with " + dbg.numNodes + " nodes and " + dbg.numEdges + " edges");
        //dbg.printGraphVizBfsAll();

        /*
          if (VERBOSE && PRINT_DEBRUIJN_IN_GRAPHVIZ) {
          System.out.println("Before tip removal, graph: ");
          dbg.printGraphVizBfsAll();
          }
        */

        // remove low coverage tips
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        tp.solve();
        /*
          if (VERBOSE && PRINT_DEBRUIJN_IN_GRAPHVIZ) {
          System.out.println("After tip removal, graph: ");
          dbg.printGraphVizBfsOrder(0);
          }
        */

        // TODO compare if we did bubble cleanup before tip removal
        // remove simple bubbles here
        DeBruijnLowCoverageBubbleRemover br = new DeBruijnLowCoverageBubbleRemover(dbg);
        br.removeSimpleBubbles();

        // 2nd tip removal
        tp = new DeBruijnLowCoverageTipRemover(dbg);
        tp.solve();

        return dbg;
    }

    /**
       Iteratively adjust k and build the corresponding de Bruijn graph for that value of k, without checking against reference but using general self-consistency checks on assembly quality, return the value of k chosen and populate the bestDbgSoFar so that getDeBruijnGraph() can be called to obtain the graph for use with DeBruijnGraphContigSolver to get contigs .
       @return k Length k of the kmers (a kmer is a substring of length k shorter than the read length for which we expect to have coverage in our samples of all unique values present in the genome being sequenced/assembled). In the de Bruijn graph the edges will represent the kmers of the reads and the nodes represent (k-1)mers.
       @throws NoSuitableKmerLengthFoundException if no value of k for which an assembly can be attempted is found due to short reads, too few reads, insufficient read coverage, or other issues.
     */
    public int solve() throws NoSuitableKmerLengthFoundException {
        int INITIAL_GUESS = Math.min(this.startingKmerSizeGuess, readLen - 1); //(int)(Math.min(Math.max(2, readLen / 10) * Math.log(reads.size())/4, readLen/10)/5) * 5;

        int INCREMENT = 5; // should also scale this with number of reads
        k = -1;

        int kWithBestEstimatedN50WithoutReference = -1;
        int bestEstimatedN50WithoutReference = 0;
        
        // TODO adaptive search, jump back to best k seen so far and explore around there, been focused on actual assembler, not choosing k which is relatively easy
        for (k = INITIAL_GUESS; k < readLen; k += INCREMENT) {
            if (VERBOSE) { System.out.println("Trying k == " + k); }
            DeBruijnGraph dbg = setupDeBruijnGraphForReadsAndK(reads, k);

            /*
            if (VERBOSE && PRINT_DEBRUIJN_IN_GRAPHVIZ) {
                System.out.println("After attempted bubble resolution, graph: ");
                dbg.printGraphVizBfsOrder(0);
            }
            */

            int n50 = getN50(dbg);
            if (n50 <= 0) { 
                if (VERBOSE) { System.out.println("n50 is zero (cannot get a contig covering enough of the sequence), moving on."); }
                //if (PRINT_DEBRUIJN_IN_GRAPHVIZ) { dbg.printGraphVizBfsAll(); }
                continue;
            }

            if (n50 > bestEstimatedN50WithoutReference) {
                if (foundKValue && kWithBestEstimatedN50WithoutReference > 0) {
                    DeBruijnGraph.clearFromDisk(BEST_SOLUTION_SO_FAR_TEMPORARY_FILENAME_PREFIX + kWithBestEstimatedN50WithoutReference + DeBruijnGraph.DBG_FILE_EXTENSION); 
                }
                kWithBestEstimatedN50WithoutReference = k;
                bestEstimatedN50WithoutReference = n50;
            } else if (n50 < bestEstimatedN50WithoutReference*0.5) {
                if (VERBOSE) {
                    System.out.println("N50 dropping by more than 50% from peak, stopping search and using best k so far kWithBestEstimatedN50WithoutReference=" + kWithBestEstimatedN50WithoutReference);
                }
                break;
            }

            if (VERBOSE) {
                System.out.println("Found a candidate. k == " + k + " has n50 == " + n50 + " (NOT computed vs a reference genome but is relative to DeBruijn graph formed from reads as this is a de novo assembler).");
            }

            if (kWithBestEstimatedN50WithoutReference == k) {
                foundKValue = true;
                bestDbgSoFar = dbg;
                // if we're configured to take the first satisfactory/satisfying solution, we can stop here
                if (this.takeFirstSatisfyingSolution) {
                    // no need to save the best solution to disk in this case, we use it immediately
                    return k;
                } else {
                    DeBruijnGraph.saveToDisk(bestDbgSoFar, BEST_SOLUTION_SO_FAR_TEMPORARY_FILENAME_PREFIX + k + DeBruijnGraph.DBG_FILE_EXTENSION);
                }

                // otherwise, we need to free up memory
                // so we can continue the search
                bestDbgSoFar = null;
            }
        }

        if (foundKValue) {
            k = kWithBestEstimatedN50WithoutReference;
            if (bestDbgSoFar == null) {
                this.bestDbgSoFar = DeBruijnGraph.loadFromDisk(BEST_SOLUTION_SO_FAR_TEMPORARY_FILENAME_PREFIX + k + DeBruijnGraph.DBG_FILE_EXTENSION);
                DeBruijnGraph.clearFromDisk(BEST_SOLUTION_SO_FAR_TEMPORARY_FILENAME_PREFIX + k + DeBruijnGraph.DBG_FILE_EXTENSION);
                // if we could not reload it from disk, we can still regenerate it
                if (this.bestDbgSoFar == null) {
                    this.bestDbgSoFar = setupDeBruijnGraphForReadsAndK(reads, k);
                }
            }
            return k;
        }

        throw new NoSuitableKmerLengthFoundException("no suitable kmer length found, currently we will not attempt to assemble if we do not think we can even assemble half of the underlying sequence, there may be an option to relax this in the future");
    }

    /**
       Return the N50 computed without a reference from the de Bruijn graph itself, that is the size of the smallest contig required to include to reach 50% coverage of the edges in the de Bruijn graph if we start with the largest contigs and work towards the smallest, a measure of how fragmented the assembly is, larger is better, equal to the sequence length of the entire assembly is best (indicating single contig assembly).
       @param dbg de Bruijn graph we are computing our N50 for
       @return n50, length in bases of smallest contigs we need to start including to cover 50% of the sequence suggested by the de Bruijn graph (50% edge coverage after cleanup)
     */
    public int getN50(DeBruijnGraph dbg) {
        // read out contigs
        DeBruijnGraphContigSolver dbgcs = new DeBruijnGraphContigSolver(dbg);
        dbgcs.solve();
        // TODO update to NOT necessarily solve for contigs (it is a relatively cheap operation but still) in the n50 call as probably we are also doing that elsewhere to print them out and this is redundant

        List<Contig> contigs = dbgcs.getContigs();

        // from largest to smallest get the cumulative sequence coverage
        // stop when >= 50% of all bases, return size of smallest contig in set

        // return -1 if including all contigs does not cover.
        double cumulativeFraction = 0.0;
        int contigIdx = 0;

        // assumes contigs are sorted in size descending order, which should already be true per DeBruijnGraphContigSolver.getContigs()
        for (; contigIdx < contigs.size(); contigIdx++) {
            Contig contig = contigs.get(contigIdx);
            cumulativeFraction += contig.getCoverageFraction();

            if (cumulativeFraction >= 0.5) {
                return contig.getEdgeLength() + dbg.k - 1;
            }
        }

        return -1;
    }

    static class NoSuitableKmerLengthFoundException extends Exception {
        public NoSuitableKmerLengthFoundException(String message) {
            super(message);
        }
    }
}
