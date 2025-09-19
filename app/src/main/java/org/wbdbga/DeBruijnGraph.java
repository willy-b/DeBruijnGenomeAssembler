/*
DeBruijnGraph for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
import java.io.Serializable;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.io.ObjectOutput;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.file.FileSystems;
import java.nio.file.Files;
/**
   The de Bruijn graph defined by (k-1)mer nodes with edges corresponding to the kmers present in the reads for some k smaller than the read length for which we have high coverage is used for genome assembly in the Eulerian method (which is linear in the number of edges), as opposed to naively trying to assemble the read overlap graph which requires finding a Hamiltonian path visiting all nodes once (instead of a Eulerian path visiting all edges once), which is NP complete.
   <p>
   Created using the de Bruijn map from the DeBruijnKmers class, solved into contigs by the DeBruijnGraphContigSolver class after cleanup in DeBruijnLowCoverageTipRemover and elsewhere.
   <p>
   Back adjacency is added to handle tip removal .
   <p>
   Based on knowledge gained from UCSD ALGS206x , can also learn the theory using the textbook from the instructors: Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach, 3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/ .
*/
public final class DeBruijnGraph implements Serializable {
    static final boolean VERBOSE = true;

    /**
       Filtered and kmer spectral corrected reads, represented as byte[] as returned by GenomeSequenceEncodingUtil.getBPBytesFromString . Unpaired data, size is number of unpaired reads. For paired data, which is not supported in this version, size is 2 * number of pairs of reads.
       @see org.wbdbga.GenomeSequenceEncodingUtil#getBPBytesFromString .
    */
    public List<byte[]> reads; // unpaired data, size is number of unpaired reads. paired data, size is 2 * number of pairs of reads.

    /**
       Length k of the kmers (a kmer is a substring of length k shorter than the read length for which we expect to have coverage in our samples of all unique values present in the genome being sequenced/assembled). In the de Bruijn graph the edges represent the kmers of the reads and the nodes represent (k-1)mers.
    */
    public int k;

    /**
       Nodes are (k-1)mers of the reads that were broken from kmers in DeBruijnKmers upstream. Each entry is a pair of integers (Integer[2]), first is index into reads array, second is the index into that read, for the first occurrence of this kmer. There should not be two entries in here which when referenced back to the actual read string have the same text. This is further encoded into byte[] by GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(int readIndex, int withinReadIndex), and decoded with GenomeSequenceEncodingUtil.getReadIndexAndWithinReadIndexFrom4ByteRepresentation(byte[] kmerBytes).
       <p>
       Values can be turned back into (k-1)mer strings using the getKMinus1MerString method.

       @see org.wbdbga.DeBruijnGraph#getKMinus1MerString
       @see org.wbdbga.GenomeSequenceEncodingUtil#get4ByteRepresentationForReadIndexAndWithinReadIndex
       @see org.wbdbga.GenomeSequenceEncodingUtil#getReadIndexAndWithinReadIndexFrom4ByteRepresentation
       @see org.wbdbga.DeBruijnKmers
    */
    public List<ValueComparedByteArray> nodeValues;

    /**
       Array of adjacency lists, index into the array using the node index in the de Bruijn graph, then can look up the outgoing directed edges in the List of Integers obtained at that index, each entry referring to a (k-1)mer suffix of this (k-1)mer which can be retrieved in nodeValues using e.g. adjacency[nodeIdx].get(suffixIdx) and transformed back into a string using the getKMinus1MerString method.
       @see org.wbdbga.DeBruijnGraph#getKMinus1MerString
    */
    public List<Integer>[] adjacency;


    /**
       Same format as the adjacency field above except each non-negative Integer entry is number of times a particular kmer edge is covered by reads.
       @see org.wbdbga.DeBruijnKmers#loadRead
       @see org.wbdbga.DeBruijnKmers#processKmers
    */
    public List<Integer>[] adjacencyCoverageCount;
    /**
       Same format as adjacency but giving the inbound edges to a node rather than the outbound, used for removing low coverage tips in the graph starting from the zero outdegree, one indegree tips and working backwards using the back adjacencies.
    */
    public List<Integer>[] backAdjacency;

    /**
       In-degree for each node by node index. Non-negative integer.
    */
    public int[] indegree;

    /**
       Out-degree for each node by node index. Non-negative integer.
    */
    public int[] outdegree;

    /**
       Number of nodes in the de Bruijn graph (number of (k-1)mers)
    */
    public int numNodes = 0;

    /**
       Number of edges in the de Bruijn graph (number of kmers)
    */
    public int numEdges = 0;

    /**
       Average read coverage is set by DeBruijnGraph and may be used as by DeBruijnLowCoverageTipRemoval to figure out which tips to remove (maybe those with 50% of the average coverage, for example).
     */
    public Double averageReadCoverage = null;

    /**
       At the time of writing, average coverage for retained edges is set by some external users of the class, not by DeBruijnGraph itself. E.g. KmerSizeSelector may prune edges and then set averageCoverageForRetainedEdges in addition to averageReadCoverage which is set here.
     */
    public Double averageCoverageForRetainedEdges = null;

    /**
       Map back from the (k-1)mer representation to its index so that we can go from its (k-1)mer to its (back) adjacency list, in degree, or out degree in constant time (O(1)).
       Note originally this was optimized (was extremely memory constrained) to use single instance of each (k-1)mer byte[] mapped by reference, but now for ease of unit testing and low fragility, any equivalent ValueComparedByteArray representing the same (k-1)mer can be used to look up its index.
    */
    Map<ValueComparedByteArray, Integer> nodeValueIndexMap = new HashMap<ValueComparedByteArray, Integer>();

    /**
       DeBruijnMap from which this de Bruijn graph was constructed, comes from DeBruijnKmers class.
       @see org.wbdbga.DeBruijnKmers
    */
    Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap;

    Integer readLen = null;

    /**
       Unpaired read data constructor, uses de Bruijn map from DeBruijnKmers.
       @param debruijnMap de Bruijn map from DeBruijnKmers.getDeBruijnMap()
       @param reads read strings encoded using GenomeSequenceEncodingUtil.getBPBytesFromString  
       @param k Length k of the kmers (a kmer is a substring of length k shorter than the read length for which we expect to have coverage in our samples of all unique values present in the genome being sequenced/assembled). In the de Bruijn graph the edges represent the kmers of the reads and the nodes represent (k-1)mers.
       @see org.wbdbga.GenomeSequenceEncodingUtil#getBPBytesFromString
       @see org.wbdbga.DeBruijnKmers#getDeBruijnMap
    */
    public DeBruijnGraph(Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap, List<byte[]> reads, int k) {
        if (debruijnMap == null || debruijnMap.keySet().size() == 0) {
            throw new IllegalArgumentException("DeBruijnGraph cannot be constructed from a null or empty debruijnMap");
        }
        this.debruijnMap = debruijnMap;
        if (reads == null || reads.size() == 0) {
            throw new IllegalArgumentException("DeBruijnGraph cannot be constructed from a null or empty reads list");
        }
        this.reads = reads;
        if (k <= 0 || k >= GenomeSequenceEncodingUtil.getDeclaredBPLength(reads.get(0))) {
            throw new IllegalArgumentException("DeBruijnGraph cannot be constructed with k values that are <= 0 or equal to the length of the read (should be substantially shorter than the read length)");
        }
        this.k = k;
        setupFromMap();
    }

    /**
       Returns node labels, not edge labels, for the de Bruijn graph, so (k-1)mers. Recall edges have length k, but nodes have length k - 1.
       @param kmerBytes kmer reference encoded using GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex (in future may support 5 byte paired format as well)
       @see org.wbdbga.GenomeSequenceEncodingUtil#get4ByteRepresentationForReadIndexAndWithinReadIndex
       @return (k-1)mer string for the node represented by kmerBytes (probably retrieved from nodeValues)
    */
    public String getKMinus1MerString(byte[] kmerBytes) {
        if (this.reads == null) {  
            throw new ReadsNullException("This DeBruijnGraph reads instance is missing, cannot decode pointers into reads into kmer strings!");
        }

        int[] kmer = GenomeSequenceEncodingUtil.decodeKmerBytesToIntArray(kmerBytes);
        int readIdx = kmer[0];
        if (readIdx < 0 || readIdx >= this.reads.size()) {
            throw new ReadIdxOutOfRangeException("DeBruijnGraph getKMinus1MerString called with invalid zero-based read index (" + readIdx + ") (# of reads: " + this.reads.size() + ")");
        }

        int offset = kmer[1];

        String read = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(readIdx));
        if (offset < 0 || offset + this.k - 1 > read.length()) {
            throw new OffsetOutOfRangeException("DeBruijnGraph getKMinus1MerString called with invalid zero-based offset (" + offset + ") (read length: " + read.length() + ")");
        }

        // -1 as for e.g. k = 5, each node has length 4 (the edge has label length 5)
        return read.substring(offset, offset + this.k - 1);
    }

    static class ReadsNullException extends RuntimeException {
        public ReadsNullException(String message) {
            super(message);
        }
    }

    static class ReadIdxOutOfRangeException extends RuntimeException {
        public ReadIdxOutOfRangeException(String message) {
            super(message);
        }
    }

    static class OffsetOutOfRangeException extends RuntimeException {
        public OffsetOutOfRangeException(String message) {
            super(message);
        }
    }

    /**
       @//deprecated
       Debugging method that can be used to generate graphviz dot language output in stdout for visualizing de Bruijn graphs.
       Marked deprecated because not currently covered by unit tests. Will have tests put in place to confirm it is working properly and then the "deprecated" tag will be removed, its name and interface may also change at this time as the author determines how best to make this available. It would be removed except it may be useful to someone as is.

       @param startNode node index to start the BFS from
    */
    /*
    @Deprecated
    public void printGraphVizBfsOrder(int startNode) {
        boolean[] visitedNodes = new boolean[adjacency.length];
        Arrays.fill(visitedNodes, false);

        StringBuilder graphviz = new StringBuilder("strict digraph {");

        Queue<Integer> nodesToVisit = new LinkedList<Integer>();
        nodesToVisit.add(startNode);

        while (!nodesToVisit.isEmpty()) {
            int nodeIdx = nodesToVisit.poll();
            if (visitedNodes[nodeIdx]) {
                continue;
            }
            visitedNodes[nodeIdx] = true;

            if (adjacency[nodeIdx] == null) { continue; }
            graphviz.append("\n    ");
            graphviz.append(getKMinus1MerString(nodeValues.get(nodeIdx).value) + " [label=\"" + getKMinus1MerString(nodeValues.get(nodeIdx).value) + "\\n(indegree: " + indegree[nodeIdx] + ", outdegree: " + outdegree[nodeIdx] + ")\"]");
            for (int j = 0; j < adjacency[nodeIdx].size(); j++) {
                int otherNodeIdx = adjacency[nodeIdx].get(j);
                graphviz.append("\n    ");
                graphviz.append(getKMinus1MerString(nodeValues.get(nodeIdx).value));
                graphviz.append(" -> ");
                graphviz.append(getKMinus1MerString(nodeValues.get(otherNodeIdx).value));
                int coverage = adjacencyCoverageCount[nodeIdx].get(j);
                graphviz.append(" [label=\"read coverage: " + coverage + "\"]");
                if (!visitedNodes[otherNodeIdx]) {
                    nodesToVisit.add(otherNodeIdx);
                }
            }
        }

        graphviz.append("\n}");
        System.out.println(graphviz);
    }
    */

    /**
       @//deprecated
       Debugging method that can be used to generate graphviz dot language output in stdout for visualizing de Bruijn graphs.
       Marked deprecated because not currently covered by unit tests. Will have tests put in place to confirm it is working properly and then the "deprecated" tag will be removed, its name and interface may also change at this time as the author determines how best to make this available. It would be removed except it may be useful to someone as is.
    */
    /*
    @Deprecated
    public void printGraphVizBfsAll() {
        boolean[] visitedNodes = new boolean[adjacency.length];
        Arrays.fill(visitedNodes, false);

        StringBuilder graphviz = new StringBuilder("strict digraph {");

        Queue<Integer> nodesToVisit = new LinkedList<Integer>();

        for (int i = 0; i < adjacency.length; i++) {
            if (!visitedNodes[i]) {
                nodesToVisit.add(i);
            }

            while (!nodesToVisit.isEmpty()) {
                int nodeIdx = nodesToVisit.poll();
                if (visitedNodes[nodeIdx]) {
                    continue;
                }
                visitedNodes[nodeIdx] = true;

                if (adjacency[nodeIdx] == null) { continue; }
                graphviz.append("\n    ");
                graphviz.append(getKMinus1MerString(nodeValues.get(nodeIdx).value) + " [label=\"" + getKMinus1MerString(nodeValues.get(nodeIdx).value) + "\\n(indegree: " + indegree[nodeIdx] + ", outdegree: " + outdegree[nodeIdx] + ")\"]");
                for (int j = 0; j < adjacency[nodeIdx].size(); j++) {
                    int otherNodeIdx = adjacency[nodeIdx].get(j);
                    graphviz.append("\n    ");
                    graphviz.append(getKMinus1MerString(nodeValues.get(nodeIdx).value));
                    graphviz.append(" -> ");
                    graphviz.append(getKMinus1MerString(nodeValues.get(otherNodeIdx).value));
                    int coverage = adjacencyCoverageCount[nodeIdx].get(j);
                    graphviz.append(" [label=\"read coverage: " + coverage + "\"]");
                    if (!visitedNodes[otherNodeIdx]) {
                        nodesToVisit.add(otherNodeIdx);
                    }
                }
            }
        }

        graphviz.append("\n}");
        System.out.println(graphviz);
    }
    */

    private void setupFromMap() {
        String firstRead = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(0));
        this.readLen = firstRead.length();

        if (VERBOSE) { System.out.println("Setting up the graph..."); }

        nodeValues = new ArrayList<ValueComparedByteArray>(debruijnMap.keySet());

        for (int lineIdx = 0; lineIdx < nodeValues.size(); lineIdx++) {
            ValueComparedByteArray nodeValue = nodeValues.get(lineIdx);
            nodeValueIndexMap.put(nodeValue, lineIdx);
        }

        int totalReadCoverage = 0;

        // first ensure we have all the node values represented (mapped to indices). not building any graph representation in this loop yet.
        for (int nodeIdx = 0; nodeIdx < nodeValues.size(); nodeIdx++) {
            if (VERBOSE && nodeIdx % 10000 == 0) {  
                System.out.println("processed " + nodeIdx + " out of " + nodeValues.size() + " nodes");
                //if (nodeIdx % 250000 == 0) { System.gc(); }
            }

            // this step can be optimized, converting the keySet into list repeatedly is an unnecessary performance hit
            if (debruijnMap.get(nodeValues.get(nodeIdx)) == null) {
                continue;
            }

            List<ValueComparedByteArray> adjacentNodes = new ArrayList<ValueComparedByteArray>(debruijnMap.get(nodeValues.get(nodeIdx)).keySet());

            for (int j = 0; j < adjacentNodes.size(); j++) {
                ValueComparedByteArray otherNodeValue = adjacentNodes.get(j);
                Integer otherNodeIdx = nodeValueIndexMap.get(otherNodeValue);
                if (otherNodeIdx == null) {
                    otherNodeIdx = nodeValues.size();
                    nodeValueIndexMap.put(otherNodeValue, otherNodeIdx);
                    nodeValues.add(otherNodeValue);
                }
            }
        }

        // now let's build the graph
        this.numNodes = nodeValues.size();
        this.numEdges = 0;
        adjacency = (List<Integer>[])new ArrayList[numNodes];
        adjacencyCoverageCount = (List<Integer>[])new ArrayList[numNodes]; // for this edge, how many times it is covered in the raw reads (ideally how many reads cover, since ideally we choose k such that kmers are unique within a read, though it is not necessary, it makes the graph less tangled) (we use this for bubble resolution where we have to choose one path as authoritative and mark the others as errors, we do this based on coverage of the edges in the path)
        backAdjacency = (List<Integer>[])new ArrayList[numNodes];

        indegree = new int[numNodes];
        Arrays.fill(indegree, 0);
        outdegree = new int[numNodes];
        Arrays.fill(outdegree, 0);

        //System.gc();
        for (int nodeIdx = 0; nodeIdx < nodeValues.size(); nodeIdx++) {
            // this step can be optimized, converting the keySet into list repeatedly is an unnecessary performance hit
            if (debruijnMap.get(nodeValues.get(nodeIdx)) == null) {
                continue;
            }

            List<ValueComparedByteArray> adjacentNodes = new ArrayList<ValueComparedByteArray>(debruijnMap.get(nodeValues.get(nodeIdx)).keySet());
            adjacency[nodeIdx] = new ArrayList<Integer>();
            adjacencyCoverageCount[nodeIdx] = new ArrayList<Integer>();

            if (backAdjacency[nodeIdx] == null) {
                backAdjacency[nodeIdx] = new ArrayList<Integer>();
            }

            for (int j = 0; j < adjacentNodes.size(); j++) {
                ValueComparedByteArray otherNodeValue = adjacentNodes.get(j);
                Integer otherNodeIdx = nodeValueIndexMap.get(otherNodeValue);
                Integer coverageCount = (int)debruijnMap.get(nodeValues.get(nodeIdx)).get(otherNodeValue);

                if (backAdjacency[otherNodeIdx] == null) {
                    backAdjacency[otherNodeIdx] = new ArrayList<Integer>();
                }

                if (adjacency[nodeIdx].contains(otherNodeIdx)) {
                    continue;
                }

                if (otherNodeIdx == null) { throw new CorruptedDeBruijnMapException("did not read map correctly"); }
                adjacency[nodeIdx].add(otherNodeIdx);
                adjacencyCoverageCount[nodeIdx].add(coverageCount);
                totalReadCoverage += coverageCount;
                if (VERBOSE && coverageCount == 0) {
                    System.out.println(getKMinus1MerString(nodeValues.get(nodeIdx).value) + " -> " + getKMinus1MerString(nodeValues.get(otherNodeIdx).value) + ": zero coverage");
                }
                backAdjacency[otherNodeIdx].add(nodeIdx);
                indegree[otherNodeIdx]++;
                //System.out.println("increased indegree for " + getKMinus1MerString(nodeValues.get(otherNodeIdx)) + " (" + otherNodeIdx + ")");  
                outdegree[nodeIdx]++;
                numEdges++;
            }

            // free up memory
            debruijnMap.remove(nodeValues.get(nodeIdx));
        }

        this.averageReadCoverage = ((double)totalReadCoverage) / ((double)numEdges);
        System.out.println("average read coverage: " + averageReadCoverage);
    }

    static class CorruptedDeBruijnMapException extends RuntimeException {
        public CorruptedDeBruijnMapException(String message) {
            super(message);
        }
    }

    static class DeBruijnGraphCorruptionException extends RuntimeException {
        public DeBruijnGraphCorruptionException(String message) {
            super(message);
        }
    }

    /**
       DBG_FILE_EXTENSION specifies the extension that de Bruijn graph files on disk are required to end in (case-insensitive) to be able to be saved, loaded, or deleted by the methods in the DeBruijnGraph class.
     */
    public static final String DBG_FILE_EXTENSION = ".dbg";

    private static boolean endsWithDbgExtension(String filename) {
        return filename != null &&
            filename.length() > DBG_FILE_EXTENSION.length() &&
            filename.substring(filename.length()-DBG_FILE_EXTENSION.length()).toLowerCase().equals(DBG_FILE_EXTENSION);
    }

    /**
       Saves a de Bruijn graph to disk as a .dbg file.
       @param dbg DeBruijnGraph instance to save to disk
       @param filename a filename that must end in ".dbg" (case insensitive)
       @throws IllegalArgumentException throws an IllegalArgumentException if the filename does not end in ".dbg"
     */
    static boolean saveToDisk(DeBruijnGraph dbg, String filename) {
        if (!endsWithDbgExtension(filename)) {
            throw new IllegalArgumentException("DeBruijnGraph.saveToDisk will only attempt to save '" + DBG_FILE_EXTENSION + "' extension files");
        }
        FileOutputStream f = null;
        ObjectOutput s = null;
        try {
            f = new FileOutputStream(filename);
            s = new ObjectOutputStream(f);
            s.writeObject(dbg);
            s.flush();
        } catch (Exception e) {
            e.printStackTrace();
            return false;
        } finally {
            if (s != null) {
                try { s.close(); } catch (Exception e) { e.printStackTrace(); }
            }
            if (f != null) {
                try { f.close(); } catch (Exception e) { e.printStackTrace(); }
            }
        }
        return true;
    }

    /**
       Removes/deletes a de Bruijn graph on disk from the disk (must be a ".dbg" file)
       @param filename a filename that must end in ".dbg" (case insensitive)
       @throws IllegalArgumentException throws an IllegalArgumentException if the filename does not end in ".dbg"
     */
    static void clearFromDisk(String filename) {
        if (!endsWithDbgExtension(filename)) {
            throw new IllegalArgumentException("DeBruijnGraph.clearFromDisk will only remove '" + DBG_FILE_EXTENSION + "' extension files");
        }
        try {
            Files.delete(FileSystems.getDefault().getPath(filename));
        } catch (Exception e) {
            e.printStackTrace();
        }        
    }

    /**
       Loads a de Bruijn graph on disk into a DeBruijnGraph object (must be a ".dbg" file)
       @param filename a filename that must end in ".dbg" (case insensitive)
       @throws IllegalArgumentException throws an IllegalArgumentException if the filename does not end in ".dbg"
     */
    static DeBruijnGraph loadFromDisk(String filename) {
        if (!endsWithDbgExtension(filename)) {
            throw new IllegalArgumentException("DeBruijnGraph.loadFromDisk will only attempt to load '" + DBG_FILE_EXTENSION + "' extension files");
        }
        FileInputStream in = null;
        ObjectInputStream s = null;
        try {
            in = new FileInputStream(filename);
            s = new ObjectInputStream(in);
            DeBruijnGraph dbg = (DeBruijnGraph)s.readObject();
            return dbg;
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (s != null) {
                try { s.close(); } catch (Exception e) { e.printStackTrace(); }
            }
            if (in != null) {
                try { in.close(); } catch (Exception e) { e.printStackTrace(); }
            }
        }
        return null;
    }
}
