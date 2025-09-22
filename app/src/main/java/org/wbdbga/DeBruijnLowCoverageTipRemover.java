/*
DeBruijnLowCoverageTipRemover for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
Erosion of low coverage tips (also called spurs) in the de Bruijn graph (variously referred to as tip erosion, tip removal, spur removal, graph shaving, or other ways) appear to be part of mainstream genome assemblers like Euler, SPAdes, and Velvet, and we use it here also as it is easy and improves assembly quality. Inspired by the UCSD ALGS207x course which also assigned a problem (without providing code) involving removing such features from a graph.

Miller, J. R., Koren, S., Sutton, G. (2010). Assembly algorithms for next-generation sequencing data. Genomics, 95(6), 315-327. https://www.sciencedirect.com/science/article/pii/S0888754310000492 discusses spurs and their removal, see their Figure 3a (this algorithm is not from that paper but they give good context to the problem.)
*/
public final class DeBruijnLowCoverageTipRemover {
    boolean VERBOSE = false;
    DeBruijnGraph dbg;
    static final double RELATIVE_COVERAGE_THRESHOLD_TO_DROP_TIP_EDGES = 0.50;

    /**
       Constructs a DeBruijnLowCoverageTipRemover from an initialized DeBruijnGraph. Make sure to pass a non-empty initialized DeBruijnGraph or you will get an IllegalArgumentException. Call solve() to perform the low coverage tip removal.
       @param dbg initialized, non-null, non-empty DeBrujnGraph to perform low coverage tip removal on
     */
    public DeBruijnLowCoverageTipRemover(DeBruijnGraph dbg) {
        if (dbg == null || dbg.nodeValues == null || dbg.nodeValues.size() == 0 || dbg.numNodes == 0) {
            throw new IllegalArgumentException("DeBruijnLowCoverageTipRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph");
        }
        this.dbg = dbg;
    }

    // TODO consider to return numTipEdgesRemoved as a returned int as well.
    /**
       Remove low coverage tips from the de Bruijn graph as they probably reflect point errors in the reads or other artifacts and will produce tiny meaningless contigs. You can find the number of tip edges removed afterwards in numTipEdgesRemoved .
     */
    public void solve() {
        // while we can
        // - find the indegree > 0, outdegree == 0 nodes, put them in a list
        List<Integer> tipNodes = null;

        int numRemovedThisCycle = 0;
        tipNodes = new ArrayList<Integer>();

        for (int i = 0; i < dbg.numNodes; i++) {
            if (isTip(i)) {
                tipNodes.add(i);
            }
        }

        System.out.println(tipNodes.size() + " tip nodes found");
        //if (tipNodes.size() == 0) { break; }

        // - for each in that list, walk backwards, removing parent edge repeatedly until indegree > 0 and outdegree > 0
        for (int j = 0; j < tipNodes.size(); j++) {
            int nodeIdx = tipNodes.remove(j--);
            while (isTip(nodeIdx)) {
                nodeIdx = removeNodeAtTipAndReturnNextNodeIdx(nodeIdx);
                if (nodeIdx != -1) {
                    numTipEdgesRemoved++;
                    numRemovedThisCycle++;
                }
            }
            //if (VERBOSE && numTipEdgesRemoved % 1000 == 0) { System.out.println("Removed " + numTipEdgesRemoved + " tip edges"); }
        }

        System.out.println(numRemovedThisCycle + " tip nodes removed");

        if (VERBOSE) {
            System.out.println("Removed " + numTipEdgesRemoved + " tip edges");
        }
    }

    /**
       solve() populates numTipEdgesRemoved with the number of tip edges removed during the cleanup process.
    */
    int numTipEdgesRemoved = 0;

    /**
       Returns whether the node at the specified node index is a tip (not necessarily low read coverage): that is zero outdegree but nonzero indegree or outdegree of one with zero indegree.

       @param nodeIdx the index of the node in the de Bruijn graph to check whether it is a tip.
       @return true if this node is a tip (not necessarily low read coverage), false otherwise.
     */
    public boolean isTip(int nodeIdx) {
        if (nodeIdx == -1) { return false; }
        return (dbg.indegree[nodeIdx] == 1 && dbg.outdegree[nodeIdx] == 0) ||
            (dbg.outdegree[nodeIdx] == 1 && dbg.indegree[nodeIdx] == 0);
    }

    /**
       Removes a tip node and returns the upstream node, call this repeatedly after checking isTip to remove a low coverage branch/tip from the de Bruijn graph (returning -1 means no further nodes should be removed in this way due to reaching the coverage threshold on the edges.)
       @param nodeIdx the index of the node to remove (must be a tip node or will generate exception!)
       @return the index of the upstream node which could be the parent or the child depending on what type of tip this is (returns -1 if tip removal was aborted, e.g. due to exceeding coverage threshold at the next node)
     */
    public int removeNodeAtTipAndReturnNextNodeIdx(int nodeIdx) {
        if (Math.max(dbg.indegree[nodeIdx], dbg.outdegree[nodeIdx]) != 1) {
            throw new AttemptToRemoveNonTipNodeException("larger of indegree/outdegree must be exactly 1 when removing a tip node");
        }

        if (Math.min(dbg.indegree[nodeIdx], dbg.outdegree[nodeIdx]) != 0) {
            throw new AttemptToRemoveNonTipNodeException("smaller of indegree/outdegree must be exactly 0 when removing tip node");
        }

        if (dbg.outdegree[nodeIdx] == 0) {
            if (dbg.backAdjacency[nodeIdx] == null || ((List<Integer>)dbg.backAdjacency[nodeIdx]).size() != 1) {
                throw new BackAdjacencyEntryCorruptedException("indegree is one but backAdjacency does not have exactly one entry for nodeIdx==" + nodeIdx);
            }

            int parentIdx = dbg.backAdjacency[nodeIdx].get(0);
            int parentAdjacencyToThisIdx = dbg.adjacency[parentIdx].indexOf(nodeIdx);
            if (parentAdjacencyToThisIdx == -1) {
                throw new ParentNodeAdjacencyEntryForChildMissingException("could not find reference to this node in parent adjacency");
            }

            double averageCoverageToUse = dbg.averageReadCoverage;
            if (dbg.adjacencyCoverageCount[parentIdx].get(parentAdjacencyToThisIdx) > averageCoverageToUse * RELATIVE_COVERAGE_THRESHOLD_TO_DROP_TIP_EDGES) {
                return -1;
            }

            dbg.outdegree[parentIdx]--;
            dbg.indegree[nodeIdx]--;
            dbg.backAdjacency[nodeIdx].remove(0);
            dbg.adjacency[parentIdx].remove(parentAdjacencyToThisIdx);
            dbg.adjacencyCoverageCount[parentIdx].remove(parentAdjacencyToThisIdx);
            dbg.numEdges--;
            if (dbg.adjacency[parentIdx].size() != dbg.adjacencyCoverageCount[parentIdx].size()) { throw new DeBruijnGraph.DeBruijnGraphCorruptionException("adjacency related data structures sizes do not match which should never happen"); }
            return parentIdx;
        } else {
            int childIdx = dbg.adjacency[nodeIdx].get(0);
            double averageCoverageToUse = dbg.averageReadCoverage;
            if (dbg.adjacencyCoverageCount[nodeIdx].get(0) > averageCoverageToUse * RELATIVE_COVERAGE_THRESHOLD_TO_DROP_TIP_EDGES) {
                return -1;
            }
            int childBackAdjacencyToThisIdx = dbg.backAdjacency[childIdx].indexOf(nodeIdx);
            dbg.backAdjacency[childIdx].remove(childBackAdjacencyToThisIdx);
            dbg.indegree[childIdx]--;
            dbg.outdegree[nodeIdx]--;
            dbg.adjacency[nodeIdx].remove(0);
            dbg.adjacencyCoverageCount[nodeIdx].remove(0);
            dbg.numEdges--;
            if (dbg.adjacency[nodeIdx].size() != dbg.adjacencyCoverageCount[nodeIdx].size()) { throw new DeBruijnGraph.DeBruijnGraphCorruptionException("adjacency related data structures sizes do not match which should never happen"); }
            return childIdx;
        }
    }

    static class AttemptToRemoveNonTipNodeException extends RuntimeException {
        public AttemptToRemoveNonTipNodeException(String message) {
            super(message);
        }
    }

    static class BackAdjacencyEntryCorruptedException extends RuntimeException {
        public BackAdjacencyEntryCorruptedException(String message) {
            super(message);
        }
    }

    static class ParentNodeAdjacencyEntryForChildMissingException extends RuntimeException {
        public ParentNodeAdjacencyEntryForChildMissingException(String message) {
            super(message);
        }
    }
}
