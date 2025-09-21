/*
DeBruijnLowCoverageBubbleRemover for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
       Remove simple bubbles in the de Bruijn graph due to point errors (erroneous base calls) in sequencer reads. There are more sophisticated techniques for bubble cleanup and some were even explored in UCSD ALGS2xx based on iterated-BFS / max flow that are NOT used here due to uneven performance on the various genomes I checked.

       Miller, J. R., Koren, S., Sutton, G. (2010). Assembly algorithms for next-generation sequencing data. Genomics, 95(6), 315-327. https://www.sciencedirect.com/science/article/pii/S0888754310000492 discusses bubbles along with other errors, see their Fig 3b (the algorithm used here is not from that paper but they discuss the problem well.)
 */
public final class DeBruijnLowCoverageBubbleRemover {
    DeBruijnGraph dbg;
    static final double RELATIVE_COVERAGE_THRESHOLD_TO_DROP_BUBBLE_EDGES = 0.80;
    boolean VERBOSE = false;

    public DeBruijnLowCoverageBubbleRemover(DeBruijnGraph dbg) {
        if (dbg == null || dbg.nodeValues == null || dbg.nodeValues.size() == 0 || dbg.numNodes == 0) {
            throw new IllegalArgumentException("DeBruijnLowCoverageBubbleRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph");
        }
        this.dbg = dbg;
    }
    
    /**
       Remove simple bubbles in the de Bruijn graph due to point errors (erroneous base calls) in sequencer reads. Note that the nodes are NOT removed from the de Bruijn graph so the erroneous base calls still show in the contigs list, just at the very bottom as tiny short single node (zero edge) or few edge contigs (resulting in contigs of approximately k bases), that would be ignored in most subsequent pipelines, like QUAST which only considers contigs of 500bp or more.

       There are more sophisticated techniques for bubble cleanup and some were even explored in UCSD ALGS2xx based on iterated-BFS / max flow that are NOT used here due to uneven performance on the various genomes I checked without going even further in the implementation (simple implementation of more sophisticated technique was not sufficient).

       Miller, J. R., Koren, S., Sutton, G. (2010). Assembly algorithms for next-generation sequencing data. Genomics, 95(6), 315-327. https://www.sciencedirect.com/science/article/pii/S0888754310000492 discusses bubbles along with other errors, see their Fig 3b (this algorithm is not from that paper but they discuss the problem well.)
     */
    public void removeSimpleBubbles() {
        Queue<Integer> disconnectedNodesToRemoveFurtherEdgesFrom = new LinkedList<Integer>();

        int numEdgesRetained = 0;
        int totalReadCoverageForRetainedEdges = 0;
        double averageCoverageForRetainedEdges = 0.0;

        int numEdgesDeleted = 0;
        int totalReadCoverageForDeletedEdges = 0;
        double averageCoverageForDeletedEdges = 0.0;

        // outgoing edge bubbles
        for (int nodeIdx = 0; nodeIdx < dbg.adjacency.length; nodeIdx++) {
            if (dbg.adjacency[nodeIdx] == null) { continue; }
            if (dbg.adjacency[nodeIdx].size() != dbg.adjacencyCoverageCount[nodeIdx].size()) { throw new DeBruijnGraph.DeBruijnGraphCorruptionException("adjacency related data structures sizes do not match which should never happen"); }
            if (dbg.outdegree[nodeIdx] <= 1) { continue; }
            if (dbg.adjacency[nodeIdx].size() <= 1) { continue; }

            int maxEdgeCoverage = Integer.MIN_VALUE;
            for (int edgeIdx = 0; edgeIdx < dbg.adjacency[nodeIdx].size(); edgeIdx++) {
                int coverage = dbg.adjacencyCoverageCount[nodeIdx].get(edgeIdx);
                if (coverage > maxEdgeCoverage) { maxEdgeCoverage = coverage; }
            }

            List<Integer> edgesToRemove = new ArrayList<Integer>();
            for (int edgeIdx = 0; edgeIdx < dbg.adjacency[nodeIdx].size(); edgeIdx++) {
                int coverage = dbg.adjacencyCoverageCount[nodeIdx].get(edgeIdx);
                if (coverage < maxEdgeCoverage * RELATIVE_COVERAGE_THRESHOLD_TO_DROP_BUBBLE_EDGES) {
                    edgesToRemove.add(edgeIdx);
                    numEdgesDeleted++;
                    totalReadCoverageForDeletedEdges += coverage;
                } else {
                    numEdgesRetained++;
                    totalReadCoverageForRetainedEdges += coverage;
                }
            }

            int offset = 0;
            for (Integer edgeIdx : edgesToRemove) {
                // recall: we are removing outgoing edge bubbles in this block. do not confuse yourself with the block below (also using variable named "edgesToRemove" removing incoming edge bubbles.
                int otherNodeIdx = dbg.adjacency[nodeIdx].get(edgeIdx - offset);
                int backAdjacencyIdx = dbg.backAdjacency[otherNodeIdx].indexOf(nodeIdx);
                if (backAdjacencyIdx != -1) {
                    dbg.backAdjacency[otherNodeIdx].remove(backAdjacencyIdx);
                    dbg.indegree[otherNodeIdx]--;
                    if (dbg.indegree[otherNodeIdx] == 0) {
                        disconnectedNodesToRemoveFurtherEdgesFrom.add(otherNodeIdx); // for outgoing edge cleanup we need to enqueue downstream disconnected nodes for cleanup. this step won't happen for incoming edge bubbles.
                    }
                }

                dbg.adjacency[nodeIdx].remove(edgeIdx - offset);
                dbg.adjacencyCoverageCount[nodeIdx].remove(edgeIdx - offset);

                dbg.outdegree[nodeIdx]--;
                dbg.numEdges--;
                offset++;
            }

            if (dbg.adjacency[nodeIdx].size() != dbg.adjacencyCoverageCount[nodeIdx].size()) { throw new DeBruijnGraph.DeBruijnGraphCorruptionException("adjacency related data structures sizes do not match which should never happen"); }
        }

        // remove incoming edge bubbles
        for (int nodeIdx = 0; nodeIdx < dbg.adjacency.length; nodeIdx++) {
            if (dbg.adjacency[nodeIdx] == null) { continue; }
            if (dbg.adjacency[nodeIdx].size() != dbg.adjacencyCoverageCount[nodeIdx].size()) { throw new DeBruijnGraph.DeBruijnGraphCorruptionException("adjacency related data structures sizes do not match which should never happen"); }
            if (dbg.indegree[nodeIdx] <= 1) { continue; }
            if (dbg.backAdjacency[nodeIdx].size() <= 1) { continue; }

            // get highest coverage incoming edge
            int maxEdgeCoverage = Integer.MIN_VALUE;
            for (int edgeIdx = 0; edgeIdx < dbg.backAdjacency[nodeIdx].size(); edgeIdx++) {
                int otherNodeIdx = dbg.backAdjacency[nodeIdx].get(edgeIdx);
                int otherEdgeIdx = dbg.adjacency[otherNodeIdx].indexOf(nodeIdx);
                if (otherEdgeIdx == -1) { continue; }

                int coverage = dbg.adjacencyCoverageCount[otherNodeIdx].get(otherEdgeIdx);
                if (coverage > maxEdgeCoverage) { maxEdgeCoverage = coverage; }
            }

            List<Integer> edgesToRemove = new ArrayList<Integer>();
            for (int edgeIdx = 0; edgeIdx < dbg.backAdjacency[nodeIdx].size(); edgeIdx++) {
                int otherNodeIdx = dbg.backAdjacency[nodeIdx].get(edgeIdx);
                int otherEdgeIdx = dbg.adjacency[otherNodeIdx].indexOf(nodeIdx);
                if (otherEdgeIdx == -1) { continue; }

                int coverage = dbg.adjacencyCoverageCount[otherNodeIdx].get(otherEdgeIdx);
                if (coverage < maxEdgeCoverage * RELATIVE_COVERAGE_THRESHOLD_TO_DROP_BUBBLE_EDGES) {
                    edgesToRemove.add(edgeIdx);
                    numEdgesDeleted++;
                    totalReadCoverageForDeletedEdges += coverage;
                } else {
                    numEdgesRetained++;
                    totalReadCoverageForRetainedEdges += coverage;
                }
            }

            int offset = 0;
            for (Integer edgeIdx : edgesToRemove) {
                // recall: we are removing incoming edge bubbles in this block. do not confuse yourself with the block above (also using variable named "edgesToRemove" removing outgoing edge bubbles.

                int otherNodeIdx = dbg.backAdjacency[nodeIdx].get(edgeIdx - offset); // the back adjacency is on this node in this case, reverse of earlier in this method which was for forward edges (this is removing an incoming edge)
                int otherEdgeIdx = dbg.adjacency[otherNodeIdx].indexOf(nodeIdx);
                if (otherEdgeIdx == -1) { continue; }

                dbg.adjacency[otherNodeIdx].remove(otherEdgeIdx);
                dbg.adjacencyCoverageCount[otherNodeIdx].remove(otherEdgeIdx);
                dbg.outdegree[otherNodeIdx]--; // it is sending the edge to us in this case

                dbg.backAdjacency[nodeIdx].remove(edgeIdx - offset);
                dbg.indegree[nodeIdx]--;
                dbg.numEdges--;
                offset++;
            }
            if (dbg.adjacency[nodeIdx].size() != dbg.adjacencyCoverageCount[nodeIdx].size()) { throw new DeBruijnGraph.DeBruijnGraphCorruptionException("adjacency related data structures sizes do not match which should never happen"); }
        }

        averageCoverageForRetainedEdges = ((double)totalReadCoverageForRetainedEdges) / ((double)numEdgesRetained);
        if (VERBOSE) {
            if (numEdgesRetained == 0 && numEdgesDeleted == 0) {
                System.out.println("No bubbles were modified.");
            } else {
                System.out.println("averageCoverageForRetainedEdges == " + averageCoverageForRetainedEdges + ", retained " + numEdgesRetained + " (the ones considered, not all edges in graph which is " + dbg.numEdges + " before removing the edges in the disconnected components as well which usually counts for a lot)");
            }
        }

        averageCoverageForDeletedEdges = ((double)totalReadCoverageForDeletedEdges) / ((double)numEdgesDeleted);
        if (VERBOSE) { 
            if (numEdgesRetained != 0 || numEdgesDeleted != 0) {
                System.out.println("averageCoverageForDeletedEdges == " + averageCoverageForDeletedEdges + ", deleted " + numEdgesDeleted + " edges"); 
            }
        }

        while (!disconnectedNodesToRemoveFurtherEdgesFrom.isEmpty()) {
            Integer nodeIdx = disconnectedNodesToRemoveFurtherEdgesFrom.poll(); // get a node we just disconnected

            // if this node has no outgoing edges, we're done
            // if this node has more than one outgoing edge, that is a more complex situation not to be resolved here
            if (dbg.adjacency[nodeIdx] == null || dbg.adjacency[nodeIdx].size() == 0 || dbg.adjacency[nodeIdx].size() > 1) {
                continue;
            }

            if (dbg.indegree[nodeIdx] != 0) { continue; } // we are removing disconnected nodes here. (this is enforced on adding to the queue so this line is redundant but being careful.)

            int coverage = dbg.adjacencyCoverageCount[nodeIdx].get(0);
            double averageCoverageToUse = (dbg.averageCoverageForRetainedEdges != null && !Double.isNaN(dbg.averageCoverageForRetainedEdges)) ? averageCoverageForRetainedEdges : dbg.averageReadCoverage;
            if (coverage > (averageCoverageToUse * RELATIVE_COVERAGE_THRESHOLD_TO_DROP_BUBBLE_EDGES)) {
                numEdgesRetained++;
                totalReadCoverageForRetainedEdges += coverage;
                continue; // this will begin a disconnected contig. in this case errors broke up our graph but we don't want to delete this fragment even though we disconnected it by cutting an error edge.
            }

            numEdgesDeleted++;
            totalReadCoverageForDeletedEdges += coverage;

            int otherNodeIdx = dbg.adjacency[nodeIdx].get(0);
            int backAdjacencyIdx = dbg.backAdjacency[otherNodeIdx].indexOf(nodeIdx);
            if (backAdjacencyIdx != -1) {
                dbg.backAdjacency[otherNodeIdx].remove(backAdjacencyIdx);
                dbg.indegree[otherNodeIdx]--;
                if (dbg.indegree[otherNodeIdx] == 0) {
                    disconnectedNodesToRemoveFurtherEdgesFrom.add(otherNodeIdx);
                }
            }
            dbg.adjacency[nodeIdx].remove(0);
            dbg.adjacencyCoverageCount[nodeIdx].remove(0);
            dbg.outdegree[nodeIdx]--;
            dbg.numEdges--;
            if (dbg.adjacency[nodeIdx].size() != dbg.adjacencyCoverageCount[nodeIdx].size()) { throw new DeBruijnGraph.DeBruijnGraphCorruptionException("adjacency related data structures sizes do not match which should never happen"); }
        }

        if (numEdgesRetained > 0) {
            averageCoverageForRetainedEdges = ((double)totalReadCoverageForRetainedEdges) / ((double)numEdgesRetained);
            dbg.averageCoverageForRetainedEdges = averageCoverageForRetainedEdges;

            if (VERBOSE) { System.out.println("after cleaning up disconnected components, averageCoverageForRetainedEdges == " + averageCoverageForRetainedEdges + ", retained " + numEdgesRetained + " (the ones considered, not all edges in graph which is " + dbg.numEdges + ")"); }
        }

        if (numEdgesDeleted > 0) {
            averageCoverageForDeletedEdges = ((double)totalReadCoverageForDeletedEdges) / ((double)numEdgesDeleted);
            if (VERBOSE) { System.out.println("after cleaning up disconnected components, averageCoverageForDeletedEdges == " + averageCoverageForDeletedEdges + ", deleted " + numEdgesDeleted + " edges"); }
        }
    }
}
