/*
DeBruijnGraphContigSolver for DeBruijnGenomeAssembler - unpaired short read (~100bp) de Bruijn graph based genome assembler
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
Find the non-branching paths in the de Bruijn graph and enumerate them as contigs ("contiguous consensus sequence" defined since Staden R "A new computer method for the storage and manipulation of DNA gel reading data", Nucleic Acids Res, 1980 doi: 10.1093/nar/8.16.3673).
<p>
Description of this algorithm (not this code but the approach) is also available in chapter 3 of the textbook: Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach,3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/ .
*/
public final class DeBruijnGraphContigSolver {
    List<List<Integer>> nonBranchingPaths = null;
    boolean[] visitedNode = null;
    List<Contig> contigs = null;
    DeBruijnGraph dbg;

    /**
       Construct a de Bruijn graph contig solver for a given de Bruijn graph (can then call solve() to get the contigs).

       @param dbg de Bruijn graph to extract contigs from
     */
    public DeBruijnGraphContigSolver(DeBruijnGraph dbg) {
        if (dbg == null || dbg.nodeValues == null || dbg.nodeValues.size() == 0 || dbg.numNodes == 0) {
            throw new IllegalArgumentException("DeBruijnGraphContigSolver must be constructed with a non-null, non-empty DeBruijnGraph");
        }
        this.dbg = dbg;
    }

    /**
       Compute contigs ("contiguous consensus sequence") from the de Bruijn graph as non-branching paths
     */
    public void solve() {
        //System.gc();
        // find the start nodes for non-branching paths
        List<Integer> nonOneInOneOutNodes = new ArrayList<Integer>();
        for (int nodeIdx = 0; nodeIdx < dbg.numNodes; nodeIdx++) {
            if (dbg.indegree[nodeIdx] != 1 || dbg.outdegree[nodeIdx] != 1) {
                //System.out.println("Starting a non-branching path from: " + dbg.getKMinus1MerString(dbg.nodeValues.get(nodeIdx)) + " (indegree: " + dbg.indegree[nodeIdx] + ", outdegree: " + dbg.outdegree[nodeIdx] + ")");
                nonOneInOneOutNodes.add(nodeIdx);
            }
        }

        visitedNode = new boolean[dbg.numNodes];
        Arrays.fill(visitedNode, false);
        nonBranchingPaths = new ArrayList<List<Integer>>();
        List<Integer> currentPath = null;
        for (int i = 0; i < nonOneInOneOutNodes.size(); i++) {
            int nodeIdx = nonOneInOneOutNodes.get(i);
            // visited will only be checked when handling the all one-in-one-out cycles
            // otherwise one-in-one-out paths can only be reached from one non-one-in-one-out node,
            // which we will visit only once in this loop
            visitedNode[nodeIdx] = true;

            if (dbg.adjacency[nodeIdx] == null || dbg.adjacency[nodeIdx].size() == 0) {
                // no outgoing edges
                if (dbg.indegree[nodeIdx] > 0) {
                    // we'll get here from an incoming edge
                    continue;
                } else {
                    // isolated node
                    // we could exclude these single node contigs since k << 500bp which is QUAST cutoff for example
                    // but I stay consistent with some sources, e.g. ALGS2xx and also include isolated nodes (user can filter them out and they show at the bottom of the list)
                    // we failed to assemble these fragments into larger contigs, but we can still account for them in our listing as unassembled minimal contigs
                    currentPath = new ArrayList<Integer>();
                    currentPath.add(nodeIdx);
                    nonBranchingPaths.add(currentPath);
                    continue;
                }
            }

            for (int j = 0; j < dbg.adjacency[nodeIdx].size(); j++) {
                currentPath = new ArrayList<Integer>();
                currentPath.add(nodeIdx);                
                Integer otherNodeIdx = dbg.adjacency[nodeIdx].get(j);
                currentPath.add(otherNodeIdx);
                visitedNode[otherNodeIdx] = true;

                // note we do not check whether we already visited here since these are one-in-one-out nodes so are excluded from nonOneInOneOutNodes
                // they cannot be visited multiple times as they are one-in-one-out so have only one parent non-one-in-one-out node which we visit once
                while (otherNodeIdx != null && dbg.outdegree[otherNodeIdx] == 1 && dbg.indegree[otherNodeIdx] == 1) {
                    if (otherNodeIdx != null && dbg.adjacency[otherNodeIdx] != null && dbg.adjacency[otherNodeIdx].size() > 0) {
                        otherNodeIdx = dbg.adjacency[otherNodeIdx].get(0);
                        currentPath.add(otherNodeIdx);
                        visitedNode[otherNodeIdx] = true;
                    } else {
                        otherNodeIdx = null;
                    }
                }

                // we intentionally allow 1 edge paths from non-one-in-one-out node to non-one-in-one-out-node
                // per definition of a non-branching path as beginning at one non-one-in-one-out node and ending at another with all interior nodes being one-in-one-out
                // that said, the short contigs resulting from single edge non branching paths (length k+1) will be excluded by QUAST as <500bp contigs are dropped
                // it would be fine to exclude them here (require length > 1 to add currentPath to nonBranchingPaths),
                // but I stay consistent with Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach,3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/ , and keep them.
                nonBranchingPaths.add(currentPath);
            }
        }

        // now just get the isolated cycles (all one-in-one-out nodes)
        // here we make use of visitedNode
        for (int nodeIdx = 0; nodeIdx < this.dbg.numNodes; nodeIdx++) {
            if (visitedNode[nodeIdx]) {
                continue;
            }

            // we already visited all non-one-in-one-out nodes earlier, but to be clear
            // we can assert it explicitly here
            if (dbg.outdegree[nodeIdx] != 1 || dbg.indegree[nodeIdx] != 1) {
                continue;
            }

            currentPath = new ArrayList<Integer>();
            // as a starting condition to enter the while loop, we set `otherNodeIdx = nodeIdx;`
            // and to add the first node to currentPath.
            Integer otherNodeIdx = nodeIdx;
            while (otherNodeIdx != null && (!visitedNode[otherNodeIdx] || otherNodeIdx.equals(nodeIdx))) {
                currentPath.add(otherNodeIdx);
                visitedNode[otherNodeIdx] = true;
                if (dbg.adjacency[otherNodeIdx] != null && dbg.adjacency[otherNodeIdx].size() > 0 
                    // it is already guaranteed that this is a one-in-one-out node, since we visited the rest in nonOneInOneOutNodes earlier
                    // but we assert it to be careful
                    && dbg.outdegree[otherNodeIdx] == 1 && dbg.indegree[otherNodeIdx] == 1
                    // we do NOT proceed in self-loops
                    && !dbg.adjacency[otherNodeIdx].get(0).equals(otherNodeIdx)) {
                    // node we only get the first child since any unvisited nodes at this point MUST be one-in-one-out nodes 
                    // since we visited the rest in nonOneInOneOutNodes earlier AND we assert it above
                    otherNodeIdx = dbg.adjacency[otherNodeIdx].get(0);
                } else {
                    otherNodeIdx = null;
                }
            }

            nonBranchingPaths.add(currentPath);
        }

        this.contigs = new ArrayList<Contig>();

        for (int pathIdx = 0; pathIdx < nonBranchingPaths.size(); pathIdx++) {
            int startNodeIdx = nonBranchingPaths.get(pathIdx).get(0);
            int endNodeIdx = nonBranchingPaths.get(pathIdx).get(nonBranchingPaths.get(pathIdx).size() - 1);

            int startIndegree = this.dbg.indegree[startNodeIdx] - (this.dbg.outdegree[startNodeIdx] - 1);
            int endOutdegree = this.dbg.outdegree[endNodeIdx] - (this.dbg.indegree[endNodeIdx] - 1);

            Contig contig = new Contig(this.dbg, nonBranchingPaths.get(pathIdx));
            this.contigs.add(contig);

            if (startIndegree == endOutdegree && startIndegree > 1) {
                for (int i = 0; i < startIndegree - 1; i++) {
                    contig = new Contig(this.dbg, nonBranchingPaths.get(pathIdx));
                    this.contigs.add(contig);
                }
            }
        }

        // sort descending order of contig length
        Collections.sort(this.contigs, new Comparator<Contig>() {
                public int compare(Contig a, Contig b) {
                    return -1 * Integer.valueOf(a.getEdgeLength()).compareTo(Integer.valueOf(b.getEdgeLength()));
                }
            });
    }

    /**
       getContigs returns the contigs ("contiguous consensus sequence") computed from the de Bruijn graph (non branchiing paths) using the solve() method, in descending order of their length (longest contig first).
       @return list of contigs, longest sorted first, that cover the edges of the de Bruijn graph and represent assembled sequence from the reads the de Bruijn graph was built from
       @throws ContigsNotAvailableUntilSolveCalledException if getContigs() is called before solve() which actually computes the contigs
       @see org.wbdbga.Contig
     */
    public List<Contig> getContigs() {
         if (this.contigs == null) {
             throw new ContigsNotAvailableUntilSolveCalledException("DeBruijnGraphContigSolver.getContigs() called before contigs have been solved for; please call .solve() first");
         }
        return new ArrayList(this.contigs);
    }

    static class ContigsNotAvailableUntilSolveCalledException extends RuntimeException {
         public ContigsNotAvailableUntilSolveCalledException(String message) {
             super(message);
         }
    }
}
