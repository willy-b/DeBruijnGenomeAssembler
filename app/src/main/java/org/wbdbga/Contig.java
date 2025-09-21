/*
Contig used by DeBruijnGenomeAssembler - unpaired short read (~100bp) de Bruijn graph based genome assembler
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
   Class representing a Contig but in the context of a de Bruijn graph, where the contig corresponds to a non-branching path in the de Bruijn graph (generally, a "contig" has been defined since 1980 as "a contiguous consensus sequence", see Staden R "A new computer method for the storage and manipulation of DNA gel reading data", Nucleic Acids Res, 1980 doi: 10.1093/nar/8.16.3673 .)

   Thus this class can be constructed by providing a path in the de Bruijn graph, after which getContigSequence can be called to give the actual sequence string.
   The class also provides a method getCoverageFraction() returning the fraction of edges in the de Bruijn graph covered by this contig.
   Typically there are multiple contigs, so it gives a sense of how fragmented the assembly was to see how many contigs must be considered to reach for example 50% coverage of de Bruijn edges (e.g. KmerSizeSelector#getN50 method).
   In an ideal assembly (which this assembler can come close to for 30x coverage 1% error synthetic 100bp reads of Carsonella ruddii for example), a single contig would cover the entire sequence, as the graph would be one long non-branching path.

   See chapter 3 of the textbook: Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach,3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/ for background on contigs derived using the de Bruijn graph approach.

   @see org.wbdbga.Contig#getContigSequence
   @see org.wbdbga.Contig#getCoverageFraction
   @see org.wbdbga.KmerSizeSelector#getN50

*/
public final class Contig {
    /**
       The de Bruijn graph on which this contig is defined.
       @see org.wbdbga.DeBruijnGraph
    */
    DeBruijnGraph dbg;

    /**
       Node indices (nonnegative integers) into the de Bruijn graph that should define a connected path.
    */
    List<Integer> nodeIndices;

    /**
       contigSequence is an output field populated by this class using the de Bruijn graph and indices provided.
    */
    String contigSequence;

    /**
       Construct a Contig using a de Bruijn graph and a list of node indices corresponding to a connected path in that graph (from solving a Eulerian walk usually). This populates the contigSequence field 
       (getContigSequence()) and allows requesting the coverage fraction (getCoverageFraction()).
       @see org.wbdbga.Contig#getContigSequence
       @see org.wbdbga.Contig#getCoverageFraction

       @param dbg DeBruijnGraph that the indices in the nodes argument index into
       @param nodes a list of connected nodes forming a contiguous non-branching path in the de Bruijn graph, which this class will represent and help to re-express as a String genome sequence
    */
    public Contig(DeBruijnGraph dbg, List<Integer> nodes) {
        if (dbg == null || dbg.nodeValues == null || dbg.nodeValues.size() == 0 || dbg.numNodes == 0) {
            throw new IllegalArgumentException("Contig must be constructed with a non-null, initialized, non-empty DeBruijnGraph");
        }
        this.dbg = dbg;
        if (nodes == null || nodes.size() == 0) {
            throw new IllegalArgumentException("Contig must be constructed with a non-null, non-empty nodes list");
        }
        for (Integer nodeIdx : nodes) {
            if (nodeIdx < 0 || nodeIdx >= this.dbg.nodeValues.size()) {
                throw new IllegalArgumentException("Contig cannot be constructed from a node index list that contains indices that are out of bounds, < 0 or greater than the largest index of the provided DeBruijnGraph nodeValues list. Found " + nodeIdx + " which is not 0 <= node index < " + this.dbg.nodeValues.size());
            }
        }
        // TODO check whether node indices are from a connected path using dbg adjacency and throw an exception if not, rather than trusting the upstream to construct this using valid data.
        this.nodeIndices = nodes;

        // generate contig sequence
        StringBuilder sb = new StringBuilder("");

        for (int i = 0; i < this.nodeIndices.size(); i++) {
            int nodeIdx = this.nodeIndices.get(i);

            ValueComparedByteArray nodeValue = this.dbg.nodeValues.get(nodeIdx);
            String kmerString = dbg.getKMinus1MerString(nodeValue.value);

            if (i == 0) {
                sb.append(kmerString);
            } else {
                sb.append(kmerString.charAt(kmerString.length() - 1));
            }
        }

        this.contigSequence = sb.toString();
    }

    /**
       Returns the fraction of edges in the de Bruijn graph covered by this contig.
       @see org.wbdbga.KmerSizeSelector#getN50 .
       @return fraction of edges in the de Bruijn graph covered by the edges in this contig
    */
    public double getCoverageFraction() {
        int numDbgEdges = this.dbg.numEdges;
        int numEdgesInContig = Math.max(this.nodeIndices.size() - 1, 0);
        return ((double)numEdgesInContig)/((double)numDbgEdges);
    }

    /**
       Returns the contig length in graph edges.
       To get equivalent sequence length add (dbg.k - 1) for the prefix of the first DeBruijnGraph node, or use getContigSequence and call .length() on that String.
       @return length in edges of the contig (non negative integer)
    */
    public int getEdgeLength() {
        return this.nodeIndices.size() - 1;
    }

    /**
       Return the ATCG genetic sequence for the contig.
       @return String containing ATCG genetic sequence of contig
    */
    public String getContigSequence() {
        return this.contigSequence;
    }
}
