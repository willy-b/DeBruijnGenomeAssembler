/*
ContigTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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

import java.util.List;
import java.util.ArrayList;

class ContigTest {
    @Test
    public void contigConstructor() {
        /*
        // "AAAAAT"; // "AAAAA", "AAAAT" // 2 kmers
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        */
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        List<Integer> nodes = new ArrayList<Integer>();
        assertTrue(dbg.nodeValues.size() >= 2);
        nodes.add(0);
        nodes.add(1);
        Contig contig = new Contig(dbg, nodes);
        assertNotNull(contig);
        assertEquals("AAAAT", contig.contigSequence);
        assertEquals("AAAAT", contig.getContigSequence());
    }

    @Test
    public void contigConstructorThrowsIllegalArgumentExceptionIfDeBruijnGraphIsNotReady() {
        List<Integer> nodes = new ArrayList<Integer>();
        nodes.add(0);
        nodes.add(1);
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                new Contig(null, nodes);
            });
        assertEquals("Contig must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception.getMessage());

        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        dbg.nodeValues = new ArrayList<ValueComparedByteArray>();
        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                new Contig(dbg, nodes);
            });
        assertEquals("Contig must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception2.getMessage());

        DeBruijnGraph dbg2 = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        dbg2.nodeValues = null;
        IllegalArgumentException exception3 = assertThrows(IllegalArgumentException.class, () -> {
                new Contig(dbg2, nodes);
            });
        assertEquals("Contig must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception3.getMessage());

        DeBruijnGraph dbg3 = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        dbg3.numNodes = 0;
        IllegalArgumentException exception4 = assertThrows(IllegalArgumentException.class, () -> {
                new Contig(dbg3, nodes);
            });
        assertEquals("Contig must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception4.getMessage());
    }

    @Test
    public void contigConstructorThrowsIllegalArgumentExceptionNodesListIsEmptyNullOrContainsOutOfBoundsIndices() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                new Contig(dbg, null);
            });
        assertEquals("Contig must be constructed with a non-null, non-empty nodes list", exception.getMessage());

        List<Integer> nodes = new ArrayList<Integer>();
        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                new Contig(dbg, nodes);
            });
        assertEquals("Contig must be constructed with a non-null, non-empty nodes list", exception2.getMessage());

        nodes.add(-1);
        IllegalArgumentException exception3 = assertThrows(IllegalArgumentException.class, () -> {
                new Contig(dbg, nodes);
            });
        assertEquals("Contig cannot be constructed from a node index list that contains indices that are out of bounds, < 0 or greater than the largest index of the provided DeBruijnGraph nodeValues list. Found -1 which is not 0 <= node index < 2", exception3.getMessage());

        nodes.remove(0);
        nodes.add(2);
        IllegalArgumentException exception4 = assertThrows(IllegalArgumentException.class, () -> {
                new Contig(dbg, nodes);
            });
        assertEquals("Contig cannot be constructed from a node index list that contains indices that are out of bounds, < 0 or greater than the largest index of the provided DeBruijnGraph nodeValues list. Found 2 which is not 0 <= node index < 2", exception4.getMessage());
    }

    @Test
    public void getCoverageFraction() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        List<Integer> nodes = new ArrayList<Integer>();
        nodes.add(0);
        nodes.add(1);
        Contig contig = new Contig(dbg, nodes);
        double numContigEdges = (contig.nodeIndices.size() - 1) + 0.0;
        assertEquals(numContigEdges/((double)dbg.numEdges), contig.getCoverageFraction());
    }

    // length in dbg edges
    @Test
    public void getEdgeLength() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        List<Integer> nodes = new ArrayList<Integer>();
        assertTrue(dbg.nodeValues.size() >= 2);
        nodes.add(0);
        nodes.add(1);
        Contig contig = new Contig(dbg, nodes);
        assertNotNull(contig);
        assertEquals(2 - 1, contig.getEdgeLength());
        assertEquals(contig.getEdgeLength() + (dbg.k - 1), contig.getContigSequence().length());
    }
}
