/*
DeBruijnLowCoverageTipRemoverTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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

class DeBruijnLowCoverageTipRemoverTest {
    @Test
    void deBruijnLowCoverageTipRemoverConstructor() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        assertNotNull(tp);
        assertEquals(dbg, tp.dbg);
    }

    @Test
    void deBruijnLowCoverageTipRemoverConstructorThrowsIllegalArgumentException() {
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(null);
        });
        assertEquals("DeBruijnLowCoverageTipRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception.getMessage());

        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
                dbg.numNodes = 0; // override and set to 0
                DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        });
        assertEquals("DeBruijnLowCoverageTipRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception2.getMessage());

        IllegalArgumentException exception3 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
                dbg.nodeValues = new ArrayList<ValueComparedByteArray>();
                DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        });
        assertEquals("DeBruijnLowCoverageTipRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception3.getMessage());

        IllegalArgumentException exception4 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
                dbg.nodeValues = null;
                DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        });
        assertEquals("DeBruijnLowCoverageTipRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception4.getMessage());
    }

    @Test
    public void isTip() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        // Above method generates a DBG for a single 6 bp read with 2 k=5 kmers:
        // "AAAAA", "AAAAT" // 2 kmers
        // thus the (k-1)mers 4mers for the De Bruijn graph nodes are:
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        // AAAA is not a tip (index 0) (as it has the cycle ...->AAAA->AAAA->AAAA->...), AAAT is a tip (index 1)
        // now it is NOT a tip we would remove, using coverage and length to check, but this method just returns whether it literally satisfies the basic tip condition
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        assertFalse(tp.isTip(0));
        assertTrue(tp.isTip(1));
    }

    @Test
    public void removeNodeAtTipAndReturnNextNodeIdxVeryBasicCase() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        // Above method generates a DBG for a single 6 bp read with 2 k=5 kmers:
        // "AAAAA", "AAAAT" // 2 kmers
        // thus the (k-1)mers 4mers for the De Bruijn graph nodes are:
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        // AAAA is not a tip (index 0) (as it has the cycle ...->AAAA->AAAA->AAAA->...), AAAT is a tip (index 1)
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        assertEquals(2, dbg.numEdges);
        assertTrue(tp.isTip(1));
        assertEquals(1, dbg.indegree[1]);
        assertEquals(0, dbg.outdegree[1]);
        assertEquals(2, dbg.outdegree[0]);
        // note to be able to remove this, we need to lower the coverage relative to the average read coverage (being a tip is not enough, it must have low relative read coverage to be removed; otherwise we would erode the real end of the sequence as if it were a tip; recall ideally the graph would just be a chain of nodes, ending in two tips, we do not want to erode those, only the weakly covered branches off the main chain.)
        // we can do this by artificially inflating the DeBruijnGraph averageReadCoverage variable.
        dbg.averageReadCoverage *= (1.0/tp.RELATIVE_COVERAGE_THRESHOLD_TO_DROP_TIP_EDGES + 0.01);
        int nextNodeIdx = tp.removeNodeAtTipAndReturnNextNodeIdx(1);
        assertEquals(0, nextNodeIdx);
        assertEquals(1, dbg.numEdges); // an edge was removed
        assertEquals(0, dbg.indegree[1]); // changes to 0
        assertEquals(0, dbg.outdegree[1]); // stays 0
        assertEquals(1, dbg.outdegree[0]); // changes to 1
    }

    @Test
    public void removeNodeAtTipAndReturnNextNodeIdxRefusesToRemoveTipNodesWhichHaveAboveThresholdReadCoverage() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        // Above method generates a DBG for a single 6 bp read with 2 k=5 kmers:
        // "AAAAA", "AAAAT" // 2 kmers
        // thus the (k-1)mers 4mers for the De Bruijn graph nodes are:
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        // AAAA is not a tip (index 0) (as it has the cycle ...->AAAA->AAAA->AAAA->...), AAAT is a tip (index 1)
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        assertEquals(2, dbg.numEdges);
        assertTrue(tp.isTip(1));
        assertEquals(1, dbg.indegree[1]);
        assertEquals(0, dbg.outdegree[1]);
        assertEquals(2, dbg.outdegree[0]);
        // note to be able to remove this, we WOULD NEED to lower the coverage relative to the average read coverage 
        // like we did in the previous test (being a tip is not enough, it must have low relative read coverage to be removed; otherwise we would erode the real end of the sequence as if it were a tip; recall ideally the graph would just be a chain of nodes, ending in two tips, we do not want to erode those, only the weakly covered branches off the main chain)
        // we can do this by artificially inflating the DeBruijnGraph averageReadCoverage variable.
        int nextNodeIdx = tp.removeNodeAtTipAndReturnNextNodeIdx(1);
        assertEquals(-1, nextNodeIdx); // edges to node were not removed and we get -1 back
        assertEquals(2, dbg.numEdges); // number of edges is unchanged
        assertEquals(1, dbg.indegree[1]); // stays at 1
        assertEquals(0, dbg.outdegree[1]); // stays 0
        assertEquals(2, dbg.outdegree[0]); // stays at 2
    }

    @Test
    public void removeNodeAtTipAndReturnNextNodeIdxRefusesToRemoveNonTipNodesAndThrowsException() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        // Above method generates a DBG for a single 6 bp read with 2 k=5 kmers:
        // "AAAAA", "AAAAT" // 2 kmers
        // thus the (k-1)mers 4mers for the De Bruijn graph nodes are:
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        // AAAA is not a tip (index 0) (as it has the cycle ...->AAAA->AAAA->AAAA->...), AAAT is a tip (index 1)
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        assertEquals(2, dbg.numEdges);
        assertFalse(tp.isTip(0));
        assertEquals(2, dbg.outdegree[0]);
        assertEquals(1, dbg.indegree[0]);

        // note to be able to remove any edges from this graph, we need to lower the coverage relative to the average read coverage (being a tip is not enough, it must have low relative read coverage to be removed; otherwise we would erode the real end of the sequence as if it were a tip; recall ideally the graph would just be a chain of nodes, ending in two tips, we do not want to erode those, only the weakly covered branches off the main chain.)
        // we can do this by artificially inflating the DeBruijnGraph averageReadCoverage variable.
        dbg.averageReadCoverage *= (1.0/tp.RELATIVE_COVERAGE_THRESHOLD_TO_DROP_TIP_EDGES + 0.01);
        // but in this case even adjusting the read coverage will not allow edges to node 0 to be removed

        DeBruijnLowCoverageTipRemover.AttemptToRemoveNonTipNodeException exception = assertThrows(DeBruijnLowCoverageTipRemover.AttemptToRemoveNonTipNodeException.class, () -> {
                int nextNodeIdx = tp.removeNodeAtTipAndReturnNextNodeIdx(0);
        });
        assertEquals("larger of indegree/outdegree must be exactly 1 when removing a tip node", exception.getMessage());
        assertEquals(2, dbg.numEdges); // number of edges is unchanged
        assertFalse(tp.isTip(0));
        assertEquals(2, dbg.outdegree[0]); // outdegree unchanged
        assertEquals(1, dbg.indegree[0]);// indegree unchanged
    }

    @Test
    public void removeNodeAtTipAndReturnNextNodeIdxRefusesToRemoveNonTipNodesAndThrowsExceptionOneInOneOut() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTestsSimpleSingleNodeCycle();
        // only contains a 1 in 1 out node self cycle
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        assertEquals(1, dbg.numEdges);
        assertFalse(tp.isTip(0));
        assertEquals(1, dbg.outdegree[0]);
        assertEquals(1, dbg.indegree[0]);

        // note to be able to remove any edges from this graph, we need to lower the coverage relative to the average read coverage (being a tip is not enough, it must have low relative read coverage to be removed; otherwise we would erode the real end of the sequence as if it were a tip; recall ideally the graph would just be a chain of nodes, ending in two tips, we do not want to erode those, only the weakly covered branches off the main chain.)
        // we can do this by artificially inflating the DeBruijnGraph averageReadCoverage variable.
        dbg.averageReadCoverage *= (1.0/tp.RELATIVE_COVERAGE_THRESHOLD_TO_DROP_TIP_EDGES + 0.01);
        // but in this case even adjusting the read coverage will not allow edges to node 0 to be removed

        DeBruijnLowCoverageTipRemover.AttemptToRemoveNonTipNodeException exception = assertThrows(DeBruijnLowCoverageTipRemover.AttemptToRemoveNonTipNodeException.class, () -> {
                int nextNodeIdx = tp.removeNodeAtTipAndReturnNextNodeIdx(0);
        });
        assertEquals("smaller of indegree/outdegree must be exactly 0 when removing tip node", exception.getMessage());
        assertEquals(1, dbg.numEdges); // number of edges is unchanged
        assertFalse(tp.isTip(0));
        assertEquals(1, dbg.outdegree[0]); // outdegree unchanged
        assertEquals(1, dbg.indegree[0]);// indegree unchanged
    }

    @Test
    public void removeNodeAtTipAndReturnNextNodeIdxThrowsBackAdjacencyEntryCorruptedExceptionIfBackAdjacencyEntryIsCorrupted() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        // Above method generates a DBG for a single 6 bp read with 2 k=5 kmers:
        // "AAAAA", "AAAAT" // 2 kmers
        // thus the (k-1)mers 4mers for the De Bruijn graph nodes are:
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        // AAAA is not a tip (index 0) (as it has the cycle ...->AAAA->AAAA->AAAA->...), AAAT is a tip (index 1)
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        assertEquals(2, dbg.numEdges);
        assertTrue(tp.isTip(1));
        assertEquals(1, dbg.indegree[1]);
        assertEquals(0, dbg.outdegree[1]);
        assertEquals(2, dbg.outdegree[0]);

        // corrupt the BackAdjacency entry on purpose
        dbg.backAdjacency[1].clear(); // should never happen but we want a meaningful error message if we break the data structure e.g. during local development
        DeBruijnLowCoverageTipRemover.BackAdjacencyEntryCorruptedException exception = assertThrows(DeBruijnLowCoverageTipRemover.BackAdjacencyEntryCorruptedException.class, () -> {
           int nextNodeIdx = tp.removeNodeAtTipAndReturnNextNodeIdx(1);
        });
        assertEquals("indegree is one but backAdjacency does not have exactly one entry for nodeIdx==1", exception.getMessage());
    }

    @Test
    public void removeNodeAtTipAndReturnNextNodeIdxThrowsParentNodeAdjacencyEntryForChildMissingException() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        // Above method generates a DBG for a single 6 bp read with 2 k=5 kmers:
        // "AAAAA", "AAAAT" // 2 kmers
        // thus the (k-1)mers 4mers for the De Bruijn graph nodes are:
        // AAAA-AAAA, AAAA-AAAT -> 1 k-1mer with two suffixes, 2 unique k-1 mers overall
        // (an example of the case where there are fewer (k-1)-mers than kmers which is not always the case (e.g. AATAAA would not have this property, and has 2 kmers that translate to two one-to-one (k-1)-mers))
        // AAAA is not a tip (index 0) (as it has the cycle ...->AAAA->AAAA->AAAA->...), AAAT is a tip (index 1)
        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        assertEquals(2, dbg.numEdges);
        assertTrue(tp.isTip(1));
        assertEquals(1, dbg.indegree[1]);
        assertEquals(0, dbg.outdegree[1]);
        assertEquals(2, dbg.outdegree[0]);

        // manually clear the parent adjacency pointing forward to the child being removed
        // should never happen but we want a meaningful error message if we break the data structure e.g. during local development
        assertEquals(2, dbg.adjacency[0].size()); // points to itself and to node 1
        assertNotEquals(-1, dbg.adjacency[0].indexOf(Integer.valueOf(1)));
        dbg.adjacency[0].remove(dbg.adjacency[0].indexOf(Integer.valueOf(1))); // remove pointer from node 0 to node 1 (which may be at index 0 or 1 of its adjacency list)
        assertEquals(1, dbg.adjacency[0].size()); // now just points to itself

        DeBruijnLowCoverageTipRemover.ParentNodeAdjacencyEntryForChildMissingException exception = assertThrows(DeBruijnLowCoverageTipRemover.ParentNodeAdjacencyEntryForChildMissingException.class, () -> {
           int nextNodeIdx = tp.removeNodeAtTipAndReturnNextNodeIdx(1);
        });
        assertEquals("could not find reference to this node in parent adjacency", exception.getMessage());
    }

    @Test
    public void solveVeryBasicCaseIncomingTip() {
        // first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
        // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
        //        T <- on one read we introduce a point error at this position within k bases of the start, which causes a low coverage incoming pair of tips in the de Bruijn graph (can also be considered an open incoming bubble and cleaned that way)
        //        *
        DeBruijnGraph dbg = DeBruijnGraphTest.getSimpleIncomingBubbleGraphBasedOnPhiX174First70bpWith1ErrorNearStart();
        // before tip/spur removal, should be multiple contigs
        DeBruijnGraphContigSolver dbgcsBeforeTipRemoval = new DeBruijnGraphContigSolver(dbg);
        dbgcsBeforeTipRemoval.solve();
        List<Contig> contigsBeforeTipRemoval = new ArrayList<Contig>();
        for (Contig c : dbgcsBeforeTipRemoval.getContigs()) {
            if (c.getEdgeLength() > 0) {
                contigsBeforeTipRemoval.add(c);
            }
        }
        int beforeTipCleanupContigCount = contigsBeforeTipRemoval.size();
        assertTrue(beforeTipCleanupContigCount > 1);

        DeBruijnLowCoverageTipRemover tp = new DeBruijnLowCoverageTipRemover(dbg);
        tp.VERBOSE = true;
        tp.solve();

        DeBruijnGraphContigSolver dbgcsAfterTipRemoval = new DeBruijnGraphContigSolver(dbg);
        dbgcsAfterTipRemoval.solve();
        List<Contig> contigsAfterTipRemoval = new ArrayList<Contig>();
        for (Contig c : dbgcsAfterTipRemoval.getContigs()) {
            if (c.getEdgeLength() > 0) {
                contigsAfterTipRemoval.add(c);
            }
        }
        int afterTipCleanupContigCount = contigsAfterTipRemoval.size();
        assertTrue(afterTipCleanupContigCount < beforeTipCleanupContigCount);
        assertEquals(1, afterTipCleanupContigCount);
    }
}
