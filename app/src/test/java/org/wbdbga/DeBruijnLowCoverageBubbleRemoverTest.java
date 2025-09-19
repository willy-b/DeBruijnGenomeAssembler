/*
DeBruijnLowCoverageBubbleRemoverTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

class DeBruijnLowCoverageBubbleRemoverTest {

    @Test
    void deBruijnLowCoverageBubbleRemoverConstructor() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        DeBruijnLowCoverageBubbleRemover br = new DeBruijnLowCoverageBubbleRemover(dbg);
        br.VERBOSE = true;
        assertNotNull(br);
        assertEquals(dbg, br.dbg);
    }

    @Test
    void deBruijnLowCoverageBubbleRemoverConstructorThrowsIllegalArgumentException() {
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnLowCoverageBubbleRemover br = new DeBruijnLowCoverageBubbleRemover(null);
        });
        assertEquals("DeBruijnLowCoverageBubbleRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception.getMessage());

        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
                dbg.numNodes = 0; // override and set to 0
                DeBruijnLowCoverageBubbleRemover br = new DeBruijnLowCoverageBubbleRemover(dbg);
        });
        assertEquals("DeBruijnLowCoverageBubbleRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception2.getMessage());

        IllegalArgumentException exception3 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
                dbg.nodeValues = new ArrayList<ValueComparedByteArray>();
                DeBruijnLowCoverageBubbleRemover br = new DeBruijnLowCoverageBubbleRemover(dbg);
        });
        assertEquals("DeBruijnLowCoverageBubbleRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception3.getMessage());

        IllegalArgumentException exception4 = assertThrows(IllegalArgumentException.class, () -> {
                DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
                dbg.nodeValues = null;
                DeBruijnLowCoverageBubbleRemover br = new DeBruijnLowCoverageBubbleRemover(dbg);
        });
        assertEquals("DeBruijnLowCoverageBubbleRemover must be constructed with a non-null, initialized, non-empty DeBruijnGraph", exception4.getMessage());
    }

    @Test
    void removeSimpleBubblesShouldCleanupLowCoverageBubbleFromPointErrorInRead() {
        // first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
        // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
        //                                     A <- on one read we introduce a point error at this position, which causes a low coverage bubble in the de Bruijn graph
        //                                     *
        DeBruijnGraph dbg = DeBruijnGraphTest.getSimpleBubbleGraphBasedOnPhiX174First70bpWith1Error();

        // before bubble removal, should be multiple contigs
        DeBruijnGraphContigSolver dbgcsBeforeBubbleRemoval = new DeBruijnGraphContigSolver(dbg);
        dbgcsBeforeBubbleRemoval.solve();
        List<Contig> contigsBeforeBubbleRemoval = new ArrayList<Contig>();
        for (Contig c : dbgcsBeforeBubbleRemoval.getContigs()) {
            if (c.getEdgeLength() > 0) {
                contigsBeforeBubbleRemoval.add(c);
            }
        }
        int beforeBubbleCleanupContigCount = contigsBeforeBubbleRemoval.size();
        assertTrue(beforeBubbleCleanupContigCount > 1);
        assertEquals(3, beforeBubbleCleanupContigCount);

        // after bubble removal, should be one contig which is not an isolated node
        DeBruijnLowCoverageBubbleRemover br = new DeBruijnLowCoverageBubbleRemover(dbg);
        br.VERBOSE = true;
        br.removeSimpleBubbles();
        DeBruijnGraphContigSolver dbgcsAfterBubbleRemoval = new DeBruijnGraphContigSolver(dbg);
        dbgcsAfterBubbleRemoval.solve();
        // bubble removal will not remove the cleaned up nodes from the DBG entirely, but will disconnect them into 0 edge contigs
        List<Contig> contigsAfterBubbleRemoval = new ArrayList<Contig>();
        for (Contig c : dbgcsAfterBubbleRemoval.getContigs()) {
            if (c.getEdgeLength() > 0) {
                contigsAfterBubbleRemoval.add(c);
            }
        }

        int afterBubbleCleanupContigCount = contigsAfterBubbleRemoval.size();
        assertTrue(afterBubbleCleanupContigCount < beforeBubbleCleanupContigCount);
        /*
        first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
        GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
                                            A <- on one read we introduce a point error at this position, which causes a low coverage bubble in the de Bruijn graph
                                            *
        Before, 3 contigs with more than 0 edges:
        GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACT <- up to point error but stopping right before
        CAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATT <- path which does not contain point error
        CAGAAGTTAACACTATCGGATATTTCTG <- path which contains point error
                      A
                      * error at "CACTA" (last A in CACTA)

        After bubble cleanup, 1 contig with more than 0 edges, representing almost entire sequence:
        GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATT
        */
        assertEquals(1, afterBubbleCleanupContigCount);
        assertEquals("GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATT", contigsAfterBubbleRemoval.get(0).getContigSequence());

        // the counts above are important, the specific fragments below for before bubble cleanup are not
        // but we assert specific sequences so that if behavior changes subtly the test breaks to signal that to force developer to acknowledge they are changing something
        assertEquals("CAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATT", contigsBeforeBubbleRemoval.get(0).getContigSequence());
        assertEquals("GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACT", contigsBeforeBubbleRemoval.get(1).getContigSequence());
        assertEquals("CAGAAGTTAACACTATCGGATATTTCTG", contigsBeforeBubbleRemoval.get(2).getContigSequence());
    }

    @Test
    void removeSimpleBubblesShouldCleanupLowCoverageIncomingBubbleFromPointErrorInRead() {
        // first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
        // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
        //        T <- on one read we introduce a point error at this position within k bases of the start, which causes a low coverage incoming bubble in the de Bruijn graph
        //        *
        DeBruijnGraph dbg = DeBruijnGraphTest.getSimpleIncomingBubbleGraphBasedOnPhiX174First70bpWith1ErrorNearStart();

        // before bubble removal, should be multiple contigs
        DeBruijnGraphContigSolver dbgcsBeforeBubbleRemoval = new DeBruijnGraphContigSolver(dbg);
        dbgcsBeforeBubbleRemoval.solve();
        List<Contig> contigsBeforeBubbleRemoval = new ArrayList<Contig>();
        for (Contig c : dbgcsBeforeBubbleRemoval.getContigs()) {
            if (c.getEdgeLength() > 0) {
                contigsBeforeBubbleRemoval.add(c);
            }
        }
        int beforeBubbleCleanupContigCount = contigsBeforeBubbleRemoval.size();
        assertTrue(beforeBubbleCleanupContigCount > 1);

        // after bubble removal, should be one contig which is not an isolated node
        DeBruijnLowCoverageBubbleRemover br = new DeBruijnLowCoverageBubbleRemover(dbg);
        br.VERBOSE = true;
        br.removeSimpleBubbles();
        DeBruijnGraphContigSolver dbgcsAfterBubbleRemoval = new DeBruijnGraphContigSolver(dbg);
        dbgcsAfterBubbleRemoval.solve();
        // bubble removal will not remove the cleaned up nodes from the DBG entirely, but will disconnect them into 0 edge contigs
        List<Contig> contigsAfterBubbleRemoval = new ArrayList<Contig>();
        for (Contig c : dbgcsAfterBubbleRemoval.getContigs()) {
            if (c.getEdgeLength() > 0) {
                contigsAfterBubbleRemoval.add(c);
            }
        }

        int afterBubbleCleanupContigCount = contigsAfterBubbleRemoval.size();
        assertTrue(afterBubbleCleanupContigCount < beforeBubbleCleanupContigCount);
        assertTrue(contigsAfterBubbleRemoval.get(0).getContigSequence().length() > contigsBeforeBubbleRemoval.get(0).getContigSequence().length());
        assertEquals("GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATT", contigsAfterBubbleRemoval.get(0).getContigSequence());
        /*
          first 70bp of NC_001422.1 Escherichia phage phiX174 (virus that infects E. coli) ( https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta )
          GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
                 T <- on one read we introduced a point error at this position within k bases of the start, which causes a low coverage incoming bubble in the de Bruijn graph
                 *

          Before incoming bubble removal, longest contig is missing the start where the error occurs on one read
                  TCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATT
          GAGTTTTATCGCTTCCATGACG <-- beginning without error path
          GAGTTTTTTCGCTTCCATGACG <-- beginning with error path

          After incoming bubble removal, sequence before error included at beginning of longest contig:
          GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATT
          GAGTTTTTTCGCTTCCATGAC <-- still a separate short contig showing the path with the error, which would be ignored due to low coverage
        */
        assertEquals("TCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATT", contigsBeforeBubbleRemoval.get(0).getContigSequence());
    }
}
