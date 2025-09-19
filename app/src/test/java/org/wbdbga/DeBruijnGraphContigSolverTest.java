/*
DeBruijnGraphContigSolverTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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

class DeBruijnGraphContigSolverTest {
    @Test
    public void deBruijnGraphContigSolverConstructor() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        DeBruijnGraphContigSolver dbgcs = new DeBruijnGraphContigSolver(dbg);
        assertNotNull(dbgcs);
        assertEquals(dbg, dbgcs.dbg);
    }

    @Test
    public void deBruijnGraphContigSolverConstructorThrowsIllegalArgumentExceptionForNullOrEmptyDeBruijnGraph() {
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                new DeBruijnGraphContigSolver(null);
        });
        assertEquals("DeBruijnGraphContigSolver must be constructed with a non-null, non-empty DeBruijnGraph", exception.getMessage());

        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        dbg.nodeValues = null;
        IllegalArgumentException exception2 = assertThrows(IllegalArgumentException.class, () -> {
                new DeBruijnGraphContigSolver(dbg);
        });
        assertEquals("DeBruijnGraphContigSolver must be constructed with a non-null, non-empty DeBruijnGraph", exception2.getMessage());

        DeBruijnGraph dbg2 = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        dbg2.nodeValues = new ArrayList<ValueComparedByteArray>();
        IllegalArgumentException exception3 = assertThrows(IllegalArgumentException.class, () -> {
                new DeBruijnGraphContigSolver(dbg2);
        });
        assertEquals("DeBruijnGraphContigSolver must be constructed with a non-null, non-empty DeBruijnGraph", exception3.getMessage());

        DeBruijnGraph dbg3 = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        dbg3.numNodes = 0;
        IllegalArgumentException exception4 = assertThrows(IllegalArgumentException.class, () -> {
                new DeBruijnGraphContigSolver(dbg3);
        });
        assertEquals("DeBruijnGraphContigSolver must be constructed with a non-null, non-empty DeBruijnGraph", exception4.getMessage());
    }

    @Test
    public void getContigs() {
         DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        DeBruijnGraphContigSolver dbgcs = new DeBruijnGraphContigSolver(dbg);
        assertNotNull(dbgcs);
         // artificially set the contigs on dbgcs for this unit test
         dbgcs.contigs = new ArrayList<Contig>();
        assertEquals(dbgcs.contigs, dbgcs.getContigs());
    }

    @Test
    public void getContigsThrowsContigsNotAvailableUntilSolveCalledExceptionIfCalledPrematurely() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        DeBruijnGraphContigSolver dbgcs = new DeBruijnGraphContigSolver(dbg);
        assertNotNull(dbgcs);
        DeBruijnGraphContigSolver.ContigsNotAvailableUntilSolveCalledException exception = assertThrows(DeBruijnGraphContigSolver.ContigsNotAvailableUntilSolveCalledException.class, () -> {
                  dbgcs.getContigs();
        });
        assertEquals("DeBruijnGraphContigSolver.getContigs() called before contigs have been solved for; please call .solve() first", exception.getMessage());
    }

    @Test
    public void solveExtremelySimpleCaseJustConfirmNoExceptionsThrownAndAllVariablesSet() {
        DeBruijnGraph dbg = DeBruijnGraphTest.setupDeBruijnGraphForTests();
        DeBruijnGraphContigSolver dbgcs = new DeBruijnGraphContigSolver(dbg);
        assertNotNull(dbgcs);
        assertNull(dbgcs.nonBranchingPaths);
        assertNull(dbgcs.visitedNode);
        assertNull(dbgcs.contigs);
        dbgcs.solve();
        assertNotNull(dbgcs.nonBranchingPaths);
        assertNotNull(dbgcs.visitedNode);
        assertNotNull(dbgcs.contigs);
        // not meaningful contigs for this toy example (too small input data)
        System.out.println("dbgcs.contigs: " + dbgcs.contigs);
    }
}
