/*
AlignmentResult used by DeBruijnGenomeAssembler - unpaired short read (~100bp) de Bruijn graph based genome assembler
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
   Class returned to score alignments, such as the overlap alignment used to remove overlapping redundant content in the tips of contigs from circular sequences. See OverlapAlignmentTest#overlapAlignTestWithoutErrors or #overlapAlignTestWithErrors for an illustration of how this result works and can be used.

   @see org.wbdbga.OverlapAlignment#align
   @see org.wbdbga.DeBruijnGenomeAssembler#removeCircularEndOverlapOnContigsAndConvertToStrings 
*/
public class AlignmentResult {

    /**
       Score of the alignment, based on the number of matching positions vs the number of mismatches or insertions/deletions required to align the two sequences.
       @see org.wbdbga.OverlapAlignment#align
    */
    public int score;

    /**
       Index within the first sequence that the aligned sequence begins.
    */
    public int originateI;

    /**
       Index within the second sequence that the aligned sequence begins. For overlap alignment, this should always be zero as we are matching a suffix of the first string with a prefix of the second string, but for other types of alignments, this value could vary.
    */
    public int originateJ;

    /**
       Constructor to record an alignment result in an AlignmentResult instance.

       @param score an integer computed from the number of matches, mismatches, and indels (insertions or deletions).
       @param originateI index within the first sequence that the aligned sequence begins.
       @param originateJ index within the second sequence that the aligned sequence begins. For overlap alignment, this should always be zero as we are matching a suffix of the first string with a prefix of the second string, but for other types of alignments, this value could vary.
    */
    public AlignmentResult(int score, int originateI, int originateJ) {
        this.score = score;
        this.originateI = originateI;
        this.originateJ = originateJ;
    }
}
