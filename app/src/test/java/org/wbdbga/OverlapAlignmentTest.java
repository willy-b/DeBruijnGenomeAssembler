/*
OverlapAlignmentTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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

/*
In this project, we do NOT use overlap alignment directly for assembly, but use the more efficient de Bruijn graph based assembly, so why is the OverlapAlignment dynamic programming based aligner here?
This overlap aligner is used by this assembler in a final cleanup/polishing step to find and remove circular overlap (less than or equal to k in length) at the end of contigs for sequences from circular genomes (it is not expensive to perform overlap alignment just on the up to k bases at the ends of each contig to see if they circularly align with themselves to trim off any such circular content.) So it is only used in the specific case where t overlaps with the end of s by < k (<< s.length(), t.length()) characters. However, it is good to test this to handle more cases.
*/

class OverlapAlignmentTest {
    @Test
    public void overlapAlignTestWithoutErrors() {
         // In overlap alignment, the end of s overlaps with the beginning of t, so we start at the beginning of t, but start anywhere in s.
         String s = "Chapter 5 of Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach, 3rd ed. Active Learning Publishers, 2018";
         String t = "an active learning approach, 3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/";
         String expectedOverlap = "an active learning approach, 3rd ed. Active Learning Publishers, 2018";
         int matchReward = 1;
         int mismatchPenalty = 10;
         int indelPenalty = 10;
         AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, s, t);
         assertEquals(expectedOverlap.length()*matchReward, result.score);
         assertEquals(s.indexOf(expectedOverlap), result.originateI);
         assertEquals(0, result.originateJ);
    }

    @Test
    public void overlapAlignTestWithErrors() {
         // In overlap alignment, the end of s overlaps with the beginning of t, so we start at the beginning of t, but start anywhere in s. Errors (mismatch and indels) are still allowed.
         String s = "Chapter 5 of Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach 3rded Active Learning Publishers (2018)";
         String t = "a active learning aproach, 3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/";
         String expectedOverlap1 = "an active learning approach- 3rd-ed- Active Learning Publishers- (2018)";
         String expectedOverlap2 = "a- active learning ap-roach, 3rd ed. Active Learning Publishers, -2018-";
         
         String overlapExpectedSpotInS = "an active learning approach";
         int matchReward = 1;
         int mismatchPenalty = 5; // adding a lot of errors here but don't want it to give up, so lower this to 5
         int indelPenalty = 5; // adding a lot of errors here but don't want it to give up, so lower this to 5
         AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, s, t);
         assertEquals((expectedOverlap1.length() - 8)*matchReward - 8*5, result.score);
         assertEquals(s.indexOf(overlapExpectedSpotInS), result.originateI);
         assertEquals(0, result.originateJ);
    }

    // some edge cases not used in the current work
    @Test
    public void overlapAlignThrowsIllegalArgumentExceptionForEmptyOrNullStrings() {
        int matchReward = 1;
        int mismatchPenalty = 5;
        int indelPenalty = 5;
        IllegalArgumentException exceptionSNull = assertThrows(IllegalArgumentException.class, () -> {
                AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, null, "abcd");
        });
        assertEquals("Aligning empty or null strings s or t is not supported", exceptionSNull.getMessage());

        IllegalArgumentException exceptionSEmpty = assertThrows(IllegalArgumentException.class, () -> {
                AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "", "abcd");
        });
        assertEquals("Aligning empty or null strings s or t is not supported", exceptionSEmpty.getMessage());

        IllegalArgumentException exceptionTNull = assertThrows(IllegalArgumentException.class, () -> {
                AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "aa", null);
        });
        assertEquals("Aligning empty or null strings s or t is not supported", exceptionTNull.getMessage());

        IllegalArgumentException exceptionTEmpty = assertThrows(IllegalArgumentException.class, () -> {
                AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "aa", "");
        });
        assertEquals("Aligning empty or null strings s or t is not supported", exceptionTEmpty.getMessage());
    }

    @Test
    public void overlapAlignTestTwoCharacterS() {
        int matchReward = 1;
        int mismatchPenalty = 5;
        int indelPenalty = 5;
        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "aa", "abcd");
        assertEquals(1, result.originateI);
        assertEquals(0, result.originateJ);
    }

    // we are not aligning single character strings,
    // which in overlap alignment for the source string, is only allowed one answer, (0,0), 
    // but want to support it
    @Test
    public void overlapAlignTestOneCharacterS() {
        int matchReward = 1;
        int mismatchPenalty = 5;
        int indelPenalty = 5;
        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "a", "abcd");
        // should be (0,0) as only possible answer for all possible s,t if s is length 1
        assertEquals(0, result.originateI);
        assertEquals(0, result.originateJ);
        assertEquals(matchReward, result.score);
    }

    // originally I had this case suggesting to start aligning t after the end of s, which conveyed the right information but it is clearer to return -1 to indicate no overlap
    @Test
    public void overlapAlignTestOneCharacterSMismatch() {
        int matchReward = 1;
        int mismatchPenalty = 5;
        int indelPenalty = 5;
        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "f", "abcd");
        assertEquals(-1, result.originateI);
        assertEquals(-1, result.originateJ);
    }

    @Test
    public void overlapAlignTestTwoCharacterSMismatch() {
        int matchReward = 1;
        int mismatchPenalty = 5;
        int indelPenalty = 5;
        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "ff", "abcd");
        assertEquals(-1, result.originateI);
        assertEquals(-1, result.originateJ);
    }

    @Test
    public void overlapAlignTestOneCharacterTMismatch() {
        int matchReward = 1;
        int mismatchPenalty = 5;
        int indelPenalty = 5;
        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "abcd", "f");
        assertEquals(-1, result.originateI);
        assertEquals(-1, result.originateJ);
    }

    @Test
    public void overlapAlignTestTwoCharacterTMismatch() {
        int matchReward = 1;
        int mismatchPenalty = 5;
        int indelPenalty = 5;
        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "abcd", "ff");
        assertEquals(-1, result.originateI);
        assertEquals(-1, result.originateJ);
    }

    @Test
    public void overlapAlignStringWithItself() {
        int matchReward = 1;
        int mismatchPenalty = 5;
        int indelPenalty = 5;
        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "abcd", "abcd");
        assertEquals(0, result.originateI);
        assertEquals(0, result.originateJ);
        assertEquals("abcd".length()*matchReward, result.score);
    }

    @Test
    public void overlapAlignTestTThatNestsWithinS() {
        int matchReward = 10;
        int mismatchPenalty = 1;
        int indelPenalty = 1;

        AlignmentResult result1 = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "abcd", "a");
        assertEquals(0, result1.originateI);
        assertEquals(0, result1.originateJ);
        assertEquals(1*matchReward - 3*indelPenalty, result1.score);

        AlignmentResult result2 = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "abcd", "b");
        assertEquals(1, result2.originateI);
        assertEquals(0, result2.originateJ);
        // skipping first character is free (originate at position 1)
        assertEquals(1*matchReward - 2*indelPenalty, result2.score);

        AlignmentResult result3 = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "abcd", "c");
        assertEquals(2, result3.originateI);
        assertEquals(0, result3.originateJ);
        assertEquals(1*matchReward - 1*indelPenalty, result3.score);

        AlignmentResult result4 = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "abcd", "d");
        assertEquals(3, result4.originateI);
        assertEquals(0, result4.originateJ);
        assertEquals(1*matchReward - 0*indelPenalty, result4.score);
    }

    // Not a case we use this class for,
    // overlap alignment is a suffix of string s with a prefix of string t so alignment should NOT begin before the start of s (originateI cannot be negative except -1 to indicate no alignment and we use originateJ as zero for overlap alignment).
    // I could imagine redefining the problem to support this but it is a different kind of alignment (see definitions at top of file), comparing overlap alignment vs local alignment, global alignment, fitting alignment, etc.
    @Test
    public void overlapAlignTestSThatNestsWithinT() {
        int matchReward = 10;
        int mismatchPenalty = 1;
        int indelPenalty = 1;

        AlignmentResult result1 = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "a", "abcd");
        assertEquals(0, result1.originateI);
        assertEquals(0, result1.originateJ);
        assertEquals(1*matchReward, result1.score);

        AlignmentResult result2 = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "b", "abcd");
        assertEquals(-1, result2.originateI); // would be 0 if allowed
        assertEquals(-1, result2.originateJ);
        // skipping first character is NOT free for overlap alignment moving s string (scoring starts from beginning of t and ends at the end of s as aligning suffix of s with prefix of t), if we did allow it would have penalty below
        assertEquals(1*matchReward - 1*indelPenalty, result2.score);

        AlignmentResult result3 = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "c", "abcd");
        assertEquals(-1, result3.originateI); // would be 0 if allowed
        assertEquals(-1, result3.originateJ);
        // skipping first characters are NOT free for overlap alignment moving s string (scoring starts from beginning of t and ends at the end of s as aligning suffix of s with prefix of t), if we did allow it would have penalty below
        assertEquals(1*matchReward - 2*indelPenalty, result3.score);

        AlignmentResult result4 = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "d", "abcd");
        assertEquals(-1, result4.originateI); // would be 0 if allowed
        assertEquals(-1, result4.originateJ);
        // skipping first characters are NOT free for overlap alignment moving s string (scoring starts from beginning of t and ends at the end of s as aligning suffix of s with prefix of t), if we did allow it would have penalty below
        assertEquals(1*matchReward - 3*indelPenalty, result4.score);
    }

    // trivial (not relevant to this project) edge case of overlapping two single character strings
    // only possible overlap is (0,0)
    @Test
    public void overlapAlignTestSandTBoth1CharacterMatch() {
        int matchReward = 5;
        int mismatchPenalty = 1;
        int indelPenalty = 1;

        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "a", "a");
        assertEquals(0, result.originateI);
        assertEquals(0, result.originateJ);
        assertEquals(matchReward, result.score);
    }

    @Test
    public void overlapAlignTestSandTBoth1CharacterMismatch() {
        int matchReward = 5;
        int mismatchPenalty = 1;
        int indelPenalty = 1;

        AlignmentResult result = OverlapAlignment.align(matchReward, mismatchPenalty, indelPenalty, "a", "b");
        assertEquals(-1, result.originateI);
        assertEquals(-1, result.originateJ);
        assertEquals(-mismatchPenalty, result.score);
    }
}
