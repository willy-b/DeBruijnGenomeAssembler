/*
OverlapAlignment for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
In overlap alignment, the end of s overlaps with the beginning of t, so we start at the beginning of t, but start anywhere in s. Errors (mismatch and indels) are still allowed.
Overlap alignment is NOT the same as edit distance, fitting alignment, local alignment, or global alignment.
<p>
In this project, we do NOT use overlap alignment directly for assembly, but use the more efficient de Bruijn graph based assembly, so why is this dynamic programming based aligner here?
This overlap aligner is used by this assembler in a final cleanup/polishing step to find and remove circular overlap (less than or equal to k in length) at the end of contigs for sequences from circular genomes (it is not expensive to perform overlap alignment just on the up to k bases at the ends of each contig to see if they circularly align with themselves to trim off any such circular content.)
<p>
See Chapter 5 of Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach, 3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/ , for a definition of overlap alignment and comparison with edit distance, fitting alignment, local alignment, and global alignment.
*/

// Credit to UCSD ALGS205x for teaching me about these alignment techniques.
// The textbook to that course contains the same theory, see Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach, 3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/.

public class OverlapAlignment {
    /**
Perform overlap alignment of two strings, s and t, and return the score, the end of s overlaps with the beginning of t, so we start at the beginning of t, but start anywhere in s.

@param m match score
@param mu mismatch score (negative sign implied)
@param sigma - indel score (negative sign implied)
@param s string that begins the aligned string, a suffix of this is matched to a prefix of string t
@param t string that ends the aligned string, a prefix of this is matched to a suffix of string s

@return AlignmentResult containing the score and origination index of the overlap alignment; if no overlap at all is better than any possible overlap then -1 will be returned for originateI and originateJ.
     */
    public static AlignmentResult align(int m, int mu, int sigma, String s, String t) {
        // edge cases we can reject or answer without processing
        if (s == null || t == null || s.length() == 0 || t.length() == 0) {
            throw new IllegalArgumentException("Aligning empty or null strings s or t is not supported");
        }
        
        s = " " + s;
        t = " " + t;

        int len1 = s.length();
        int len2 = t.length();

        int[][] mat = new int[len1][len2];
        char[][] choiceChars = new char[len1][len2];
        int terminateI = -1;
        int terminateJ = -1;

        final char INSERTION = 0;
        final char DELETION = 1;
        final char MATCH = 2;
        final char MISMATCH = 3;
        final char ORIGINATE = 4; // unlike local alignment we can only originate anywhere j == 0 (we can start anywhere in s but must start at beginning of t)
        final char TERMINATE = 5; // unlike local alignment we can only terminate anywhere i = len1 - 1  (we can terminate anywhere in t but must terminate at the end of s)
        final char UNSET_INVALID = 6;

        // m - match score
        // mu - mismatch score (negative sign implied)
        // sigma - indel score (negative sign implied)
        // originate - free
        // terminate - free

        mat[0][0] = 0;
        // go column by column
        for (int j = 0; j < len2; j++) {
            for (int i = 0; i < len1; i++) {
                if (i == 0 && j == 0) {
                    mat[i][j] = 0;
                    continue;
                }

                char choiceChar = UNSET_INVALID;
                int highestScore = Integer.MIN_VALUE;
                if (j > 0) {
                    highestScore = mat[i][j - 1] - sigma;
                    choiceChar = INSERTION;
                }

                if (i > 0) {
                    if (mat[i - 1][j] - sigma > highestScore) {
                        highestScore = mat[i - 1][j] - sigma;
                        choiceChar = DELETION;
                    }
                }

                if (i > 0 && j > 0) {
                    // consider match
                    if (s.charAt(i) == t.charAt(j) && mat[i - 1][j - 1] + m > highestScore) {
                        highestScore = mat[i - 1][j - 1] + m;
                        choiceChar = MATCH;
                    }
                    // consider mismatch
                    if (s.charAt(i) != t.charAt(j) && mat[i - 1][j - 1] - mu > highestScore) {
                        highestScore = mat[i - 1][j - 1] - mu;
                        choiceChar = MISMATCH;
                    }
                }

                // consider originate
                // unlike local alignment we can only originate anywhere j == 0 (we can start anywhere in s but must start at beginning of t)
                if (j == 0 && mat[0][0] > highestScore) {
                    highestScore = mat[0][0];
                    choiceChar = ORIGINATE;
                }

                if (i == len1 - 1 && j == len2 - 1) {
                    // check terminate inputs
                    // unlike local alignment we can only terminate anywhere i = len1 - 1  (we can terminate anywhere in t but must terminate at the end of s)
                    int ii = len1 - 1;
                    for (int jj = 0; jj < len2; jj++) {
                        /*
                        // `ii = len1 - 1` cannot be 0 because len1 is (" " + s).length() and s.length() > 0 , so len1 >= 2.
                        // leaving this comment here so someone does not edit the above lines inadvertently to make this untrue,
                        // and because for other types of alignments where ii was not fixed at len1 - 1, we would need this.
                        if (ii == 0 && jj == 0) {
                            continue;
                        }
                        */

                        if (ii == i && jj == j) {
                            continue;
                        }

                        if (mat[ii][jj] > highestScore && choiceChars[ii][jj] != ORIGINATE) {
                            highestScore = mat[ii][jj];
                            terminateI = ii;
                            terminateJ = jj;
                            choiceChar = TERMINATE;
                        }
                    }
                }

                mat[i][j] = highestScore;
                choiceChars[i][j] = choiceChar;
            }
        }

        int i = len1 - 1;
        int j = len2 - 1;

        // sb1, sb2 are left in here in case someone wants to uncomment
        // the lines later down to print out the alignment details
        StringBuilder sb1 = new StringBuilder("");
        StringBuilder sb2 = new StringBuilder("");

        int originateI = -1;
        int originateJ = -1;
        while (j > 0 || i > 0) {
            char choiceChar = choiceChars[i][j];
            if (choiceChar == INSERTION) {
                sb1.insert(0, "-");
                sb2.insert(0, t.charAt(j));
                j--;
            } else if (choiceChar == DELETION) {
                sb1.insert(0, s.charAt(i));
                sb2.insert(0, "-");
                i--;
            } else if (choiceChar == MATCH || choiceChar == MISMATCH) {
                sb1.insert(0, s.charAt(i));
                sb2.insert(0, t.charAt(j));
                i--; j--;
                if (i == 0 && j == 0) {
                    originateI = i;
                    originateJ = j;
                }
            } else if (choiceChar == ORIGINATE) {
                //System.out.println("Originate: (" + i + ", " + j + ")");
                originateI = i;
                originateJ = j;
                i = 0;
                j = 0;
                continue;
            } else if (choiceChar == TERMINATE) {
                //System.out.println("Terminate: (" + terminateI + ", " + terminateJ + ")");
                i = terminateI;
                j = terminateJ;
                continue;
            } else {
                throw new InvalidAlignmentStepChoice("Error unexpected choiceChar " + choiceChar + " at i==" + i + ", j==" + j +".");
            }
        }

        // If debugging or curious about a particular alignment, it may be useful to uncomment the following two lines:
        // ```
        // System.out.println("aligned s: " + sb1);
        // System.out.println("aligned t: " + sb2);
        // ```

        // For e.g. s shorter than t, alignment beginning after the end of s means it did not align
        // subtract extra 1 for " " appended in beginning
        if (originateI > s.length() - 1 - 1) {
            // originally I had this case suggesting to start aligning t after the end of s to indicate no overlap
            // which conveyed the right information, but it is clearer to return -1 (an impossible value) to indicate no overlap
            originateI = -1;
            originateJ = -1; // instead of 0
        }
        return new AlignmentResult(mat[len1 - 1][len2 - 1], originateI, originateJ);
    }

    /**
       Exception class that will only be invoked due to coding errors made in updates to this class that introduce unsupported alignment step choices.
     */
    static class InvalidAlignmentStepChoice extends RuntimeException {
        public InvalidAlignmentStepChoice(String message) {
            super(message);
        }
    }
}
