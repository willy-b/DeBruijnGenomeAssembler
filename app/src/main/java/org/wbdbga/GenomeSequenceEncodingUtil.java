/*
GenomeSequenceEncodingUtil for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
Series of utilities to transform genome sequence into a compressed representation closer to the ideal 2 bits per base pair, and to translate kmers into 4 byte pointers (3 bytes specifying a read where the kmer sequence occurred and 1 byte specifying the offset within that read) and back again.

The author was originally developing this on a budget laptop with a couple of gigabytes of ram and thus found it necessary to approach the optimal 2 bits per basepair for representing the genetic sequence, and also avoid explicitly representing sequences, e.g. for kmers, more than once by working with bytes that contain the read index and offset of the first place we saw that kmer. These utilities were developed to achieve these goals.
*/
public class GenomeSequenceEncodingUtil {
    static char[] ALPHABET = new char[]{'A', 'C', 'G', 'T'};

    /**
       Accepts a base-pair byte array produced by getBPBytesFromString and reads out the length in base pairs declared in the first byte.
       Should match lineBytes.length * 4 - 4 .
       @see org.wbdbga.GenomeSequenceEncodingUtil#getBPBytesFromString

       @param lineBytes byte[] representation of a read (or other line of genetic sequence less than or equal to 255 bases long) encoded by GenomeSequenceEncodingUtil#getBPBytesFromString
       @return declared length in bases of the encoded sequence
     */
    public static int getDeclaredBPLength(byte[] lineBytes) {
        // note: forces unsigned
        return (int)(lineBytes[0] & 0b1111111); // should match lineBytes.length * 4 - 4
    }

    /**
       Returns the base pair sequence string encoded in the provided byte array (byte[]) which should have been encoded by getBPBytesFromString .
       @see org.wbdbga.GenomeSequenceEncodingUtil#getBPBytesFromString

       @param lineBytes byte[] representation of a read (or other line of genetic sequence less than or equal to 255 bases long) encoded by GenomeSequenceEncodingUtil#getBPBytesFromString
       @return the original read or line in ATCG string format that had been encoded into the byte[] passed in lineBytes
     */
    public static String getStringFromBPBytes(byte[] lineBytes) {
        int numCharsInLine = getDeclaredBPLength(lineBytes); // should match lineBytes.length * 4 - 4
        if (numCharsInLine + 4 > 4 * lineBytes.length) {
            throw new IllegalArgumentException("Malformed BP bytes array: first byte indicates a number of base pairs that could not fit in the array provided");
        }

        // currently unreachable. but just in case the previous lines change.
        if (numCharsInLine <= 0) {
            throw new IllegalArgumentException("Invalid BP bytes array, first byte indicates a negative number of base pairs in the string");
        }

        char[] chars = new char[numCharsInLine];

        for (int i = 0; i < numCharsInLine; i++) {
            byte b = lineBytes[i / 4 + 1]; // add one to skip lineBytes[0] which resolves ambiguity in last byte of whether it has 1,2,3 or 4 characters in case of L % 4 != 0
            b = (byte)((b >> ((i % 4) * 2)) & 0b11);
            chars[i] = (ALPHABET[(int)b]);
        }

        return new String(chars);
    }

    /**
    Encodes the provided base pair sequence string (must contain only 'A','C','T','G' characters or InvalidCharacterInBPStringException will be thrown; must be less than or equal to 255 base pairs long) as a byte array (byte[]) with the first byte encoding the number of base pairs and subsequent bytes encoding up to 4 bases each.
    The resulting byte[] can be deooded using getStringFromBPBytes .
    @see org.wbdbga.GenomeSequenceEncodingUtil#getStringFromBPBytes

    @param line String (ATCG characters only) representation of a read (or other line of genetic sequence less than or equal to 255 bases long)
    @return byte[] encoded version of the read/line using 1 byte to encode the length, and then 2 bits per base after that (up to 4 bases per byte)
    */
    public static byte[] getBPBytesFromString(String line) {
        //if (line.length() % 4 != 0) { throw new IllegalArgumentException("while we can support read lengths not divisible by 4, it requires using extra bits to track how many characters we used, so we are not going to do that. typical read lengths are divisible by 4. we can truncate when this is not the case."); }

        if (line.length() > 255) {
            throw new TooLongBPReadStringToEncodeUsingBPBytesFromStringException("getBPBytesFromString is intended to encode reads (usually 100bp or so) and kmers (shorter than reads) and so currently only supports individual strings up to 255 bases long.");
        }
        int numBytes = ((line.length() * 2) / 8) + ((line.length() * 2) % 8 != 0 ? 1 : 0);
        numBytes++; // add the byte to track the number of characters to eliminate the ambiguity in the (L-1)th byte (does it have 1,2,3 or 4 characters if L % 4 != 0?)
        byte[] lineBytes = new byte[numBytes];
        Arrays.fill(lineBytes, (byte)0);
        lineBytes[0] = (byte)line.length(); 
        // [A,C,G,T] == 00->11
        for (int i = 0; i < line.length(); i++) {
            byte b;
            switch (line.charAt(i)) {
            case 'A':
                b = 0;
                break;
            case 'C':
                b = 1;
                break;
            case 'G':
                b = 2;
                break;
            case 'T':
                b = 3;
                break;
            default:
                throw new InvalidCharacterInBPStringException("Unexpected character '" + line.charAt(i) + "', must be one of 'A','C','G','T'");
            }

            lineBytes[i/4 + 1] |= (b << ((i % 4)*2)); // add one since lineBytes[0] is used to capture number of characters for decoding side.
        }

        return lineBytes;
    }

    /**
    Refer to genome sequence strings using 4 bytes: assuming we don't support more than 2^24 (3 bytes) reads and length 255 reads (1 byte), we can save a factor of 2 on integer pair representation, which is already a large savings vs keeping copies of the actual kmers or whatever genome sequence strings we are referring to. Decode from returned byte array using getReadIndexAndWithinReadIndexFrom4ByteRepresentation . Note that there is no length stored, it is implied as this method is used to get pointers for kmers and (k-1)mers where k is known and does not need to be stored repeatedly for each encoded k/(k-1)mer.
    @see org.wbdbga.GenomeSequenceEncodingUtil#getReadIndexAndWithinReadIndexFrom4ByteRepresentation

    @param readIndex the index of the read which provides an example of the sequence being represented that we will reference
    @param withinReadIndex the index/offset within the read for the beginning of the fixed length sequence being referenced
    @return a byte[] encoding a reference to a sequence in a particular read at a particular offset (3 bytes for read index, 1 byte for within-read offset)
    */
    public static byte[] get4ByteRepresentationForReadIndexAndWithinReadIndex(int readIndex, int withinReadIndex) {
        if (readIndex < 0 || readIndex >= (1 << 24)) {
            throw new ReadIndexOutOfSupportedRangeFor4ByteSequenceInReadReferenceEncoding("get4ByteRepresentationForReadIndexAndWithinReadIndex only supports non-negative read indices, and since 3 byte encoding for read index is used, 2^24 reads is the maximum");
        }
        if (withinReadIndex < 0 || withinReadIndex > 255) {
            throw new WithinReadIndexOutOfSupportedRangeFor4ByteSequenceInReadReferenceEncoding("get4ByteRepresentationForReadIndexAndWithinReadIndex only supports non-negative within read indices, and since 1 byte encoding for the within read index is used, only 0-255 is supported for within read (consistent with getBPBytesFromString which only encodes up to 255 base long reads)");
        }
        byte[] result = new byte[4];

        // store the read index
        result[0] = (byte)(readIndex & 0xff);
        result[1] = (byte)((readIndex >> 8) & 0xff);
        result[2] = (byte)((readIndex >> 16) & 0xff);

        // store the within read index
        result[3] = (byte)(withinReadIndex & 0xff);

        return result;
    }

    /**
    We refer to genome sequence strings using 4 bytes which is already a large savings vs keeping copies of the actual kmers or whatever genome sequence strings we are referring to, this method decodes the byte array returned from get4ByteRepresentationForReadIndexAndWithinReadIndex into a read index and within read index for looking up some fixed length sequence (say a kmer of length k or a de Bruijn node (k-1)mer.

    @see org.wbdbga.GenomeSequenceEncodingUtil#get4ByteRepresentationForReadIndexAndWithinReadIndex

    @param bytes a reference to a sequence in a particular read at a particular offset encoded as byte[] by get4ByteRepresentationForReadIndexAndWithinReadIndex
    @return int[] containing a pair of ints, the read index and the within read index for an instance of some sequence of interest, like a kmer or (k-1)mer.
    */
    public static int[] getReadIndexAndWithinReadIndexFrom4ByteRepresentation(byte[] bytes) {
        int readIndex = 0;
        int withinReadIndex = 0;

        readIndex |= ((int)bytes[0]) & 0xff;
        readIndex |= (((int)bytes[1]) & 0xff) << 8;
        readIndex |= (((int)bytes[2]) & 0xff) << 16;

        withinReadIndex |= ((int)bytes[3]) & 0xff;

        return new int[]{readIndex, withinReadIndex};
    }

    /**
       Forward compatibility wrapper for upcoming 5 byte paired read support to abstract the representation of read and within read index. Paired reads are not in this version so this just forwards to getReadIndexAndWithinReadIndexFrom4ByteRepresentation .
       @see org.wbdbga.GenomeSequenceEncodingUtil#getReadIndexAndWithinReadIndexFrom4ByteRepresentation

       @param kmerBytes a reference to a sequence in a particular read at a particular offset encoded as byte[] by get4ByteRepresentationForReadIndexAndWithinReadIndex or another encoding method using a different number of bytes like a forthcoming one for paired read representation
       @return int[] containing a pair of ints, the read index and the within read index for an instance of some sequence of interest, like a kmer or (k-1)mer.
    */
    public static int[] decodeKmerBytesToIntArray(byte[] kmerBytes) {
        int[] result;
        if (kmerBytes.length == 4) {
            result = GenomeSequenceEncodingUtil.getReadIndexAndWithinReadIndexFrom4ByteRepresentation(kmerBytes);
        } else {
            throw new IllegalArgumentException("Paired reads or other usecases for kmer bytes of length 5 or otherwise not length 4 are not currently supported!");
        }
        return result;
    }

    static class TooLongBPReadStringToEncodeUsingBPBytesFromStringException extends RuntimeException {
        public TooLongBPReadStringToEncodeUsingBPBytesFromStringException(String message) {
            super(message);
        }
    }

    static class InvalidCharacterInBPStringException extends RuntimeException {
        public InvalidCharacterInBPStringException(String message) {
            super(message);
        }
    }

    static class ReadIndexOutOfSupportedRangeFor4ByteSequenceInReadReferenceEncoding extends RuntimeException {
        public ReadIndexOutOfSupportedRangeFor4ByteSequenceInReadReferenceEncoding(String message) {
            super(message);
        }
    }

    static class WithinReadIndexOutOfSupportedRangeFor4ByteSequenceInReadReferenceEncoding extends RuntimeException {
        public WithinReadIndexOutOfSupportedRangeFor4ByteSequenceInReadReferenceEncoding(String message) {
            super(message);
        }
    }
}
