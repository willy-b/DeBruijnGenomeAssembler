/*
ValueComparedByteArray for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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

import java.util.Arrays;
import java.io.Serializable;

/**
Safely use byte arrays in sets and as keys in maps with .equals comparing values instead of references (when necessary).
<p>
To avoid wasting precious memory on low-powered laptops and approach 2-bits per base pair for the first copy of the sequence and efficiently point into that single copy via 4 byte pointers (3 bytes for read index, 1 byte for offset within the read), we use byte arrays (see GenomeSequenceEncodingUtil) to represent kmer pointers, reads, etc.
<p>
But so we can safely use byte arrays in sets and as keys in maps with .equals comparing values instead of references (when necessary).
<p>
In an early version of this I just assidiuously avoided creating multiple objects for the same byte array so that the reference comparison was safe, but this is very fragile, easy to break, and hard to write unit tests for, so we wrap byte[] in ValueComparedByteArray below to make it simple and safe.
*/
public final class ValueComparedByteArray implements Serializable {
    byte[] value;

    /**
       Construct a ValueComparedByteArray from a byte[] , wrapping it so it can safely be used as a Map key or Set entry without being forced to use only a single instance (by default the .equals would compare only references), making it safe for multiple copies at the cost of extra overhead.
       @param value byte[] to wrap in a ValueComparedByteArray, subsequent .equals will return true for and .hashCode results will be equal for all other instances constructed from any byte[] value equivalent to the value argument
     */
    public ValueComparedByteArray(byte[] value) {
        if (value == null) {
            throw new IllegalArgumentException("ValueComparedByteArray should not be used to wrap null byte arrays");
        }
        this.value = value;
    }

    public boolean equals(Object o) {
        if (o instanceof ValueComparedByteArray) {
            return Arrays.equals(this.value, ((ValueComparedByteArray)o).value);
        }
        return false;
    }

    /**
       @see java.util.Arrays#hashCode
       @return int hashCode computed by Arrays.hashCode from the byte[] this ValueWrappedByteArray was constructed from
     */
    public int hashCode() {
        return Arrays.hashCode(this.value);
    }

    /**
       Static method alternative to using the constructor for getting a ValueComparedByteArray from a byte[] , wrapping it so it can safely be used as a Map key or Set entry without being forced to use only a single instance (by default the .equals would compare only references), making it safe for multiple copies at the cost of extra overhead.
       @param b byte[] to wrap in a ValueComparedByteArray, subsequent .equals calls will return true for and .hashCode will be equal for all other instances constructed from any byte[] with the same values as b
       @return ValueComparedByteArray wrapped version of the byte[] passed in (.value field set to b argument)
     */
    public static ValueComparedByteArray from(byte[] b) {
        return new ValueComparedByteArray(b);
    }

    /**
       returns a string showing the value of the byte[] instead of the reference (uses Arrays.toString)
       @see java.util.Arrays#toString
     */
    public String toString() {
        return Arrays.toString(this.value);
    }
}
