/*
ValueComparedByteArrayTest for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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

import java.util.Arrays;

class ValueComparedByteArrayTest {
    @Test
    public void valueComparedByteArrayConstructor() {
        byte[] data = new byte[]{0,0,0,1};
        ValueComparedByteArray vcba = new ValueComparedByteArray(data);
        assertNotNull(vcba);
        assertEquals(data, vcba.value);
    }

    @Test
    public void valueComparedByteArrayConstructorThrowsIllegalArgumentExceptionIfAttemptingToWrapANullArray() {
        byte[] data = null;
        IllegalArgumentException exception = assertThrows(IllegalArgumentException.class, () -> {
                ValueComparedByteArray vcba = new ValueComparedByteArray(data);
            });
        assertEquals("ValueComparedByteArray should not be used to wrap null byte arrays", exception.getMessage());
    }

    @Test
    public void equalsWithValueComparedByteArray() {
        byte[] data = new byte[]{0,1,0,1};
        byte[] separateDataWithSameValue = new byte[]{0,1,0,1};
        byte[] dataWithDifferentValue = new byte[]{0,1,0,0};
        assertTrue((new ValueComparedByteArray(data)).equals(new ValueComparedByteArray(data)));
        assertEquals(new ValueComparedByteArray(data), new ValueComparedByteArray(data));
        assertTrue((new ValueComparedByteArray(data)).equals(new ValueComparedByteArray(separateDataWithSameValue)));
        assertEquals(new ValueComparedByteArray(data), new ValueComparedByteArray(separateDataWithSameValue));
        assertFalse((new ValueComparedByteArray(data)).equals(new ValueComparedByteArray(dataWithDifferentValue)));
        assertNotEquals(new ValueComparedByteArray(data), new ValueComparedByteArray(dataWithDifferentValue));
    }

    @Test
    public void equalsWithNonValueComparedByteArrayReturnsFalse() {
        byte[] data = new byte[]{0,1,0,1};
        byte[] separateDataWithSameValue = new byte[]{0,1,0,1};
        byte[] dataWithDifferentValue = new byte[]{0,1,0,0};
        assertFalse((new ValueComparedByteArray(data)).equals(data));
        assertNotEquals(new ValueComparedByteArray(data), data);
        assertFalse((new ValueComparedByteArray(data)).equals(separateDataWithSameValue));
        assertNotEquals(new ValueComparedByteArray(data), separateDataWithSameValue);
        assertFalse((new ValueComparedByteArray(data)).equals(dataWithDifferentValue));
        assertNotEquals(new ValueComparedByteArray(data), dataWithDifferentValue);
    }

    @Test
    public void valueComparedByteArrayHashCode() {
        byte[] data = new byte[]{0,1,0,1};
        byte[] separateDataWithSameValue = new byte[]{0,1,0,1};
        byte[] dataWithDifferentValue = new byte[]{0,1,0,0};
        // hashCodes for same values should be same
        assertEquals((new ValueComparedByteArray(data)).hashCode(),(new ValueComparedByteArray(data)).hashCode());
        assertEquals((new ValueComparedByteArray(data)).hashCode(),(new ValueComparedByteArray(separateDataWithSameValue)).hashCode());
        assertNotEquals((new ValueComparedByteArray(data)).hashCode(),(new ValueComparedByteArray(dataWithDifferentValue)).hashCode());

        // override of hashCode should make it equal to the Arrays.hashCode of underlying data (not related to a hash of the object reference)
        assertEquals(Arrays.hashCode(data), (new ValueComparedByteArray(data)).hashCode());
        assertEquals(Arrays.hashCode(dataWithDifferentValue), (new ValueComparedByteArray(dataWithDifferentValue)).hashCode());
    }

    @Test
    public void valueComparedByteArrayToString() {
        byte[] data = new byte[]{0,1,0,1};
        assertEquals(Arrays.toString(data), (new ValueComparedByteArray(data)).toString());
    }
}
