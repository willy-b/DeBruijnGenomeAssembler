/*
KmerTrieFilterReads for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
Iterative repair of reads using the kmer spectrum before building the de Bruijn graph.
<p>
As part of new iterative read repair process, we build a trie from the kmers we read out after dropping low coverage kmers. We use this class to process each read and find if there is a single position we can update to get a large increase in spectral alignment against the filtered kmer spectrum.
<p>
If there is, this class will recommend that update, then, outside this class, the read can be updated, and then the kmers will be recomputed from the filtered reads. That way low coverage kmers which were dropped which are not easily resolved by a single point error correction will be restored but all low coverage kmers which can be avoided by single point error corrections will be dropped.
<p>
This approach was inspired by the discussion of pre-assembly read error correction mentioned as used by the EULER assembler in P.A. Pevzner, H. Tang, M.S. Waterman, An Eulerian path approach to DNA fragment assembly, Proc. Natl. Acad. Sci. U.S.A. 98 (17) 9748-9753, https://doi.org/10.1073/pnas.171285098 (2001), under the heading Error Correction and Data Corruption .
*/
public final class KmerTrieFilterReads {
    List<byte[]> reads; // modifies the reads

    // note that no .contains or operations relying on .equals of byte[] is used with this set (Arrays.equals must be used to compare byte[] , but I wrap it in ValueComparedByteArray anyway to avoid concern and in case we want to do that in the future)
    // it is just a set to iterate over
    Set<ValueComparedByteArray> kmers; // does NOT modify this set

    int k;
    int MISALIGNED_KMER_THRESHOLD;
    int MIN_KMER_BREAKS_REQUIRED_CORRECTED_PER_UPDATE; // k is a safe value
    int numFakeReadsAddedForOrphanedKmers;

    List<Map<Character, Integer>> trie;

    Set<Integer> spectrallyMisalignedReadsToIgnore = new HashSet<Integer>();
    Set<ReadCorrection> readCorrectionsToMake = new HashSet<ReadCorrection>();

    static final class ReadCorrection {
        int readIdx;
        int withinReadIdx;
        char newCharacter;

        public ReadCorrection(int readIdx, int withinReadIdx, char newCharacter) {
            this.readIdx = readIdx;
            this.withinReadIdx = withinReadIdx;
            // could also be an enum or some other way to restrict the alphabet
            this.newCharacter = newCharacter;
            if (newCharacter != 'A' && newCharacter != 'T' && newCharacter != 'C' && newCharacter != 'G') {
                throw new InvalidCharacterException("ReadCorrections must specify a new character of 'A','T','C', or 'G', found '" + newCharacter + "'");
            }
        }
        static class InvalidCharacterException extends RuntimeException {
            public InvalidCharacterException(String message) {
                super(message);
            }
        }
    }
    
    /**
      Construct KmerTrieFilterReads from a list of reads, a set of kmers, the k value, and the number of fake reads.
      It can then be used to recommend ReadCorrections for the iterative repair of reads using the kmer spectrum before building the de Bruijn graph.

      @param reads List of GenomeSequenceEncodingUtil.getBPBytesFromString encoded reads 
      @param kmers Set of ValueComparedByteArray-wrapped 4 byte kmer references as returned by GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex
      @param k Length k of the kmers (a kmer is a substring of length k shorter than the read length for which we expect to have coverage in our samples of all unique values present in the genome being sequenced/assembled). In the de Bruijn graph the edges will represent the kmers of the reads and the nodes represent (k-1)mers.
      @param numFakeReadsAddedForOrphanedKmers number of fake reads to preserve GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex references in the case where a ReadCorrection removes the kmer from the originally referenced read but we know it is still present elsewhere (but don't want to slow down to find another place in the reads it occurs), see DeBruijnKmers.unloadRead
      @see org.wbdbga.GenomeSequenceEncodingUtil#getBPBytesFromString
      @see org.wbdbga.GenomeSequenceEncodingUtil#get4ByteRepresentationForReadIndexAndWithinReadIndex
      @see org.wbdbga.DeBruijnKmers#unloadRead
     */
    public KmerTrieFilterReads(List<byte[]> reads, Set<ValueComparedByteArray> kmers, int k, int numFakeReadsAddedForOrphanedKmers) {
        if (reads == null) {
            throw new IllegalArgumentException("reads argument must not be null");
        }
        if (k <= 0) {
            throw new IllegalArgumentException("Invalid value of k");
        }
        if (kmers == null) {
            throw new IllegalArgumentException("kmers argument must not be null");
        }

        this.reads = reads;
        // note that no .contains or operations relying on .equals of byte[] is used with this set (Arrays.equals must be used to compare byte[] , or they should be wrapped in ValueComparedByteArray)
        // it is just a set to iterate over
        this.kmers = kmers;
        this.k = k;

        MISALIGNED_KMER_THRESHOLD = 4 * k;
        MIN_KMER_BREAKS_REQUIRED_CORRECTED_PER_UPDATE = k / 4;

        trie = new ArrayList<Map<Character, Integer>>();
        this.numFakeReadsAddedForOrphanedKmers = numFakeReadsAddedForOrphanedKmers;
    }

    /**
       Build a trie of depth k from all solid (well covered) kmers for recommending ReadCorrections in recommendReadCorrections()/checkRead(). This will be used to determine at which positions in reads that are not well aligned to the kmer spectrum could have minimal substitutions to align them with the kmer spectrum, by walking the trie from each offset and recording the positions at which the trie traversal fails and the alternate paths that were available in that trie node for that positions. To understand the structure of the trie, see KmerTrieFilterReadsTest#buildTrie .

       @see org.wbdbga.KmerTrieFilterReads#recommendReadCorrections
       @see org.wbdbga.KmerTrieFilterReads#checkRead
     */
    public void buildTrie() {
        // root node
        trie.add(new HashMap<Character, Integer>());

        for (ValueComparedByteArray kmer : kmers) {
            int[] readIndexAndWithinReadIndex = GenomeSequenceEncodingUtil.decodeKmerBytesToIntArray(kmer.value);
            byte[] readBytes = this.reads.get(readIndexAndWithinReadIndex[0]);
            String read = GenomeSequenceEncodingUtil.getStringFromBPBytes(readBytes);
            String kmerString = read.substring(readIndexAndWithinReadIndex[1], readIndexAndWithinReadIndex[1] + k);

            int currentNodeIdx = 0;
            for (int charIdx = 0; charIdx < k; charIdx++) {
                char c = kmerString.charAt(charIdx);
                Map<Character, Integer> currentNode = trie.get(currentNodeIdx);
                if (currentNode.get(c) == null) {
                    // add it
                    Map<Character, Integer> newNode = new HashMap<Character, Integer>();
                    int newNodeIdx = trie.size();
                    trie.add(newNode);
                    currentNode.put(c, newNodeIdx);
                }
                currentNodeIdx = currentNode.get(c);
            }
        }

        System.out.println("Built a trie with " + trie.size() + " nodes");
    }

    /**
       Process each read and find if there is a single position we can update to get a large increase in kmer spectral alignment against the filtered kmer spectrum.

       Ensure buildTrie() is called before attempting to repair reads or a NotInitializedException will be thrown.

       @return Set of ReadCorrections specifying which bases in which reads to update to improve the kmer spectral alignment of the reads in this dataset.
     */
    public Set<ReadCorrection> recommendReadCorrections() {
        // for each read
        for (int readIdx = 0; readIdx < this.reads.size() - this.numFakeReadsAddedForOrphanedKmers; readIdx++) {
            ReadSpectralAlignmentResult result = checkRead(readIdx);

            // if the read is perfectly spectrally aligned, move on
            if (result.numSitesCausingMisalignedKmers == 0) {
                continue;
            }

            // determine whether this read just generally has very poor spectral alignment (if so we can mark it for being ignored)
            if (result.numMisalignedKmers > MISALIGNED_KMER_THRESHOLD) {
                spectrallyMisalignedReadsToIgnore.add(readIdx);
            }

            // determine whether there is a character we can replace which will correct kmer breaks for this read
            if (result.siteCausingMostBrokenKmers != null) {
                int charIdx = result.siteCausingMostBrokenKmers;
                //System.out.println("Read idx " + readIdx + ", charIdx " + charIdx + " is causing " + result.brokenKmersPerWithinReadIdx.get(charIdx) + " kmer breaks");
                Set<Character> potentiallyResolvingCharacters = result.resolvingCharactersByPosition.get(charIdx);
                if (potentiallyResolvingCharacters == null || potentiallyResolvingCharacters.size() == 0) {
                    // we could check the other sites causing broken kmers but if this read's main broken position is not resolvable let's move on
                    System.out.println("Read idx " + readIdx + ", charIdx " + charIdx + " is not resolvable, skipping");
                    continue;
                }

                byte[] oldReadBytes = this.reads.get(readIdx);
                String read = GenomeSequenceEncodingUtil.getStringFromBPBytes(oldReadBytes);
                StringBuilder sb = new StringBuilder(read);
                Character bestCharSubstitution = null;
                Integer bestNumMisalignedKmers = null;

                for (Character c : potentiallyResolvingCharacters) {
                    // temporarily fix the read
                    sb.setCharAt(charIdx, c);
                    String newRead = sb.toString();
                    byte[] newReadBytes = GenomeSequenceEncodingUtil.getBPBytesFromString(newRead);
                    this.reads.set(readIdx, newReadBytes);

                    ReadSpectralAlignmentResult newResult = checkRead(readIdx);
                    //System.out.println("previously read idx " + readIdx + ", charIdx " + charIdx + " had " + result.numMisalignedKmers + " now it has " + newResult.numMisalignedKmers + " misaligned kmers");
                    if (newResult.numMisalignedKmers < result.numMisalignedKmers - MIN_KMER_BREAKS_REQUIRED_CORRECTED_PER_UPDATE || (newResult.numMisalignedKmers == 0 && result.numMisalignedKmers != 0)) {
                        bestCharSubstitution = c;
                        bestNumMisalignedKmers = newResult.numMisalignedKmers;
                    } else {
                        //System.out.println("read idx " + readIdx + ", charIdx " + charIdx + " -> not enough of an improvement to make a correction");
                    }
                }

                this.reads.set(readIdx, oldReadBytes);

                if (bestCharSubstitution != null) {
                    readCorrectionsToMake.add(new ReadCorrection(readIdx, charIdx, bestCharSubstitution));

                    // we are moving this out to DeBruijnKmers.setupFromReadStringsAndSpecifiedK
                    // so that we can correct kmers as we update reads without reloading all of the reads
                    /**
                       sb = new StringBuilder(read);
                       sb.setCharAt(charIdx, bestCharSubstitution);
                       String newRead = sb.toString();
                       byte[] newReadBytes = GenomeSequenceEncodingUtil.getBPBytesFromString(newRead);
                       this.reads.set(readIdx, newReadBytes);
                    **/
                    //System.out.println("read idx " + readIdx + ", charIdx " + charIdx + " -> found a correction");
                }
            }
        }

        System.out.println("number spectrally misaligned reads (not amenable to point correction - generally low alignment): " + spectrallyMisalignedReadsToIgnore.size());
        System.out.println("number corrections recommended: " + readCorrectionsToMake.size());

        /**
         // we need to update our kmer references for this update if the kmer was first encountered on this particular read and is referenced elsewhere. we have not been doing these read removals at this time.
         int offset = 0;
         List<Integer> readsToRemove = new ArrayList<Integer>(spectrallyMisalignedReadsToIgnore);
         Collections.sort(readsToRemove);

         for (Integer misalignedReadToRemove : readsToRemove) {
         this.reads.remove(misalignedReadToRemove - offset++);
         }
        **/

        return readCorrectionsToMake;
    }

    /**
       Calling this method removes references to all the data structures used by this class to free up memory.
     */
    public void cleanup() {
        this.trie = null;
        this.kmers = null;
        this.reads = null;
        this.spectrallyMisalignedReadsToIgnore = null;
        this.readCorrectionsToMake = null;
        //System.gc();
    }

    /**
       Results of kmer spectral alignment for a read, used to recommend a ReadCorrection if the number of misaligned kmers (numMisalignedKmers) is high enough and there is a single site causing enough of those problems (see numSitesCausingMisalignedKmers, siteCausingMostBrokenKmers), and there is a single resolving character substituton for enough of those misaligned kmers (resolvingCharactersByPosition).
       @see org.wbdbga.KmerTrieFilterReads.ReadCorrection
       @see org.wbdbga.KmerTrieFilterReads#checkRead
       @see org.wbdbga.KmerTrieFilterReads#recommendReadCorrections
     */
    static final class ReadSpectralAlignmentResult {
        int readIdx;
        int numSitesCausingMisalignedKmers = 0;
        int numMisalignedKmers = 0;

        Integer siteCausingMostBrokenKmers = null;
        Map<Integer, Integer> brokenKmersPerWithinReadIdx;
        Map<Integer, Set<Character>> resolvingCharactersByPosition;

        /**
           Construct a ReadSpectralAlignmentResult by specifying the read index, the broken kmers by index within that read (at which position does kmer alignment fail), and the resoving characters for each position within the read to first order avoid kmer breaks.
           @param readIdx non-negative integer index into the reads in the input data specifying which read this spectral alignment result pertains to
           @param brokenKmersPerWithinReadIdx a map of the number of kmer alignments that failed by index within this read (offset)
           @param resolvingCharactersByPosition a map by the offset within this read to the characters that would have first-order fixed the kmer spectral alignments that failed at that position
           @see org.wbdbga.KmerTrieFilterReads#checkRead
        */
        public ReadSpectralAlignmentResult(int readIdx, Map<Integer, Integer> brokenKmersPerWithinReadIdx, Map<Integer, Set<Character>> resolvingCharactersByPosition) {
            if (readIdx < 0) {
                throw new IllegalArgumentException("readIdx argument must be non-negative");
            }
            if (brokenKmersPerWithinReadIdx == null) {
                throw new IllegalArgumentException("brokenKmersPerWithinReadIdx argument must be non-null");
            }
            if (resolvingCharactersByPosition == null) {
                throw new IllegalArgumentException("resolvingCharactersByPosition argument must be non-null");
            }
            this.readIdx = readIdx;
            this.brokenKmersPerWithinReadIdx = brokenKmersPerWithinReadIdx;
            this.resolvingCharactersByPosition = resolvingCharactersByPosition;
        }
    }

    /**
       Check a read against the set of solid kmers for alignment, returning a ReadSpectralAlignmentResult, including information that recommendReadCorrections can use to recommend a ReadCorrection if the number of misaligned kmers (numMisalignedKmers) is high enough and there is a single site causing enough of those problems (see numSitesCausingMisalignedKmers, siteCausingMostBrokenKmers), and there is a single resolving character substituton for enough of those misaligned kmers (resolvingCharactersByPosition).

       @see org.wbdbga.KmerTrieFilterReads#buildTrie
       @see org.wbdbga.KmerTrieFilterReads#recommendReadCorrections
       @see org.wbdbga.KmerTrieFilterReads.ReadCorrection

       @param readIdx the index of the read to check for kmer spectral alignment and return a ReadSpectralAlignmentResult for.
       @return ReadSpectralAlignmentResult with including information used to recommend a ReadCorrection like numMisalignedKmers, numSitesCausingMisalignedKmers, siteCausingMostBrokenKmers, and resolvingCharactersByPosition
     */
    public ReadSpectralAlignmentResult checkRead(int readIdx) {
        if (this.trie == null || this.trie.size() == 0) {
            throw new NotInitializedException("ensure buildTrie is called before attempting to check or repair reads");
        }
        if (readIdx < 0) {
            throw new IllegalArgumentException("readIdx must be non-negative");
        }
        if (readIdx >= this.reads.size()) {
            throw new IllegalArgumentException("readIdx must be within the bounds of this.reads");
        }

        // initialize a map of index to count of broken kmers
        Map<Integer, Integer> brokenKmersPerWithinReadIdx = new HashMap<Integer, Integer>();
        // initialize a map of index to characters which would avoid broken kmers (first-order, not guaranteed to match all the way down the rest of the string yet)
        Map<Integer, Set<Character>> resolvingCharactersByPosition = new HashMap<Integer, Set<Character>>();
        int numMisalignedKmers = 0;
        Integer siteCausingMostBrokenKmers = null;

        // populate those maps by walking along the read and trying to match the kmer trie from each position
        String read = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(readIdx));
        for (int withinReadIdx = 0; withinReadIdx < read.length() - (k - 1); withinReadIdx++) {
            int currentNodeIdx = 0;
            Map<Character, Integer> trieNode; // = trie.get(currentNodeIdx);
            for (int wrIdx2 = 0; wrIdx2 < k; wrIdx2++) {
                trieNode = trie.get(currentNodeIdx);
                int charIdx = withinReadIdx + wrIdx2;
                char c = read.charAt(charIdx);
                Integer nextNodeIdx = trieNode.get(c);

                // we have a kmer break
                if (nextNodeIdx == null) {
                    Integer numBrokenKmers = brokenKmersPerWithinReadIdx.get(charIdx);
                    Set<Character> allowableNextCharacters = trieNode.keySet();
                    if (numBrokenKmers == null) { 
                        numBrokenKmers = 0;                        
                    }

                    brokenKmersPerWithinReadIdx.put(charIdx, ++numBrokenKmers);
                    if (siteCausingMostBrokenKmers == null || numBrokenKmers.compareTo(brokenKmersPerWithinReadIdx.get(siteCausingMostBrokenKmers)) > 0) {
                        siteCausingMostBrokenKmers = charIdx;
                    }
                    Set<Character> rcbp = resolvingCharactersByPosition.get(charIdx);
                    if (rcbp == null) {
                        rcbp = new HashSet<Character>();
                        rcbp.addAll(allowableNextCharacters);
                        resolvingCharactersByPosition.put(charIdx, rcbp);
                    } else {
                        // take intersection with existing
                        for (char pc : GenomeSequenceEncodingUtil.ALPHABET) {
                            if (!allowableNextCharacters.contains(pc) && rcbp.contains(pc)) {
                                rcbp.remove(pc);
                            }
                        }
                    }
                    numMisalignedKmers++;
                    break;
                } else {
                    // we don't have a kmer break, continue the kmer
                    currentNodeIdx = nextNodeIdx;
                }
            }
        }

        ReadSpectralAlignmentResult result = new ReadSpectralAlignmentResult(readIdx, brokenKmersPerWithinReadIdx, resolvingCharactersByPosition);
        result.numSitesCausingMisalignedKmers = brokenKmersPerWithinReadIdx.keySet().size();
        result.numMisalignedKmers = numMisalignedKmers;
        result.siteCausingMostBrokenKmers = siteCausingMostBrokenKmers;
        return result;
    }
    static class NotInitializedException extends RuntimeException {
        public NotInitializedException(String message) {
            super(message);
        }
    }
}
