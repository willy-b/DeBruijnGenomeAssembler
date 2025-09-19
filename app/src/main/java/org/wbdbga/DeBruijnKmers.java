/*
DeBruijnKmers for DeBruijnGenomeAssembler - an unpaired short read (~100bp) de Bruijn graph based genome assembler
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
Critical class that builds a de Bruijn map from the raw reads and a specified k value and also performs cleaning of the reads using kmer spectral alignment and drops low coverage kmers - the output of this class is passed to the constructor of DeBruijnGraph.

Started with an approach similar to what is discussed in UCSD ALGS206x,
and discussed in chapter 3 of the textbook: Compeau, Phillip and Pevzner, Pavel. Bioinformatics algorithms : an active learning approach,3rd ed. Active Learning Publishers, 2018, available at https://www.bioinformaticsalgorithms.org/ .
Then updated to support read coverage tracking and cutoff, then updated to not explicitly represent kmer strings, use int pairs,
and finally for kmer spectral alignment (see also KmerTrieFilterReads.java for that).

@see org.wbdbga.DeBruijnKmers#setupFromReadStringsAndSpecifiedK
@see org.wbdbga.DeBruijnKmers#processKmers
@see org.wbdbga.DeBruijnGraph
@see org.wbdbga.KmerTrieFilterReads
*/
public class DeBruijnKmers {
    static final boolean VERBOSE = true;

    List<byte[]> reads = null; // Shared reference to reads represented once upstream. Should not be modified. For unpaired data, size is number of unpaired reads. For paired data (not implemented in this version), size would be 2 * number of pairs of reads.
    Map<Integer, List<ValueComparedByteArray>> kmerHashCodeToReadIndexAndSubstringIndexMap = new HashMap<Integer, List<ValueComparedByteArray>>();
    Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> debruijnMap = new HashMap<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>>();

    // We are trying to avoid representing millions of redundant substrings from the reads in the HashMap keysets and lists and elsewhere, which is not feasible (only works to explicitly represent kmers up to ~500K bp genomes with a few GB of ram)
    /**
       Length k of the kmers (a kmer is a substring of length k shorter than the read length for which we expect to have coverage in our samples of all unique values present in the genome being sequenced/assembled). In the de Bruijn graph the edges will represent the kmers of the reads and the nodes represent (k-1)mers.
     */
    public int k;

    /**
       @deprecated
       boolean flag indicating direct evidence that kmers are not unique by seeing them repeated within individual reads has been observed, used to recommend increasing k. Deprecated because the tangling of the de Bruijn graph caused by this can be observed directly downstream (sometimes even very non-unique kmers can be handled by the de Bruijn graph as we know coverage ratios) and recommend increasing k from there. May bring this back but tracking it was causing a lot of memory churn and is unnecessary.
    */
    @Deprecated
    boolean kmersNotUniqueInIndividualReads = false;

    /**
       Track the number of times each kmer has been seen across reads.

       @see org.wbdbga.DeBruijnKmers#loadRead
       @see org.wbdbga.DeBruijnKmers#unloadRead
     */
    Map<ValueComparedByteArray, Byte> kmerCoverageMap = null;

    /**
       The length in bases of a read in the input data. Fixed length is required in this version but variable length reads will be supported.
     */
    int readLen;

    /**
       total number of observed kmers (even counting repeats for the same kmers).
     */
    int totalReadCoverageForKmers = 0;

    /**
       average read coverage (number of times observed) for kmers.
    */
    Double averageReadCoverageForKmers = null;

    /**
       Since we handle kmers as references into a read in which they were observed, but correct reads, it is possible for the original read that a kmer is stored as a reference into to no longer contain that kmer, in which case we introduce a fake read entry used only to contain the orphaned kmer and otherwise excluded from analysis. We count the number of such fake reads added to the end so we can stop early when iterating through the reads and avoid them in most places where we are only interested in real reads (data).
     */
    int numFakeReadsAddedForOrphanedKmers = 0;

    /**
       This setting makes assembly a few fold slower but mostly avoids the appearance of crashing one's computer for up to twenty minutes at a time on megabase genomes with decent read coverage while this runs (or potentially hours if searching multiple k values is enabled). When this is true, instead for an equivalent load previously making the computer unresponsive for twenty minutes (will vary with machine specifications and genome size), it might run for over an hour per k value, but without interrupting use of the computer completely. Set to false if you just want it to go quickly and don't need to interact with the machine in the meantime.
     */
    public boolean PRIORITIZE_NOT_BLOCKING_OS_OVER_SPEED = true;

    /**
       Large genomes with distinct kmer count greater than this are treated differently due to assumptions of lower coverage and potential problems keeping them in memory.
     */
    final int LARGE_GENOME_KMER_THRESHOLD = 5000000;

    // two flags to control during unit testing
    // sometimes we need to disable trie corrections to allow errors through when setting up unit tests for other parts of the code
    boolean disableTrieCorrections = false;
    // in unit tests we may want to use load/unload read after setupFromReadStringsAndSpecifiedK has completed (unrealistic, normally we free up the memory)
    boolean freeKmerHashCodeToReadIndexAndSubstringIndexMapAtEndOfSetup = true;

    // for unit tests where we want trie corrections enabled
    boolean useTrieCorrectionsForSmallGenomes = false;

    // TODO consider to combine the constructor with setupFromReadStringsAndSpecifiedK
    /**
       Construct a DeBruijnKmers, after which it can be set up using setupFromReadStringsAndSpecifiedK.
    */
    public DeBruijnKmers() {}

    /**
       Setup DeBruijnKmers with a k value and the input data in the form of reads (encoded using GenomeSequenceEncodingUtil.getBPBytesFromString).
       @see org.wbdbga.GenomeSequenceEncodingUtil#getBPBytesFromString
       @param k Length k of the kmers (a kmer is a substring of length k shorter than the read length for which we expect to have coverage in our samples of all unique values present in the genome being sequenced/assembled). In the de Bruijn graph the edges will represent the kmers of the reads and the nodes represent (k-1)mers.
       @param reads the input data in the form of short contiguous sequencing reads which potentially contain errors to be corrected, as a list of byte[] encoded using GenomeSequenceEncodingUtil.getBPBytesFromString
    */
    public void setupFromReadStringsAndSpecifiedK(int k, List<byte[]> reads) {
        if (reads == null) {
            throw new IllegalArgumentException("reads list must not be null");
        }

        if (reads.size() == 0) {
            throw new IllegalArgumentException("reads list must not be empty");
        }
        if (reads instanceof LinkedList) {
            reads = new ArrayList<byte[]>(reads);
        }

        this.reads = reads;

        String firstRead = GenomeSequenceEncodingUtil.getStringFromBPBytes(reads.get(0));

        if (k > firstRead.length()) {
            throw new IllegalArgumentException("k (==" + k + ") must be <= read length (== " + firstRead.length() + ")");
        }
        this.k = k;

        boolean noTrieCorrectionsAchieved = false;

        readLen = firstRead.length();
        kmerCoverageMap = new HashMap<ValueComparedByteArray, Byte>(this.reads.size() * readLen);
        totalReadCoverageForKmers = 0;
        averageReadCoverageForKmers = null;

        loadKmerCoverageMapFromReads();

        int numDistinctKmers = kmerCoverageMap.keySet().size();

        if (VERBOSE) {
            System.out.println("Generated kmers from reads for k == " + k + ", processing...");
        }

        /*
           estimate size of genome before doing spectral correction
           for decent sized k, this is approximately the number of kmers left after removing low coverage kmers
        */
        Integer estimatedGenomeSizeBp = null; // we will set this after removing low coverage kmers

        /*
           estimate from previous the read coverage, classify has low or high read coverage
           for low read coverage or estimated very large genomes, do kmer spectral alignment
           for high read coverage or estimated small genomes, skip kmer spectral alignment and handle it in the de bruijn graph

           this will save approximately 45 seconds on n. deltocephalinicola and carsonella ruddii sized genomes and 20 seconds on phi x174 where this step is just unnecessary.
        */
        Double estimatedAverageCoverageOfGenome = null;

        boolean largeGenomeKmerThresholdReached = false;
        averageReadCoverageForKmers = ((double)totalReadCoverageForKmers)/((double)numDistinctKmers);
        int numRoundsSpectralCorrection = 0;

        if (numDistinctKmers > LARGE_GENOME_KMER_THRESHOLD) {
            largeGenomeKmerThresholdReached = true;
        }

        while (!noTrieCorrectionsAchieved) {
            boolean removedKmers = true;
            while (removedKmers) {
                removedKmers = false;
                int removedCount = removeLowCoverageKmers();
                if (removedCount > 0) { removedKmers = true; }
                //System.gc();
                if (VERBOSE) {
                    System.out.println("Removed " + removedCount + " extremely low coverage kmers (kmerCoverageMap.keySet().size() == " + kmerCoverageMap.keySet().size() + ")");
                }
            }

            estimatedGenomeSizeBp = kmerCoverageMap.keySet().size();
            estimatedAverageCoverageOfGenome = ((double)totalReadCoverageForKmers)/((double)estimatedGenomeSizeBp);

            if (VERBOSE) {
                System.out.println("estimatedGenomeSizeBp: " + estimatedGenomeSizeBp);
                System.out.println("estimatedAverageCoverageOfGenome: " + estimatedAverageCoverageOfGenome);
            }

            if (disableTrieCorrections || (!useTrieCorrectionsForSmallGenomes && (estimatedGenomeSizeBp < 200000 || estimatedAverageCoverageOfGenome > 20))) {
                break;
            }

            int kToUseForTrie = k;
            if (kmerCoverageMap.keySet().size() > 1000000) {
                kToUseForTrie = Math.min(k / 2, 8);
            }

            if (VERBOSE) {
                System.out.println("Building kmer trie (to do iterative spectral alignment of the reads) using k==" + kToUseForTrie);
            }

            KmerTrieFilterReads ktfr = new KmerTrieFilterReads(reads, kmerCoverageMap.keySet(), kToUseForTrie, numFakeReadsAddedForOrphanedKmers);
            ktfr.buildTrie();

            Set<KmerTrieFilterReads.ReadCorrection> s = ktfr.recommendReadCorrections();

            if (s == null || s.size() == 0) {
                noTrieCorrectionsAchieved = true;
            } else {
                // for each read correction (at most 1 per read in current implementation)
                // we want to apply it in place (updating kmers) without reloading all the reads.
                for (KmerTrieFilterReads.ReadCorrection rc : s) {
                    totalReadCoverageForKmers += applyReadCorrection(rc);
                }
            }
            numRoundsSpectralCorrection++;

            System.out.println("Performed " + s.size() + " corrections to reads. Reprocessing.");
            s = null;
            ktfr.cleanup();
            ktfr = null;

            recomputeTotalAndAverageKmerReadCoverage();
        }

        if (!largeGenomeKmerThresholdReached) {
            loadKmerCoverageMapFromReads(); // we do this only for small genomes where we can afford to bring back in all low coverage kmers not spectrally corrected.
        }

        recomputeTotalAndAverageKmerReadCoverage();

        if (freeKmerHashCodeToReadIndexAndSubstringIndexMapAtEndOfSetup) {
            // free up memory, once read coverage is computed we don't need it any longer
            kmerHashCodeToReadIndexAndSubstringIndexMap = null;
            kmerHashCodeToReadIndexAndSubstringIndexMap = new HashMap<Integer, List<ValueComparedByteArray>>();
        }
    }

    private void recomputeTotalAndAverageKmerReadCoverage() {
        totalReadCoverageForKmers = 0;

        for (ValueComparedByteArray kmer : kmerCoverageMap.keySet()) {
            totalReadCoverageForKmers += (int)kmerCoverageMap.get(kmer);
        }

        int numDistinctKmers = kmerCoverageMap.keySet().size();
        averageReadCoverageForKmers = ((double)totalReadCoverageForKmers)/((double)numDistinctKmers);
    }

    private int removeLowCoverageKmers() {
        return removeLowCoverageKmers(null);
    }

    private int removeLowCoverageKmers(Integer overrideThreshold) {
        int removedCount = 0;
        int numToRemove = 0;
        List<ValueComparedByteArray> kmersToRemove = new LinkedList<ValueComparedByteArray>();

        boolean largeGenomeKmerThresholdReached = false;
        if (kmerCoverageMap.keySet().size() > LARGE_GENOME_KMER_THRESHOLD) {
            largeGenomeKmerThresholdReached = true;
        }

        for (ValueComparedByteArray kmer : kmerCoverageMap.keySet()) {
            // TODO make the minimum threshold here configurable as this may be something that prevents certain datasets from being assembled due to being too aggressive
            // someone who is not getting decent sized contigs should be able to turn this off / make this less aggressive to see if that fixes their problem.
            int threshold = (int)Math.min(averageReadCoverageForKmers / 4, 6); // was experimenting with being very aggressive on dropping kmers before the final round (when the reduced kmer pool is used for kmer spectral alignment). however it caused a disconnected de bruijn graph (11K nodes out of 5,062,927 nodes), per my current threshold. I ran it with this aggressive mode again and raised the imbalance node threshold and checked the contigs. There was not an improvement there, so I set this threshold lower and for large genomes just use threshold 1.
            if (overrideThreshold != null) {
                threshold = overrideThreshold;
            }

            if (largeGenomeKmerThresholdReached) {
                // unfortunately in this case we cannot afford to add back in these low coverage nodes that we could not correct via kmer spectrum and correct the related issues in the de bruijn graph, it will exceed our memory budget. so we still need to do some filtering.
                // however we will be conservative and only drop low coverage on an absolute scale.

                // for the e coli genome, I am seeing that since we cannot keep all the kmers that were not spectrally aligned, doing a mix of approaches is worse than just using a low threshold and dropping the kmers only observed once (almost always an error as we have > 1x coverage).
                threshold = 1;
            }

            if (kmerCoverageMap.get(kmer) <= Byte.valueOf((byte)Math.max(1, threshold))) {
                kmersToRemove.add(kmer);
                numToRemove++;
            }
            if (numToRemove >= 500000 && kmerCoverageMap.keySet().size() > 1000000) {
                break; // otherwise we are getting too close to the memory limit
            }
        }

        for (ValueComparedByteArray kmer : kmersToRemove) {
            if (kmerCoverageMap.containsKey(kmer)) {
                kmerCoverageMap.remove(kmer);
                // also remove reference to this kmer byte array from the string-keyed data structure
                int[] readIndexAndOffset = GenomeSequenceEncodingUtil.getReadIndexAndWithinReadIndexFrom4ByteRepresentation(kmer.value);
                String kmerString = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(readIndexAndOffset[0])).substring(readIndexAndOffset[1], readIndexAndOffset[1] + k);
                List<ValueComparedByteArray> savedKmers = kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmerString.hashCode());
                if (savedKmers != null) {
                    if (savedKmers.contains(kmer)) {
                        savedKmers.remove(kmer);
                    }
                    if (savedKmers.size() == 0) {
                        kmerHashCodeToReadIndexAndSubstringIndexMap.remove(kmerString.hashCode());
                    }
                }
                // end cleanup in other data structure
                removedCount++;
            }
        }

        kmersToRemove = null;
        return removedCount;
    }

    private void loadKmerCoverageMapFromReads() {
        kmerCoverageMap.clear();
        if (VERBOSE) {
            System.out.println("started loadKmerCoverageMapFromReads");
        }

        // TODO assert expected bounds on numFakeReadsAddedForOrphanedKmers
        int upperBound = this.reads.size() - numFakeReadsAddedForOrphanedKmers;

        for (int readIdx = 0; readIdx < upperBound; readIdx++) {
            totalReadCoverageForKmers += loadRead(readIdx);

            if (readIdx % 1000 == 0) {
                if (VERBOSE) {
                    System.out.println("Processed read " + (readIdx+1) + " of " + this.reads.size() + " for k == " + k + " keyset: " + kmerHashCodeToReadIndexAndSubstringIndexMap.keySet().size());
                }

                // makes it a few fold slower but mostly avoids the appearance of crashing the end user's computer for up to twenty minutes on megabase genomes with decent read coverage while this runs (or potentially hours if searching multiple k values is enabled) (instead e.g. runs for over an hour, instead of twenty minutes, per k value but without interrupting use of the computer completely)
                if (PRIORITIZE_NOT_BLOCKING_OS_OVER_SPEED) {
                    // try to avoid long blocking GC pauses
                    try {
                        Thread.sleep(1);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    if (readIdx % 15000 == 0) {
                        System.gc();
                    }
                }
            }
        }
    }

    /**
       Decrement read coverage for kmers in read and or delete the kmers. If this read is the defining read for that kmer and read coverage is greater than 1, create a new fake read containing just that kmer to refer to.
       @param readIdx the index of the read to unload
       @return total read coverage removed across kmers
     */
    int unloadRead(int readIdx) {
        int totalReadCoverageRemoved = 0;
        String read = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(readIdx));

        // the first steps are the opposite of loadRead
        int readLen = read.length();

        // same loop structure as in loadRead(readIdx)
        for (int offset = 0; offset < readLen - k + 1; offset++) {
            String kmer;

            kmer = read.substring(offset, offset + k);
            byte[] kmerReadAndIndexPair = getKmerReadAndIndexPairForKmerString(kmer);
            int[] kmerReadAndIndexPairInts;

            boolean wasInTable = kmerReadAndIndexPair != null;

            if (!wasInTable) {
                // This is actually expected behavior due to removeLowCoverageKmers . Some kmers were already removed due to not being part of the solid kmer spectrum.
                // Thus when we go to unload them from reads, they are already gone.

                kmerReadAndIndexPairInts = new int[]{readIdx, offset};
                kmerReadAndIndexPair = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(kmerReadAndIndexPairInts[0], kmerReadAndIndexPairInts[1]);
            } else {
                kmerReadAndIndexPairInts = GenomeSequenceEncodingUtil.decodeKmerBytesToIntArray(kmerReadAndIndexPair);
            }

            int count = ((int)kmerCoverageMap.getOrDefault(ValueComparedByteArray.from(kmerReadAndIndexPair), (byte)0)) - 1;
            totalReadCoverageRemoved++;

            if (count <= 0) {
                // The negative count case is to be expected if the kmer was dropped intentionally due to low coverage in removeLowCoverageKmers step.
                // Kmer spectral filtering is performed after dropping low coverage reads to align with the solid kmers
                // (alternatively we could just ignore them when building the kmer trie for spectral alignment below the "solid" threshold, without removing them in the map but we do want them out of the de Bruijn graph as well).
                
                // This is the easy case. We can remove the kmer.
                List<ValueComparedByteArray> savedKmers = kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmer.hashCode());

                if (savedKmers != null) {
                    if (savedKmers.contains(ValueComparedByteArray.from(kmerReadAndIndexPair))) {
                        savedKmers.remove(ValueComparedByteArray.from(kmerReadAndIndexPair));
                    }
                    if (savedKmers.size() == 0) {
                        kmerHashCodeToReadIndexAndSubstringIndexMap.remove(kmer.hashCode());
                    }
                }
                kmerCoverageMap.remove(ValueComparedByteArray.from(kmerReadAndIndexPair));
            } else if (kmerReadAndIndexPairInts[0] == readIdx) {
                // it is OK to do nothing for count > 0, so long as this read is not the defining read for this kmer
                // if THIS IS the defining read in this kmer AND coverage is still > 0 then we need define a new fake read and update the kmer to point to that
                int newKmerFakeReadIdx = this.reads.size();
                byte[] newKmerReadAndIndexPair = null;
                this.reads.add(GenomeSequenceEncodingUtil.getBPBytesFromString(kmer));

                numFakeReadsAddedForOrphanedKmers++;

                newKmerReadAndIndexPair = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(newKmerFakeReadIdx, 0);

                // update kmer coverage map with new key
                kmerCoverageMap.remove(ValueComparedByteArray.from(kmerReadAndIndexPair));
                kmerCoverageMap.put(ValueComparedByteArray.from(newKmerReadAndIndexPair), (byte)count);

                // update kmerHashCodeToReadIndexAndSubstringIndexMap
                List<ValueComparedByteArray> savedKmers = kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmer.hashCode());
                if (savedKmers == null) {
                    savedKmers = new ArrayList<ValueComparedByteArray>();
                    kmerHashCodeToReadIndexAndSubstringIndexMap.put(kmer.hashCode(), savedKmers);
                } else {
                    savedKmers.remove(ValueComparedByteArray.from(kmerReadAndIndexPair));
                }
                ValueComparedByteArray newVCBA = ValueComparedByteArray.from(newKmerReadAndIndexPair);
                // it should not ever be the case that the new read and index pair is already there,
                // (because only one read and index pair defines a kmer at a time)
                // but if it is, we don't need to add it
                if (!savedKmers.contains(newVCBA)) {
                    savedKmers.add(newVCBA);
                }

            } else { // just update the count
                kmerCoverageMap.put(ValueComparedByteArray.from(kmerReadAndIndexPair), (byte)count);
            }
        }

        return totalReadCoverageRemoved;
    }

    // returns net change to overall read coverage across kmers (dividing this tracked overall number by number of kmers gives average)
    public int applyReadCorrection(KmerTrieFilterReads.ReadCorrection rc) {
        // unload the existing read (update read coverage for all kmers in it)
        int readCoverageRemoved = unloadRead(rc.readIdx);
        //System.out.println("Read coverage removed: " + readCoverageRemoved);

        // edit the read to apply the correction
        String read = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(rc.readIdx));
        StringBuilder sb = new StringBuilder(read);
        sb.setCharAt(rc.withinReadIdx, rc.newCharacter);
        String newRead = sb.toString();
        byte[] newReadBytes = GenomeSequenceEncodingUtil.getBPBytesFromString(newRead);
        this.reads.set(rc.readIdx, newReadBytes);

        // reload the read
        int readCoverageAdded = loadRead(rc.readIdx);
        //System.out.println("Read coverage added: " + readCoverageAdded);

        return readCoverageAdded - readCoverageRemoved;
    }

    /**
       Increment read coverage for kmers in read and / or load new kmers defined by it.
       @param readIdx the index of the read to load
       @return total read coverage added across kmers
     */
    int loadRead(int readIdx) throws NotInitializedException {
        int totalReadCoverageAdded = 0;
        String read;

        try {
            read = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(readIdx));
        } catch (Exception e) {
            System.err.println("Was not able to process a read. 0-based line index: " + readIdx + ". continuing.");
            e.printStackTrace();
            return 0;
        }

        int readLen = read.length();

        //Set<String> kmersInRead = new HashSet<String>(read.length());

        for (int offset = 0; offset < readLen - k + 1; offset++) {
            String kmer;
            kmer = read.substring(offset, offset + k);
            /*
            if (kmersInRead.contains(kmer)) {
               kmersNotUniqueInIndividualReads = true;
               System.out.println("kmer \"" + kmer + "\" occurs more than once within read \"" + read + "\", this is nonideal and suggests k should be increased, but can be natural for highly repetitive genomes, which may be better assembled using paired reads (not supported in this version)");
            }
            kmersInRead.add(kmer);
            */

            byte[] kmerReadAndIndexPair;
            //System.out.println("Called putKmerReadAndIndexPairForKmerString with " + kmer + ")");
            kmerReadAndIndexPair = putKmerReadAndIndexPairForKmerString(kmer, readIdx, offset);

            int[] kmerReadAndIndexPairInt = GenomeSequenceEncodingUtil.decodeKmerBytesToIntArray(kmerReadAndIndexPair);

            int count = ((int)kmerCoverageMap.getOrDefault(ValueComparedByteArray.from(kmerReadAndIndexPair), (byte)0)) + 1;

            totalReadCoverageAdded++;

            kmerCoverageMap.put(ValueComparedByteArray.from(kmerReadAndIndexPair), (byte)count);
        }

        return totalReadCoverageAdded;
    }

    int largestNumberOfValuesPerKey = 0;
    byte[] putKmerReadAndIndexPairForKmerString(String kmer, int readIndex, int offset) throws NotInitializedException {
        String originalKmer = kmer;

        List<ValueComparedByteArray> entries = kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmer.hashCode());
        if (entries == null) {
            entries = new LinkedList<ValueComparedByteArray>();
            kmerHashCodeToReadIndexAndSubstringIndexMap.put(kmer.hashCode(), entries);
        }

        byte[] readAndIndexPair2 = getKmerReadAndIndexPairForKmerString(originalKmer);
        byte[] readAndIndexPair;
        if (readAndIndexPair2 == null) {
            //readAndIndexPair = new Integer[]{readIndex, offset};
            if (this.reads == null) {
                throw new NotInitializedException("putKmerReadAndIndexPairForKmerString has been called on a DeBruijnKmers instance that has not been properly initialized; please ensure setupFromReadStringsAndSpecifiedK is called first and that DeBruijnKmers.reads is not being set to null");
            }
            if (readIndex >= this.reads.size()) {
                throw new OutOfRangeReadIndexException("putKmerReadAndIndexPairForKmerString has been called with an invalid entry referring to a read index that is outside the range of this.reads, this should not ever happen.");
            }
            // TODO: this validation step is slow
            if (offset < 0 || offset + kmer.length() > GenomeSequenceEncodingUtil.getDeclaredBPLength(this.reads.get(readIndex))) {
                throw new OutOfRangeReadOffsetException("putKmerReadAndIndexPairForKmerString has been called with an invalid entry referring to a offset that denotes a kmer extending beyond the range of the read, this should not ever happen.");
            }
            readAndIndexPair = GenomeSequenceEncodingUtil.get4ByteRepresentationForReadIndexAndWithinReadIndex(readIndex, offset);
            //System.out.println("For kmer=\"" + kmer + "\" set readIndex=" + readIndex + ", within read index=" + offset + ", and hashCode=" + kmer.hashCode());
            entries.add(ValueComparedByteArray.from(readAndIndexPair));
            //System.out.println("Added to kmer \""+ kmer +"\" with hashCode=" + kmer.hashCode() + ", entries: " + entries);
            if (entries.size() > largestNumberOfValuesPerKey) {
                largestNumberOfValuesPerKey = entries.size();
                if (VERBOSE) {
                    System.out.println("kmer \"" + kmer + "\" shares a hashCode with " + (largestNumberOfValuesPerKey - 1) + " other kmers");
                }
            }
        } else {
            // let me note it is intended behavior that if this kmer is aleady stored,
            // the readIndex and offset arguments are just ignored and we use the existing values
            // that is because for this representation we intend to store the value of the kmer string as a pointer into some consistent instance of it somewhere
            // not all instances, to avoid string comparisons by allowing direct comparisons of the byte pointer representations
            readAndIndexPair = readAndIndexPair2;
        }

        return readAndIndexPair;
    }

    byte[] getKmerReadAndIndexPairForKmerString(String kmer) throws NotInitializedException {
        byte[] kmerReadAndIndexPair = null;
        // yes we look up the hashcode instead of string. this is because there are millions (for >1M bp genomes) of highly overlapping kmer strings and we cannot afford to have our hashmap keeping these strings in the keyset when we can keep indices instead.
        // but multiple indices can match the same kmer and we want to only return the earliest result.

        //System.out.println("Retrieve from kmerHashCodeToReadIndexAndSubstringIndexMap for kmer " + kmer + " with hashCode " + kmer.hashCode());
        List<ValueComparedByteArray> entries = kmerHashCodeToReadIndexAndSubstringIndexMap.get(kmer.hashCode());
        //System.out.println(entries);
        if (entries != null) {
            for (ValueComparedByteArray entryBytes : entries) {
                int[] entry = GenomeSequenceEncodingUtil.decodeKmerBytesToIntArray(entryBytes.value);

                String otherKmer;
                if (this.reads == null) {
                    throw new NotInitializedException("getKmerReadAndIndexPairForKmerString has been called on a DeBruijnKmers instance that has not been properly initialized; please ensure setupFromReadStringsAndSpecifiedK is called first and that DeBruijnKmers.reads is not being set to null");
                }

                if (entry[0] >= this.reads.size()) {
                    throw new CorruptedReadIndexAndSubstringIndexMapException("getKmerReadAndIndexPairForKmerString found an invalid entry referring to a read index that is outside the range of this.reads, this should not ever happen.");
                }

                String otherRead = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(entry[0]));
                if (entry[1] < 0 || entry[1] + kmer.length() > otherRead.length()) {
                    throw new CorruptedReadIndexAndSubstringIndexMapException("getKmerReadAndIndexPairForKmerString found an invalid entry, with an offset which overflows the read length, this should not ever happen.");
                }

                otherKmer = otherRead.substring(entry[1], entry[1] + kmer.length()); // for a value of k, say, 10, in the initial processing we deal with k==10 length strings. in DeBruijnGraph they get split into 2 debruijn graph nodes each of length (k-1). so this preprocessing should be dealing with length k whereas in DeBruijnGraph we need to make sure it is length k-1.

                //System.out.println("compare \"" + kmer + "\" with \"" + otherKmer + "\"");
                if (kmer.equals(otherKmer)) {
                    kmerReadAndIndexPair = entryBytes.value;
                    break;
                } else {
                    //if (VERBOSE) { System.out.println("Compared two kmers with same hash code which were found not to be equal: \"" + kmer + "\" and \"" + otherKmer + "\""); }
                }
            }
        }

        /*
         // done in calling method
         if (kmerReadAndIndexPair == null) {
         kmerReadAndIndexPair = new Integer[]{readIdx, offset};
         }
        */

        return kmerReadAndIndexPair;
    }

    /**
    Splits kmers to prefix and suffix (k-1)-mers for each kmer with counts of transitions from each prefix to suffix.
    This is necessary because then the edges connecting the prefix and suffix (k-1)-mers represent the original kmer,
    which means the read/fragment assembly problem is finding a path that visits all edges in a graph once, NOT visits all nodes in a graph once.
    Finding the longest path visiting every node once (a Hamiltonian path) is NP-complete, whereas finding a path that visits every edge in a graph once (an Eulerian path) can be solved in linear time.
    See "How to apply de Bruijn graphs to genome assembly" by Compeau PE, Pevzner PA, Tesler G at https://pmc.ncbi.nlm.nih.gov/articles/PMC5531759 or "Bioinformatics algorithms : an active learning approach, 3rd edition" by Phillip Compeau and Pavel Pevzner (2018) from Active Learning Publishers, available at https://www.bioinformaticsalgorithms.org/ for more details.
    The prefix to suffix (k-1)-mer to (k-1)-mer counts map created here in debruijnMap can then be retrieved with DeBruijnKmers.getDeBruijnMap() and used to construct a DeBruijnGraph instance to solve the problem.
    @see org.wbdbga.DeBruijnKmers#getDeBruijnMap
    @see org.wbdbga.DeBruijnGraph
    @throws NotInitializedException if setupFromReadStringsAndSpecifiedK has not been called first to populate the kmerCoverageMap
    */
    public void processKmers() {
        if (kmerCoverageMap == null || kmerCoverageMap.keySet().size() == 0) {
            throw new NotInitializedException("kmerCoverageMap is not yet loaded. setupFromReadStringsAndSpecifiedK must be called before this method.");
        }

        int numDistinctKmers = kmerCoverageMap.keySet().size();

        kmerHashCodeToReadIndexAndSubstringIndexMap = null; // no longer needed! (we will replace it with (k-1)-mers in the map as we break kmers observed into nodes for prefix and suffix)
        //System.gc();
        kmerHashCodeToReadIndexAndSubstringIndexMap = new HashMap<Integer, List<ValueComparedByteArray>>();

        // similar to what I did after splitting the text into K-mers in DeBruijnString.java ("4. Construct the de Bruijn Graph of a String")
        // except here we do not handle kmer strings explicitly (redundant; requires too much memory) so bit more complicated to handle indices only (and replace indices with first place it occurs so that indices are unique even if the kmer occurs across many reads).
        int numKmersChecked = 0;
        int numKmersToCheck = kmerCoverageMap.keySet().size();

        Queue<Map.Entry<ValueComparedByteArray, Byte>> kmersToProcess = new LinkedList<Map.Entry<ValueComparedByteArray, Byte>>(kmerCoverageMap.entrySet());
        kmerCoverageMap = null;
        //System.gc();

        if (VERBOSE) {
            System.out.println("Starting to generate De Bruijn map [(k-1)mer to (k-1)mer edges] from kmers");
        }

        while (!kmersToProcess.isEmpty()) {
            Map.Entry<ValueComparedByteArray, Byte> kmerWithCount = kmersToProcess.poll();
            byte[] kmerBytes = kmerWithCount.getKey().value;
            int kmerCount = (int)(kmerWithCount.getValue());

            int[] kmer = GenomeSequenceEncodingUtil.decodeKmerBytesToIntArray(kmerBytes);

            String read = GenomeSequenceEncodingUtil.getStringFromBPBytes(this.reads.get(kmer[0]));
            String otherRead = null;

            String kmerPrefix = read.substring(kmer[1], kmer[1] + (k - 1));

            byte[] prefix;
            prefix = putKmerReadAndIndexPairForKmerString(kmerPrefix, kmer[0], kmer[1]); // note these kmers are implicitly *k-1 characters long*. for De Bruijn graphs each node represents k-1 characters, with the EDGE from one k-1 node to another k-1 node representing the actual length k kmer.

            String kmerSuffix = read.substring(kmer[1] + 1, kmer[1] + (k - 1) + 1);

            byte[] suffix;
            suffix = putKmerReadAndIndexPairForKmerString(kmerSuffix, kmer[0], kmer[1] + 1); // note these kmers are implicitly *k-1 characters long*. for De Bruijn graphs each node represents k-1 characters, with the EDGE from one k-1 node to another k-1 node representing the actual length k kmer.

            Map<ValueComparedByteArray, Byte> suffixesForPrefix = debruijnMap.get(ValueComparedByteArray.from(prefix));
            if (suffixesForPrefix == null) {
                suffixesForPrefix = new HashMap<ValueComparedByteArray, Byte>();
                debruijnMap.put(ValueComparedByteArray.from(prefix), suffixesForPrefix);
            }

            //Integer kmerCount = kmerCoverageMap.get(kmerBytes);
            int count = ((int)suffixesForPrefix.getOrDefault(ValueComparedByteArray.from(suffix), Byte.valueOf((byte)0))) + kmerCount;
            suffixesForPrefix.put(ValueComparedByteArray.from(suffix), Byte.valueOf((byte)count));

            if (numKmersChecked % 10000 == 0 && VERBOSE) {
                System.out.println("Processed kmer " + (numKmersChecked+1) + " out of " + numKmersToCheck);
                //if (numKmersChecked % 250000 == 0) {
                    //System.gc();
                //}
            }
            numKmersChecked++;
        }

        // necessary: free up the memory now otherwise VM will hang on to this stuff
        kmerHashCodeToReadIndexAndSubstringIndexMap = null;
        //System.gc();
        kmerHashCodeToReadIndexAndSubstringIndexMap = new HashMap<Integer, List<ValueComparedByteArray>>();
    }

    /**
       After setupFromReadStringsAndSpecifiedK and processKmers are complete, call this method to get the de Bruijn map to pass to the constructor of DeBruijnGraph .

       @see org.wbdbga.DeBruijnKmers#setupFromReadStringsAndSpecifiedK
       @see org.wbdbga.DeBruijnKmers#processKmers
       @see org.wbdbga.DeBruijnGraph
       @return de Bruijn map to pass to the constructor of the DeBruijnGraph, mapping (k-1)mers to their (k-1)mer suffixes with coverage counts
       @throws NotInitializedException if debruijnMap is not yet available due to setupFromReadStringsAndSpecifiedK and processKmers not both being called before this method
     */
    public Map<ValueComparedByteArray, Map<ValueComparedByteArray, Byte>> getDeBruijnMap() {
        if (debruijnMap == null || debruijnMap.keySet().size() == 0) {
            throw new NotInitializedException("debruijnMap not yet created. call setupFromReadStringsAndSpecifiedK and processKmers before calling getDebruijnMap()");
        }

        return debruijnMap;
    }

    class NotInitializedException extends RuntimeException {
        public NotInitializedException(String message) {
            super(message);
        }
    }

    class CorruptedReadIndexAndSubstringIndexMapException extends RuntimeException {
        public CorruptedReadIndexAndSubstringIndexMapException(String message) {
            super(message);
        }
    }

    class OutOfRangeReadIndexException extends RuntimeException {
        public OutOfRangeReadIndexException(String message) {
            super(message);
        }
    }

    class OutOfRangeReadOffsetException extends RuntimeException {
        public OutOfRangeReadOffsetException(String message) {
            super(message);
        }
    }
}
