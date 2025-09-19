This folder contains reference genomes to generate fake error-prone reads (using GenerateNoisyReadsFromWholeGenomeFasta) from and practice using the assembler in the parent folder on. The assemblies can then be graded against the original source genomes using [QUAST](https://quast.sourceforge.net/) .

## GenerateNoisyReadsFromWholeGenomeFasta utility for generating fake data for assembly

Usage:
```
javac GenerateNoisyReadsFromWholeGenomeFasta.java
java GenerateNoisyReadsFromWholeGenomeFasta --input_fasta carsonella_ruddii_complete_genome_AP009180.1.fasta > ../examples/noisyCarsonellaRuddiiReads.txt
```
The noisy reads text file can then be used as input to DeBruijnGenomeAssembler in the parent folder to check ability to reassemble the genome from the noisy reads (fragments).

```
Supported arguments:
  --input_fasta: specify an input FASTA file containing a whole genome to generate synthetic reads from
 (e.g. carsonella_ruddii_complete_genome_AP009180.1.fasta from https://www.ncbi.nlm.nih.gov/nuccore/AP009180.1?report=fasta&log$=seqview&format=text )
  --error_probability: probability of single base pair error in synthetic reads generated (double; 0 <= error_probability <= 1; defaults to 0.01)
    e.g. "--error_probability 0.01" for 1% errors
  --read_length: number of base pairs long each synthetic read should be (integer 25 <= read length <= 200; defaults to 100 bp)
  --read_coverage: how many synthetic reads to generate based on target coverage == (number of reads)*(read length)/(genome length), average number of times we expect each base to be sampled (integer: 1 <= coverage; defaults to 30x coverage)
  --random_seed: some long, using the same values for the same input fasta gives a reproducible random sampling
```

To obtain `carsonella_ruddii_complete_genome_AP009180.1.fasta` and other example FASTA files from NIH NCBI, see "Example genomes" section below. I am tempted to include them as DNA cannot be copyrighted in the US, but it might be able to be in the UK or other countries and in any case these genomes are available for download from the US government at high speed for free so why not link to the copy over there.

## Example genomes from NIH NCBI

### Escherichia phage phiX174, complete genome (tiny)

**5,386 bp circular DNA, tiny; can be re-assembled in seconds.**

Command to download FASTA file which can be used with GenerateNoisyReadsFromWholeGenomeFasta to create synthetic reads for the DeBruijnGenomeAssembler:
```
wget 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=9626372&extrafeat=null&conwithfeat=on&hide-cdd=on' -O escherichia_phage_phiX174_complete_genome.fasta
```

To generate noisy synthetic 100bp reads with 1% errors 
(like a sequencer might output) to practice reassembling,
you can then run
```
javac GenerateNoisyReadsFromWholeGenomeFasta.java
java GenerateNoisyReadsFromWholeGenomeFasta --input_fasta escherichia_phage_phiX174_complete_genome.fasta > ../examples/noisyPhiX174Reads.txt
```

Then from the parent folder
```
./gradlew build
java -jar ./app/build/libs/app.jar < examples/noisyPhiX174Reads.txt
```

The QUAST ready output will be in `output.fasta` by default. Then, follow the [instructions regarding QUAST in the examples README](../examples/README.md#run-quast-to-check-our-assembly-quality-against-the-reference) to see how to check the quality of the assembly.

URL to read more in the browser or if the above FASTA download does not work: https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1

### AP009180.1 Carsonella ruddii complete genome

**159,662 bp circular DNA, small; can be re-assembled in seconds for specified k value, to minutes if searching multiple k values.**

Command to download FASTA file which can be used with GenerateNoisyReadsFromWholeGenomeFasta to create synthetic reads for the DeBruijnGenomeAssembler:
```
wget 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=116235183&extrafeat=null&conwithfeat=on&hide-cdd=on' -O carsonella_ruddii_complete_genome_AP009180.1.fasta
```

To generate noisy synthetic 100bp reads with 1% errors 
(like a sequencer might output) to practice reassembling,
you can then run
```
javac GenerateNoisyReadsFromWholeGenomeFasta.java
java GenerateNoisyReadsFromWholeGenomeFasta --input_fasta carsonella_ruddii_complete_genome_AP009180.1.fasta > ../examples/noisyCarsonellaRuddiiReads.txt
```

Then from the parent folder
```
./gradlew build
java -jar ./app/build/libs/app.jar < examples/noisyCarsonellaRuddiiReads.txt
```

The QUAST ready output will be in `output.fasta` by default. Then, follow the [instructions regarding QUAST in the examples README](../examples/README.md#run-quast-to-check-our-assembly-quality-against-the-reference) to see how to check the quality of the assembly.

(Note the challenge problem on p164 of Compeau, P. E., Pevzner, P. A. (2018). Bioinformatics algorithms : an active learning approach, 3rd ed. Active Learning Publishers. https://www.bioinformaticsalgorithms.org/ suggests that even with error free reads, using unpaired reads for de Bruijn graph assembly of Carsonella ruddii (only 159662bp) will fail and paired will be required but we can get a perfect assembly in the 1st contig above using unpaired noisy reads.)

URL to read more in the browser or if the above FASTA download does not work: https://www.ncbi.nlm.nih.gov/nuccore/AP009180.1

### NZ_CABFNN010000001.1 Escherichia coli J53 isolate whole genome sequence

**4,705,562 bp linear DNA, large; as fast as 10 minutes to reassemble on a fast laptop with single k value, but hours on default setting searching many k values.**

Command to download FASTA file which can be used with GenerateNoisyReadsFromWholeGenomeFasta to create synthetic reads for the DeBruijnGenomeAssembler:
```
wget 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1692789651&extrafeat=null&conwithfeat=on&hide-cdd=on' -O escherichia_coli_J53_isolate_whole_genome_sequence.fasta
```

To generate noisy synthetic 100bp reads with 1% errors 
(like a sequencer might output) to practice reassembling,
you can then run
```
javac GenerateNoisyReadsFromWholeGenomeFasta.java
java GenerateNoisyReadsFromWholeGenomeFasta --input_fasta escherichia_coli_J53_isolate_whole_genome_sequence.fasta --read_coverage 40 > ../examples/noisyEColiReads.txt
```

Then from the parent folder
```
./gradlew build
java -jar ./app/build/libs/app.jar < examples/noisyEColiReads.txt
```
For these large megabase genomes, it will take a long time to let the program try different values of k, you may skip that by running above with `--k 45` (which is expected to be a decent value):
```
./gradlew build
java -jar ./app/build/libs/app.jar --k 45 < examples/noisyEColiReads.txt
```

You may also want to increase the available memory (I increase to 12GB below assuming 16GB available, adjust as appropriate) and make the most of it by adjusting the GC settings and to wrap in a `time` command to see how long it takes. To do that, add `time java -verbosegc -XX:-UseGCOverheadLimit -Xmx12288M -Xms12288M -XX:+UseSerialGC` so that the command becomes:
```
time java -verbosegc -XX:-UseGCOverheadLimit -Xmx12288M -Xms12288M -XX:+UseSerialGC -jar ./app/build/libs/app.jar --k 45 < examples/noisyEColiReads.txt
```

The QUAST ready output will be in `output.fasta` by default. Then, follow the [instructions regarding QUAST in the examples README](../examples/README.md#run-quast-to-check-our-assembly-quality-against-the-reference) to see how to check the quality of the assembly.

URL to read more in the browser or if the above FASTA download does not work: https://www.ncbi.nlm.nih.gov/nuccore/NZ_CABFNN010000001.1

### CP023474.1 Brevibacillus brevis X23 chromosome, complete genome (stretch goal, the Bioinformatics Algorithm textbook says assembly with unpaired reads should struggle with C. ruddii at 160K genome size and recommends comparing with paired read assembly; but we do C. ruddii at 160Kbp perfectly and E. coli at much larger 4.7Mbp fairly well with unpaired, so we try this; but we do finally hit barriers here at 6,643,437 bp with short unpaired reads!)

**6,643,437 bp circular DNA; large, hours to reassemble.**

B. brevis produces many medicinally and agriculturally important antibiotics, e.g. [tyrocidine](https://en.wikipedia.org/wiki/Tyrocidine), [gramicidin](https://en.wikipedia.org/wiki/Gramicidin), gratisin, and edeine. Note: these peptide antibiotics are not directly encoded in the DNA (to complicate things they are not translated from mRNA transcripts by ribosomes), they are indirectly encoded via nonribosomal peptide (NRP) synthetases; only the code describing the NRP synthetases is in the DNA (not their products).

See the associated paper for this strain/genome at Chen, W., Wang, Y., Li, D., Li, L., Xiao, Q., & Zhou, Q. (2012). Draft genome sequence of Brevibacillus brevis strain X23, a biocontrol agent against bacterial wilt. https://pmc.ncbi.nlm.nih.gov/articles/PMC3497470/ .

Command to download FASTA file which can be used with GenerateNoisyReadsFromWholeGenomeFasta to create synthetic reads for the DeBruijnGenomeAssembler:
```
wget 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1247079690&extrafeat=null&conwithfeat=on&hide-cdd=on' -O CP023474.1_Brevibacillus_brevis_X23_complete_genome.fasta
```

To generate noisy synthetic 100bp reads with 1% errors 
(like a sequencer might output) to practice reassembling,
you can then run
```
javac GenerateNoisyReadsFromWholeGenomeFasta.java
java GenerateNoisyReadsFromWholeGenomeFasta --input_fasta CP023474.1_Brevibacillus_brevis_X23_complete_genome.fasta --read_coverage 40 > ../examples/noisyBBrevisReads.txt
```

Then from the parent folder, we run the standard jar, but with additional memory and serial GC to make the most of the memory we have:
(>= 16GB memory recommended if doing ths locally, below we show allocating ~12GB out of 16GB; and we run `time` to keep track of how long it takes, this is much larger than E. coli so we will skip the slow kmer size grid search and specify the kmer fragment size we found to be good on E. coli, k==45)
```
./gradlew build
time java -verbosegc -XX:-UseGCOverheadLimit -Xmx12288M -Xms12288M -XX:+UseSerialGC -jar ./app/build/libs/app.jar --k 45 < examples/noisyBBrevisReads.txt
```
Note it saves potentially hours at small quality cost to add `--k 45` to the above to skip searching for the best k value and use an established decent value from E. coli and other megabase genomes, but that can be omitted if you can wait many hours for such a search (e.g. if you are doing this on a dedicated Google cloud machine.)

The QUAST ready output will be in `output.fasta` by default. Then, follow the [instructions regardng QUAST in the examples README](../examples/README.md#run-quast-to-check-our-assembly-quality-against-the-reference) to see how to check the quality of the assembly.

URL to read more in the browser or if the above FASTA download does not work: https://www.ncbi.nlm.nih.gov/nuccore/CP023474.1

### License

Copyright (C) 2025 William Bruns

This documentation comes with ABSOLUTELY NO WARRANTY and is provided "AS IS".
This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See [the full license](../LICENSE) .