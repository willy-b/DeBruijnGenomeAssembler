### Generate 1% error 100bp 30x coverage noisy reads of carsonella_ruddii_complete_genome_AP009180.1.fasta

`GenerateNoisyReadsFromWholeGenomeFasta` can generate noisy reads from any genome FASTA (run with `--help` for details) 
but we can start by generating noisy reads of the Carsonella Ruddii genome by default, using https://www.ncbi.nlm.nih.gov/nuccore/AP009180.1?report=fasta&log$=seqview&format=text as a source.

(See the [GenerateNoisyReadsFromWholeGenomeFasta README](../util/README.md) for examples of how to download other genomes (the famous E. coli and the antibiotic producing B. brevis for example) .)

```
cd DeBruijnGenomeAssembler/util

# download the genome to generate the fake sequencer reads from 
# and that we can grade ourselves against at the end using QUAST
wget 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=116235183&extrafeat=null&conwithfeat=on&hide-cdd=on' -O carsonella_ruddii_complete_genome_AP009180.1.fasta

# generate the fake sequencer reads
javac GenerateNoisyReadsFromWholeGenomeFasta.java
java GenerateNoisyReadsFromWholeGenomeFasta --input_fasta carsonella_ruddii_complete_genome_AP009180.1.fasta --error_probability 0.01 --read_length 100 --read_coverage 30 --random_seed 1 > ../examples/noisyCarsonellaRuddiiReads.txt
```

### Assemble Carsonella ruddii reads into perfect single contig assembly (kmer size selector will search based on reference-free n50 estimate and choose k=45)

```
[liveuser@localhost-live ~]$ cd DeBruijnGenomeAssembler
[liveuser@localhost-live DeBruijnGenomeAssembler]$ ./gradlew build # or `gradle build` if a compatible gradle is already installed
[...]
[liveuser@localhost-live DeBruijnGenomeAssembler]$ java -jar ./app/build/libs/app.jar < examples/noisyCarsonellaRuddiiReads.txt
numReadsSkipped: 0
Trying k == 30
started loadKmerCoverageMapFromReads
kmer "AATTGCAACTAAACATTAATATTATTTTTT" shares a hashCode with 0 other kmers
kmer "ATTTTAAAAAAATATTGATAAGATCTTGCT" shares a hashCode with 1 other kmers
Generated kmers from reads for k == 30, processing...
Removed 500000 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 527529)
Removed 345183 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 182346)
Removed 0 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 182346)
estimatedGenomeSizeBp: 182346
estimatedAverageCoverageOfGenome: 18.650027968806555
started loadKmerCoverageMapFromReads
Starting to generate De Bruijn map [(k-1)mer to (k-1)mer edges] from kmers
Processed kmer 10000 out of 1027529
Processed kmer 20000 out of 1027529
Processed kmer 30000 out of 1027529
Processed kmer 40000 out of 1027529
Processed kmer 50000 out of 1027529
Processed kmer 60000 out of 1027529
Processed kmer 70000 out of 1027529
Processed kmer 80000 out of 1027529
Processed kmer 90000 out of 1027529
Processed kmer 100000 out of 1027529
Processed kmer 110000 out of 1027529
Processed kmer 120000 out of 1027529
Processed kmer 130000 out of 1027529
Processed kmer 140000 out of 1027529
Processed kmer 150000 out of 1027529
Processed kmer 160000 out of 1027529
Processed kmer 170000 out of 1027529
Processed kmer 180000 out of 1027529
Processed kmer 190000 out of 1027529
Processed kmer 200000 out of 1027529
Processed kmer 210000 out of 1027529
Processed kmer 220000 out of 1027529
Processed kmer 230000 out of 1027529
Processed kmer 240000 out of 1027529
Processed kmer 250000 out of 1027529
Processed kmer 260000 out of 1027529
Processed kmer 270000 out of 1027529
Processed kmer 280000 out of 1027529
Processed kmer 290000 out of 1027529
Processed kmer 300000 out of 1027529
Processed kmer 310000 out of 1027529
Processed kmer 320000 out of 1027529
Processed kmer 330000 out of 1027529
Processed kmer 340000 out of 1027529
Processed kmer 350000 out of 1027529
Processed kmer 360000 out of 1027529
Processed kmer 370000 out of 1027529
Processed kmer 380000 out of 1027529
Processed kmer 390000 out of 1027529
Processed kmer 400000 out of 1027529
Processed kmer 410000 out of 1027529
Processed kmer 420000 out of 1027529
Processed kmer 430000 out of 1027529
Processed kmer 440000 out of 1027529
Processed kmer 450000 out of 1027529
Processed kmer 460000 out of 1027529
Processed kmer 470000 out of 1027529
Processed kmer 480000 out of 1027529
Processed kmer 490000 out of 1027529
Processed kmer 500000 out of 1027529
Processed kmer 510000 out of 1027529
Processed kmer 520000 out of 1027529
Processed kmer 530000 out of 1027529
Processed kmer 540000 out of 1027529
Processed kmer 550000 out of 1027529
Processed kmer 560000 out of 1027529
Processed kmer 570000 out of 1027529
Processed kmer 580000 out of 1027529
Processed kmer 590000 out of 1027529
Processed kmer 600000 out of 1027529
Processed kmer 610000 out of 1027529
Processed kmer 620000 out of 1027529
Processed kmer 630000 out of 1027529
Processed kmer 640000 out of 1027529
Processed kmer 650000 out of 1027529
Processed kmer 660000 out of 1027529
Processed kmer 670000 out of 1027529
Processed kmer 680000 out of 1027529
Processed kmer 690000 out of 1027529
Processed kmer 700000 out of 1027529
Processed kmer 710000 out of 1027529
Processed kmer 720000 out of 1027529
Processed kmer 730000 out of 1027529
Processed kmer 740000 out of 1027529
Processed kmer 750000 out of 1027529
Processed kmer 760000 out of 1027529
Processed kmer 770000 out of 1027529
Processed kmer 780000 out of 1027529
Processed kmer 790000 out of 1027529
Processed kmer 800000 out of 1027529
Processed kmer 810000 out of 1027529
Processed kmer 820000 out of 1027529
Processed kmer 830000 out of 1027529
Processed kmer 840000 out of 1027529
Processed kmer 850000 out of 1027529
Processed kmer 860000 out of 1027529
Processed kmer 870000 out of 1027529
Processed kmer 880000 out of 1027529
Processed kmer 890000 out of 1027529
Processed kmer 900000 out of 1027529
Processed kmer 910000 out of 1027529
Processed kmer 920000 out of 1027529
Processed kmer 930000 out of 1027529
Processed kmer 940000 out of 1027529
Processed kmer 950000 out of 1027529
Processed kmer 960000 out of 1027529
Processed kmer 970000 out of 1027529
Processed kmer 980000 out of 1027529
Processed kmer 990000 out of 1027529
Processed kmer 1000000 out of 1027529
Processed kmer 1010000 out of 1027529
Processed kmer 1020000 out of 1027529
Setting up the graph...
processed 0 out of 1002324 nodes
processed 10000 out of 1002446 nodes
processed 20000 out of 1002569 nodes
processed 30000 out of 1002674 nodes
processed 40000 out of 1002778 nodes
processed 50000 out of 1002880 nodes
processed 60000 out of 1003002 nodes
processed 70000 out of 1003114 nodes
processed 80000 out of 1003206 nodes
processed 90000 out of 1003321 nodes
processed 100000 out of 1003433 nodes
processed 110000 out of 1003549 nodes
processed 120000 out of 1003676 nodes
processed 130000 out of 1003777 nodes
processed 140000 out of 1003883 nodes
processed 150000 out of 1004003 nodes
processed 160000 out of 1004122 nodes
processed 170000 out of 1004239 nodes
processed 180000 out of 1004350 nodes
processed 190000 out of 1004462 nodes
processed 200000 out of 1004572 nodes
processed 210000 out of 1004694 nodes
processed 220000 out of 1004803 nodes
processed 230000 out of 1004915 nodes
processed 240000 out of 1005018 nodes
processed 250000 out of 1005150 nodes
processed 260000 out of 1005261 nodes
processed 270000 out of 1005373 nodes
processed 280000 out of 1005470 nodes
processed 290000 out of 1005594 nodes
processed 300000 out of 1005709 nodes
processed 310000 out of 1005823 nodes
processed 320000 out of 1005942 nodes
processed 330000 out of 1006069 nodes
processed 340000 out of 1006187 nodes
processed 350000 out of 1006304 nodes
processed 360000 out of 1006402 nodes
processed 370000 out of 1006523 nodes
processed 380000 out of 1006654 nodes
processed 390000 out of 1006765 nodes
processed 400000 out of 1006865 nodes
processed 410000 out of 1006973 nodes
processed 420000 out of 1007108 nodes
processed 430000 out of 1007205 nodes
processed 440000 out of 1007324 nodes
processed 450000 out of 1007450 nodes
processed 460000 out of 1007549 nodes
processed 470000 out of 1007660 nodes
processed 480000 out of 1007781 nodes
processed 490000 out of 1007894 nodes
processed 500000 out of 1007991 nodes
processed 510000 out of 1008106 nodes
processed 520000 out of 1008209 nodes
processed 530000 out of 1008335 nodes
processed 540000 out of 1008439 nodes
processed 550000 out of 1008559 nodes
processed 560000 out of 1008680 nodes
processed 570000 out of 1008801 nodes
processed 580000 out of 1008904 nodes
processed 590000 out of 1009018 nodes
processed 600000 out of 1009140 nodes
processed 610000 out of 1009261 nodes
processed 620000 out of 1009379 nodes
processed 630000 out of 1009481 nodes
processed 640000 out of 1009604 nodes
processed 650000 out of 1009717 nodes
processed 660000 out of 1009846 nodes
processed 670000 out of 1009960 nodes
processed 680000 out of 1010056 nodes
processed 690000 out of 1010187 nodes
processed 700000 out of 1010296 nodes
processed 710000 out of 1010402 nodes
processed 720000 out of 1010516 nodes
processed 730000 out of 1010638 nodes
processed 740000 out of 1010742 nodes
processed 750000 out of 1010849 nodes
processed 760000 out of 1010970 nodes
processed 770000 out of 1011066 nodes
processed 780000 out of 1011169 nodes
processed 790000 out of 1011287 nodes
processed 800000 out of 1011397 nodes
processed 810000 out of 1011516 nodes
processed 820000 out of 1011645 nodes
processed 830000 out of 1011752 nodes
processed 840000 out of 1011866 nodes
processed 850000 out of 1011976 nodes
processed 860000 out of 1012080 nodes
processed 870000 out of 1012183 nodes
processed 880000 out of 1012296 nodes
processed 890000 out of 1012400 nodes
processed 900000 out of 1012529 nodes
processed 910000 out of 1012647 nodes
processed 920000 out of 1012755 nodes
processed 930000 out of 1012867 nodes
processed 940000 out of 1012976 nodes
processed 950000 out of 1013098 nodes
processed 960000 out of 1013229 nodes
processed 970000 out of 1013343 nodes
processed 980000 out of 1013449 nodes
processed 990000 out of 1013566 nodes
processed 1000000 out of 1013677 nodes
processed 1010000 out of 1013702 nodes
average read coverage: 5.796754154870568
Built DeBruijn graph (of k-1 length strings for nodes, kmer strings as edges) with 1013702 nodes and 1027529 edges
22930 tip nodes found
420462 tip nodes removed
758 tip nodes found
7619 tip nodes removed
Found a candidate. k == 30 has n50 == 60727 (NOT computed vs a reference genome but is relative to DeBruijn graph formed from reads as this is a de novo assembler).
Trying k == 35
started loadKmerCoverageMapFromReads
kmer "AATTGCAACTAAACATTAATATTATTTTTTTTCCT" shares a hashCode with 0 other kmers
kmer "TTATAAATGTAAATATTAAATACAGCAAAATAGAT" shares a hashCode with 1 other kmers
Generated kmers from reads for k == 35, processing...
Removed 500000 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 582269)
Removed 401833 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 180436)
Removed 0 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 180436)
estimatedGenomeSizeBp: 180436
estimatedAverageCoverageOfGenome: 17.520162273603937
started loadKmerCoverageMapFromReads
Starting to generate De Bruijn map [(k-1)mer to (k-1)mer edges] from kmers
Processed kmer 10000 out of 1082269
Processed kmer 20000 out of 1082269
Processed kmer 30000 out of 1082269
Processed kmer 40000 out of 1082269
Processed kmer 50000 out of 1082269
Processed kmer 60000 out of 1082269
Processed kmer 70000 out of 1082269
Processed kmer 80000 out of 1082269
Processed kmer 90000 out of 1082269
Processed kmer 100000 out of 1082269
Processed kmer 110000 out of 1082269
Processed kmer 120000 out of 1082269
Processed kmer 130000 out of 1082269
Processed kmer 140000 out of 1082269
Processed kmer 150000 out of 1082269
Processed kmer 160000 out of 1082269
Processed kmer 170000 out of 1082269
Processed kmer 180000 out of 1082269
Processed kmer 190000 out of 1082269
Processed kmer 200000 out of 1082269
Processed kmer 210000 out of 1082269
Processed kmer 220000 out of 1082269
Processed kmer 230000 out of 1082269
Processed kmer 240000 out of 1082269
Processed kmer 250000 out of 1082269
Processed kmer 260000 out of 1082269
Processed kmer 270000 out of 1082269
Processed kmer 280000 out of 1082269
Processed kmer 290000 out of 1082269
Processed kmer 300000 out of 1082269
Processed kmer 310000 out of 1082269
Processed kmer 320000 out of 1082269
Processed kmer 330000 out of 1082269
Processed kmer 340000 out of 1082269
Processed kmer 350000 out of 1082269
Processed kmer 360000 out of 1082269
Processed kmer 370000 out of 1082269
Processed kmer 380000 out of 1082269
Processed kmer 390000 out of 1082269
Processed kmer 400000 out of 1082269
Processed kmer 410000 out of 1082269
Processed kmer 420000 out of 1082269
Processed kmer 430000 out of 1082269
Processed kmer 440000 out of 1082269
Processed kmer 450000 out of 1082269
Processed kmer 460000 out of 1082269
Processed kmer 470000 out of 1082269
Processed kmer 480000 out of 1082269
Processed kmer 490000 out of 1082269
Processed kmer 500000 out of 1082269
Processed kmer 510000 out of 1082269
Processed kmer 520000 out of 1082269
Processed kmer 530000 out of 1082269
Processed kmer 540000 out of 1082269
Processed kmer 550000 out of 1082269
Processed kmer 560000 out of 1082269
Processed kmer 570000 out of 1082269
Processed kmer 580000 out of 1082269
Processed kmer 590000 out of 1082269
Processed kmer 600000 out of 1082269
Processed kmer 610000 out of 1082269
Processed kmer 620000 out of 1082269
Processed kmer 630000 out of 1082269
Processed kmer 640000 out of 1082269
Processed kmer 650000 out of 1082269
Processed kmer 660000 out of 1082269
Processed kmer 670000 out of 1082269
Processed kmer 680000 out of 1082269
Processed kmer 690000 out of 1082269
Processed kmer 700000 out of 1082269
Processed kmer 710000 out of 1082269
Processed kmer 720000 out of 1082269
Processed kmer 730000 out of 1082269
Processed kmer 740000 out of 1082269
Processed kmer 750000 out of 1082269
Processed kmer 760000 out of 1082269
Processed kmer 770000 out of 1082269
Processed kmer 780000 out of 1082269
Processed kmer 790000 out of 1082269
Processed kmer 800000 out of 1082269
Processed kmer 810000 out of 1082269
Processed kmer 820000 out of 1082269
Processed kmer 830000 out of 1082269
Processed kmer 840000 out of 1082269
Processed kmer 850000 out of 1082269
Processed kmer 860000 out of 1082269
Processed kmer 870000 out of 1082269
Processed kmer 880000 out of 1082269
Processed kmer 890000 out of 1082269
Processed kmer 900000 out of 1082269
Processed kmer 910000 out of 1082269
Processed kmer 920000 out of 1082269
Processed kmer 930000 out of 1082269
Processed kmer 940000 out of 1082269
Processed kmer 950000 out of 1082269
Processed kmer 960000 out of 1082269
Processed kmer 970000 out of 1082269
Processed kmer 980000 out of 1082269
Processed kmer 990000 out of 1082269
Processed kmer 1000000 out of 1082269
Processed kmer 1010000 out of 1082269
Processed kmer 1020000 out of 1082269
Processed kmer 1030000 out of 1082269
Processed kmer 1040000 out of 1082269
Processed kmer 1050000 out of 1082269
Processed kmer 1060000 out of 1082269
Processed kmer 1070000 out of 1082269
Processed kmer 1080000 out of 1082269
Setting up the graph...
processed 0 out of 1059951 nodes
processed 10000 out of 1060082 nodes
processed 20000 out of 1060210 nodes
processed 30000 out of 1060338 nodes
processed 40000 out of 1060451 nodes
processed 50000 out of 1060576 nodes
processed 60000 out of 1060699 nodes
processed 70000 out of 1060820 nodes
processed 80000 out of 1060936 nodes
processed 90000 out of 1061070 nodes
processed 100000 out of 1061181 nodes
processed 110000 out of 1061304 nodes
processed 120000 out of 1061436 nodes
processed 130000 out of 1061572 nodes
processed 140000 out of 1061682 nodes
processed 150000 out of 1061809 nodes
processed 160000 out of 1061932 nodes
processed 170000 out of 1062065 nodes
processed 180000 out of 1062188 nodes
processed 190000 out of 1062315 nodes
processed 200000 out of 1062441 nodes
processed 210000 out of 1062560 nodes
processed 220000 out of 1062694 nodes
processed 230000 out of 1062816 nodes
processed 240000 out of 1062933 nodes
processed 250000 out of 1063054 nodes
processed 260000 out of 1063182 nodes
processed 270000 out of 1063315 nodes
processed 280000 out of 1063445 nodes
processed 290000 out of 1063558 nodes
processed 300000 out of 1063678 nodes
processed 310000 out of 1063815 nodes
processed 320000 out of 1063939 nodes
processed 330000 out of 1064065 nodes
processed 340000 out of 1064193 nodes
processed 350000 out of 1064333 nodes
processed 360000 out of 1064463 nodes
processed 370000 out of 1064587 nodes
processed 380000 out of 1064705 nodes
processed 390000 out of 1064831 nodes
processed 400000 out of 1064969 nodes
processed 410000 out of 1065093 nodes
processed 420000 out of 1065211 nodes
processed 430000 out of 1065341 nodes
processed 440000 out of 1065474 nodes
processed 450000 out of 1065600 nodes
processed 460000 out of 1065714 nodes
processed 470000 out of 1065841 nodes
processed 480000 out of 1065970 nodes
processed 490000 out of 1066087 nodes
processed 500000 out of 1066198 nodes
processed 510000 out of 1066333 nodes
processed 520000 out of 1066445 nodes
processed 530000 out of 1066561 nodes
processed 540000 out of 1066679 nodes
processed 550000 out of 1066793 nodes
processed 560000 out of 1066924 nodes
processed 570000 out of 1067041 nodes
processed 580000 out of 1067162 nodes
processed 590000 out of 1067301 nodes
processed 600000 out of 1067432 nodes
processed 610000 out of 1067544 nodes
processed 620000 out of 1067670 nodes
processed 630000 out of 1067800 nodes
processed 640000 out of 1067923 nodes
processed 650000 out of 1068056 nodes
processed 660000 out of 1068178 nodes
processed 670000 out of 1068292 nodes
processed 680000 out of 1068429 nodes
processed 690000 out of 1068557 nodes
processed 700000 out of 1068698 nodes
processed 710000 out of 1068808 nodes
processed 720000 out of 1068918 nodes
processed 730000 out of 1069058 nodes
processed 740000 out of 1069171 nodes
processed 750000 out of 1069283 nodes
processed 760000 out of 1069405 nodes
processed 770000 out of 1069532 nodes
processed 780000 out of 1069654 nodes
processed 790000 out of 1069775 nodes
processed 800000 out of 1069900 nodes
processed 810000 out of 1070025 nodes
processed 820000 out of 1070135 nodes
processed 830000 out of 1070247 nodes
processed 840000 out of 1070387 nodes
processed 850000 out of 1070513 nodes
processed 860000 out of 1070651 nodes
processed 870000 out of 1070784 nodes
processed 880000 out of 1070905 nodes
processed 890000 out of 1071028 nodes
processed 900000 out of 1071153 nodes
processed 910000 out of 1071272 nodes
processed 920000 out of 1071394 nodes
processed 930000 out of 1071509 nodes
processed 940000 out of 1071628 nodes
processed 950000 out of 1071768 nodes
processed 960000 out of 1071902 nodes
processed 970000 out of 1072023 nodes
processed 980000 out of 1072143 nodes
processed 990000 out of 1072258 nodes
processed 1000000 out of 1072397 nodes
processed 1010000 out of 1072544 nodes
processed 1020000 out of 1072665 nodes
processed 1030000 out of 1072785 nodes
processed 1040000 out of 1072903 nodes
processed 1050000 out of 1073037 nodes
processed 1060000 out of 1073161 nodes
processed 1070000 out of 1073161 nodes
average read coverage: 5.0086466488460815
Built DeBruijn graph (of k-1 length strings for nodes, kmer strings as edges) with 1073161 nodes and 1082269 edges
26539 tip nodes found
573039 tip nodes removed
775 tip nodes found
6365 tip nodes removed
Found a candidate. k == 35 has n50 == 159696 (NOT computed vs a reference genome but is relative to DeBruijn graph formed from reads as this is a de novo assembler).
Trying k == 40
started loadKmerCoverageMapFromReads
kmer "AATTGCAACTAAACATTAATATTATTTTTTTTCCTACACA" shares a hashCode with 0 other kmers
kmer "TAAAATCATTACAAAAACTAAAATGTAAATTAAATTGATT" shares a hashCode with 1 other kmers
Generated kmers from reads for k == 40, processing...
Removed 500000 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 614606)
Removed 436493 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 178113)
Removed 0 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 178113)
estimatedGenomeSizeBp: 178113
estimatedAverageCoverageOfGenome: 16.404069326775698
started loadKmerCoverageMapFromReads
Starting to generate De Bruijn map [(k-1)mer to (k-1)mer edges] from kmers
Processed kmer 10000 out of 1114606
Processed kmer 20000 out of 1114606
Processed kmer 30000 out of 1114606
Processed kmer 40000 out of 1114606
Processed kmer 50000 out of 1114606
Processed kmer 60000 out of 1114606
Processed kmer 70000 out of 1114606
Processed kmer 80000 out of 1114606
Processed kmer 90000 out of 1114606
Processed kmer 100000 out of 1114606
Processed kmer 110000 out of 1114606
Processed kmer 120000 out of 1114606
Processed kmer 130000 out of 1114606
Processed kmer 140000 out of 1114606
Processed kmer 150000 out of 1114606
Processed kmer 160000 out of 1114606
Processed kmer 170000 out of 1114606
Processed kmer 180000 out of 1114606
Processed kmer 190000 out of 1114606
Processed kmer 200000 out of 1114606
Processed kmer 210000 out of 1114606
Processed kmer 220000 out of 1114606
Processed kmer 230000 out of 1114606
Processed kmer 240000 out of 1114606
Processed kmer 250000 out of 1114606
Processed kmer 260000 out of 1114606
Processed kmer 270000 out of 1114606
Processed kmer 280000 out of 1114606
Processed kmer 290000 out of 1114606
Processed kmer 300000 out of 1114606
Processed kmer 310000 out of 1114606
Processed kmer 320000 out of 1114606
Processed kmer 330000 out of 1114606
Processed kmer 340000 out of 1114606
Processed kmer 350000 out of 1114606
Processed kmer 360000 out of 1114606
Processed kmer 370000 out of 1114606
Processed kmer 380000 out of 1114606
Processed kmer 390000 out of 1114606
Processed kmer 400000 out of 1114606
Processed kmer 410000 out of 1114606
Processed kmer 420000 out of 1114606
Processed kmer 430000 out of 1114606
Processed kmer 440000 out of 1114606
Processed kmer 450000 out of 1114606
Processed kmer 460000 out of 1114606
Processed kmer 470000 out of 1114606
Processed kmer 480000 out of 1114606
Processed kmer 490000 out of 1114606
Processed kmer 500000 out of 1114606
Processed kmer 510000 out of 1114606
Processed kmer 520000 out of 1114606
Processed kmer 530000 out of 1114606
Processed kmer 540000 out of 1114606
Processed kmer 550000 out of 1114606
Processed kmer 560000 out of 1114606
Processed kmer 570000 out of 1114606
Processed kmer 580000 out of 1114606
Processed kmer 590000 out of 1114606
Processed kmer 600000 out of 1114606
Processed kmer 610000 out of 1114606
Processed kmer 620000 out of 1114606
Processed kmer 630000 out of 1114606
Processed kmer 640000 out of 1114606
Processed kmer 650000 out of 1114606
Processed kmer 660000 out of 1114606
Processed kmer 670000 out of 1114606
Processed kmer 680000 out of 1114606
Processed kmer 690000 out of 1114606
Processed kmer 700000 out of 1114606
Processed kmer 710000 out of 1114606
Processed kmer 720000 out of 1114606
Processed kmer 730000 out of 1114606
Processed kmer 740000 out of 1114606
Processed kmer 750000 out of 1114606
Processed kmer 760000 out of 1114606
Processed kmer 770000 out of 1114606
Processed kmer 780000 out of 1114606
Processed kmer 790000 out of 1114606
Processed kmer 800000 out of 1114606
Processed kmer 810000 out of 1114606
Processed kmer 820000 out of 1114606
Processed kmer 830000 out of 1114606
Processed kmer 840000 out of 1114606
Processed kmer 850000 out of 1114606
Processed kmer 860000 out of 1114606
Processed kmer 870000 out of 1114606
Processed kmer 880000 out of 1114606
Processed kmer 890000 out of 1114606
Processed kmer 900000 out of 1114606
Processed kmer 910000 out of 1114606
Processed kmer 920000 out of 1114606
Processed kmer 930000 out of 1114606
Processed kmer 940000 out of 1114606
Processed kmer 950000 out of 1114606
Processed kmer 960000 out of 1114606
Processed kmer 970000 out of 1114606
Processed kmer 980000 out of 1114606
Processed kmer 990000 out of 1114606
Processed kmer 1000000 out of 1114606
Processed kmer 1010000 out of 1114606
Processed kmer 1020000 out of 1114606
Processed kmer 1030000 out of 1114606
Processed kmer 1040000 out of 1114606
Processed kmer 1050000 out of 1114606
Processed kmer 1060000 out of 1114606
Processed kmer 1070000 out of 1114606
Processed kmer 1080000 out of 1114606
Processed kmer 1090000 out of 1114606
Processed kmer 1100000 out of 1114606
Processed kmer 1110000 out of 1114606
Setting up the graph...
processed 0 out of 1094962 nodes
processed 10000 out of 1095102 nodes
processed 20000 out of 1095237 nodes
processed 30000 out of 1095382 nodes
processed 40000 out of 1095508 nodes
processed 50000 out of 1095643 nodes
processed 60000 out of 1095781 nodes
processed 70000 out of 1095908 nodes
processed 80000 out of 1096048 nodes
processed 90000 out of 1096189 nodes
processed 100000 out of 1096322 nodes
processed 110000 out of 1096449 nodes
processed 120000 out of 1096584 nodes
processed 130000 out of 1096725 nodes
processed 140000 out of 1096845 nodes
processed 150000 out of 1096977 nodes
processed 160000 out of 1097114 nodes
processed 170000 out of 1097247 nodes
processed 180000 out of 1097391 nodes
processed 190000 out of 1097531 nodes
processed 200000 out of 1097656 nodes
processed 210000 out of 1097801 nodes
processed 220000 out of 1097936 nodes
processed 230000 out of 1098079 nodes
processed 240000 out of 1098206 nodes
processed 250000 out of 1098332 nodes
processed 260000 out of 1098456 nodes
processed 270000 out of 1098602 nodes
processed 280000 out of 1098739 nodes
processed 290000 out of 1098881 nodes
processed 300000 out of 1099009 nodes
processed 310000 out of 1099144 nodes
processed 320000 out of 1099294 nodes
processed 330000 out of 1099431 nodes
processed 340000 out of 1099564 nodes
processed 350000 out of 1099713 nodes
processed 360000 out of 1099861 nodes
processed 370000 out of 1100003 nodes
processed 380000 out of 1100142 nodes
processed 390000 out of 1100284 nodes
processed 400000 out of 1100406 nodes
processed 410000 out of 1100548 nodes
processed 420000 out of 1100697 nodes
processed 430000 out of 1100836 nodes
processed 440000 out of 1100963 nodes
processed 450000 out of 1101102 nodes
processed 460000 out of 1101243 nodes
processed 470000 out of 1101361 nodes
processed 480000 out of 1101503 nodes
processed 490000 out of 1101655 nodes
processed 500000 out of 1101782 nodes
processed 510000 out of 1101904 nodes
processed 520000 out of 1102035 nodes
processed 530000 out of 1102189 nodes
processed 540000 out of 1102320 nodes
processed 550000 out of 1102442 nodes
processed 560000 out of 1102569 nodes
processed 570000 out of 1102705 nodes
processed 580000 out of 1102838 nodes
processed 590000 out of 1102964 nodes
processed 600000 out of 1103100 nodes
processed 610000 out of 1103246 nodes
processed 620000 out of 1103391 nodes
processed 630000 out of 1103520 nodes
processed 640000 out of 1103655 nodes
processed 650000 out of 1103799 nodes
processed 660000 out of 1103935 nodes
processed 670000 out of 1104071 nodes
processed 680000 out of 1104199 nodes
processed 690000 out of 1104326 nodes
processed 700000 out of 1104467 nodes
processed 710000 out of 1104609 nodes
processed 720000 out of 1104753 nodes
processed 730000 out of 1104883 nodes
processed 740000 out of 1105021 nodes
processed 750000 out of 1105161 nodes
processed 760000 out of 1105294 nodes
processed 770000 out of 1105418 nodes
processed 780000 out of 1105550 nodes
processed 790000 out of 1105686 nodes
processed 800000 out of 1105815 nodes
processed 810000 out of 1105962 nodes
processed 820000 out of 1106095 nodes
processed 830000 out of 1106234 nodes
processed 840000 out of 1106366 nodes
processed 850000 out of 1106488 nodes
processed 860000 out of 1106627 nodes
processed 870000 out of 1106772 nodes
processed 880000 out of 1106911 nodes
processed 890000 out of 1107063 nodes
processed 900000 out of 1107205 nodes
processed 910000 out of 1107337 nodes
processed 920000 out of 1107471 nodes
processed 930000 out of 1107612 nodes
processed 940000 out of 1107750 nodes
processed 950000 out of 1107883 nodes
processed 960000 out of 1108020 nodes
processed 970000 out of 1108140 nodes
processed 980000 out of 1108295 nodes
processed 990000 out of 1108450 nodes
processed 1000000 out of 1108580 nodes
processed 1010000 out of 1108718 nodes
processed 1020000 out of 1108840 nodes
processed 1030000 out of 1108980 nodes
processed 1040000 out of 1109129 nodes
processed 1050000 out of 1109280 nodes
processed 1060000 out of 1109419 nodes
processed 1070000 out of 1109552 nodes
processed 1080000 out of 1109687 nodes
processed 1090000 out of 1109830 nodes
processed 1100000 out of 1109894 nodes
average read coverage: 4.402509048040294
Built DeBruijn graph (of k-1 length strings for nodes, kmer strings as edges) with 1109894 nodes and 1114606 edges
29948 tip nodes found
712586 tip nodes removed
2616 tip nodes found
12335 tip nodes removed
Found a candidate. k == 40 has n50 == 159701 (NOT computed vs a reference genome but is relative to DeBruijn graph formed from reads as this is a de novo assembler).
Trying k == 45
started loadKmerCoverageMapFromReads
kmer "AATTGCAACTAAACATTAATATTATTTTTTTTCCTACACATTCAA" shares a hashCode with 0 other kmers
kmer "GTTTTACTATCAAATTTAATATTATTTTTAATACAAAAATTTAAT" shares a hashCode with 1 other kmers
Generated kmers from reads for k == 45, processing...
Removed 500000 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 625465)
Removed 449891 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 175574)
Removed 0 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 175574)
estimatedGenomeSizeBp: 175574
estimatedAverageCoverageOfGenome: 15.27725061797305
started loadKmerCoverageMapFromReads
Starting to generate De Bruijn map [(k-1)mer to (k-1)mer edges] from kmers
Processed kmer 10000 out of 1125465
Processed kmer 20000 out of 1125465
Processed kmer 30000 out of 1125465
Processed kmer 40000 out of 1125465
Processed kmer 50000 out of 1125465
Processed kmer 60000 out of 1125465
Processed kmer 70000 out of 1125465
Processed kmer 80000 out of 1125465
Processed kmer 90000 out of 1125465
Processed kmer 100000 out of 1125465
Processed kmer 110000 out of 1125465
Processed kmer 120000 out of 1125465
Processed kmer 130000 out of 1125465
Processed kmer 140000 out of 1125465
Processed kmer 150000 out of 1125465
Processed kmer 160000 out of 1125465
Processed kmer 170000 out of 1125465
Processed kmer 180000 out of 1125465
Processed kmer 190000 out of 1125465
Processed kmer 200000 out of 1125465
Processed kmer 210000 out of 1125465
Processed kmer 220000 out of 1125465
Processed kmer 230000 out of 1125465
Processed kmer 240000 out of 1125465
Processed kmer 250000 out of 1125465
Processed kmer 260000 out of 1125465
Processed kmer 270000 out of 1125465
Processed kmer 280000 out of 1125465
Processed kmer 290000 out of 1125465
Processed kmer 300000 out of 1125465
Processed kmer 310000 out of 1125465
Processed kmer 320000 out of 1125465
Processed kmer 330000 out of 1125465
Processed kmer 340000 out of 1125465
Processed kmer 350000 out of 1125465
Processed kmer 360000 out of 1125465
Processed kmer 370000 out of 1125465
Processed kmer 380000 out of 1125465
Processed kmer 390000 out of 1125465
Processed kmer 400000 out of 1125465
Processed kmer 410000 out of 1125465
Processed kmer 420000 out of 1125465
Processed kmer 430000 out of 1125465
Processed kmer 440000 out of 1125465
Processed kmer 450000 out of 1125465
Processed kmer 460000 out of 1125465
Processed kmer 470000 out of 1125465
Processed kmer 480000 out of 1125465
Processed kmer 490000 out of 1125465
Processed kmer 500000 out of 1125465
Processed kmer 510000 out of 1125465
Processed kmer 520000 out of 1125465
Processed kmer 530000 out of 1125465
Processed kmer 540000 out of 1125465
Processed kmer 550000 out of 1125465
Processed kmer 560000 out of 1125465
Processed kmer 570000 out of 1125465
Processed kmer 580000 out of 1125465
Processed kmer 590000 out of 1125465
Processed kmer 600000 out of 1125465
Processed kmer 610000 out of 1125465
Processed kmer 620000 out of 1125465
Processed kmer 630000 out of 1125465
Processed kmer 640000 out of 1125465
Processed kmer 650000 out of 1125465
Processed kmer 660000 out of 1125465
Processed kmer 670000 out of 1125465
Processed kmer 680000 out of 1125465
Processed kmer 690000 out of 1125465
Processed kmer 700000 out of 1125465
Processed kmer 710000 out of 1125465
Processed kmer 720000 out of 1125465
Processed kmer 730000 out of 1125465
Processed kmer 740000 out of 1125465
Processed kmer 750000 out of 1125465
Processed kmer 760000 out of 1125465
Processed kmer 770000 out of 1125465
Processed kmer 780000 out of 1125465
Processed kmer 790000 out of 1125465
Processed kmer 800000 out of 1125465
Processed kmer 810000 out of 1125465
Processed kmer 820000 out of 1125465
Processed kmer 830000 out of 1125465
Processed kmer 840000 out of 1125465
Processed kmer 850000 out of 1125465
Processed kmer 860000 out of 1125465
Processed kmer 870000 out of 1125465
Processed kmer 880000 out of 1125465
Processed kmer 890000 out of 1125465
Processed kmer 900000 out of 1125465
Processed kmer 910000 out of 1125465
Processed kmer 920000 out of 1125465
Processed kmer 930000 out of 1125465
Processed kmer 940000 out of 1125465
Processed kmer 950000 out of 1125465
Processed kmer 960000 out of 1125465
Processed kmer 970000 out of 1125465
Processed kmer 980000 out of 1125465
Processed kmer 990000 out of 1125465
Processed kmer 1000000 out of 1125465
Processed kmer 1010000 out of 1125465
Processed kmer 1020000 out of 1125465
Processed kmer 1030000 out of 1125465
Processed kmer 1040000 out of 1125465
Processed kmer 1050000 out of 1125465
Processed kmer 1060000 out of 1125465
Processed kmer 1070000 out of 1125465
Processed kmer 1080000 out of 1125465
Processed kmer 1090000 out of 1125465
Processed kmer 1100000 out of 1125465
Processed kmer 1110000 out of 1125465
Processed kmer 1120000 out of 1125465
Setting up the graph...
processed 0 out of 1108366 nodes
processed 10000 out of 1108520 nodes
processed 20000 out of 1108673 nodes
processed 30000 out of 1108831 nodes
processed 40000 out of 1108971 nodes
processed 50000 out of 1109119 nodes
processed 60000 out of 1109266 nodes
processed 70000 out of 1109418 nodes
processed 80000 out of 1109572 nodes
processed 90000 out of 1109714 nodes
processed 100000 out of 1109860 nodes
processed 110000 out of 1110007 nodes
processed 120000 out of 1110149 nodes
processed 130000 out of 1110310 nodes
processed 140000 out of 1110438 nodes
processed 150000 out of 1110589 nodes
processed 160000 out of 1110741 nodes
processed 170000 out of 1110886 nodes
processed 180000 out of 1111044 nodes
processed 190000 out of 1111196 nodes
processed 200000 out of 1111333 nodes
processed 210000 out of 1111489 nodes
processed 220000 out of 1111629 nodes
processed 230000 out of 1111782 nodes
processed 240000 out of 1111928 nodes
processed 250000 out of 1112070 nodes
processed 260000 out of 1112219 nodes
processed 270000 out of 1112360 nodes
processed 280000 out of 1112519 nodes
processed 290000 out of 1112676 nodes
processed 300000 out of 1112817 nodes
processed 310000 out of 1112955 nodes
processed 320000 out of 1113116 nodes
processed 330000 out of 1113266 nodes
processed 340000 out of 1113426 nodes
processed 350000 out of 1113568 nodes
processed 360000 out of 1113726 nodes
processed 370000 out of 1113886 nodes
processed 380000 out of 1114047 nodes
processed 390000 out of 1114197 nodes
processed 400000 out of 1114341 nodes
processed 410000 out of 1114490 nodes
processed 420000 out of 1114645 nodes
processed 430000 out of 1114801 nodes
processed 440000 out of 1114949 nodes
processed 450000 out of 1115090 nodes
processed 460000 out of 1115250 nodes
processed 470000 out of 1115406 nodes
processed 480000 out of 1115535 nodes
processed 490000 out of 1115690 nodes
processed 500000 out of 1115847 nodes
processed 510000 out of 1115993 nodes
processed 520000 out of 1116126 nodes
processed 530000 out of 1116282 nodes
processed 540000 out of 1116432 nodes
processed 550000 out of 1116583 nodes
processed 560000 out of 1116732 nodes
processed 570000 out of 1116862 nodes
processed 580000 out of 1117008 nodes
processed 590000 out of 1117153 nodes
processed 600000 out of 1117287 nodes
processed 610000 out of 1117442 nodes
processed 620000 out of 1117586 nodes
processed 630000 out of 1117745 nodes
processed 640000 out of 1117901 nodes
processed 650000 out of 1118046 nodes
processed 660000 out of 1118190 nodes
processed 670000 out of 1118356 nodes
processed 680000 out of 1118506 nodes
processed 690000 out of 1118636 nodes
processed 700000 out of 1118775 nodes
processed 710000 out of 1118936 nodes
processed 720000 out of 1119091 nodes
processed 730000 out of 1119250 nodes
processed 740000 out of 1119405 nodes
processed 750000 out of 1119553 nodes
processed 760000 out of 1119709 nodes
processed 770000 out of 1119844 nodes
processed 780000 out of 1119976 nodes
processed 790000 out of 1120137 nodes
processed 800000 out of 1120277 nodes
processed 810000 out of 1120425 nodes
processed 820000 out of 1120578 nodes
processed 830000 out of 1120723 nodes
processed 840000 out of 1120875 nodes
processed 850000 out of 1121020 nodes
processed 860000 out of 1121159 nodes
processed 870000 out of 1121308 nodes
processed 880000 out of 1121471 nodes
processed 890000 out of 1121616 nodes
processed 900000 out of 1121779 nodes
processed 910000 out of 1121935 nodes
processed 920000 out of 1122084 nodes
processed 930000 out of 1122233 nodes
processed 940000 out of 1122377 nodes
processed 950000 out of 1122530 nodes
processed 960000 out of 1122676 nodes
processed 970000 out of 1122827 nodes
processed 980000 out of 1122963 nodes
processed 990000 out of 1123110 nodes
processed 1000000 out of 1123286 nodes
processed 1010000 out of 1123428 nodes
processed 1020000 out of 1123579 nodes
processed 1030000 out of 1123723 nodes
processed 1040000 out of 1123858 nodes
processed 1050000 out of 1124023 nodes
processed 1060000 out of 1124186 nodes
processed 1070000 out of 1124336 nodes
processed 1080000 out of 1124477 nodes
processed 1090000 out of 1124628 nodes
processed 1100000 out of 1124776 nodes
processed 1110000 out of 1124911 nodes
processed 1120000 out of 1124911 nodes
average read coverage: 3.922543126618775
Built DeBruijn graph (of k-1 length strings for nodes, kmer strings as edges) with 1124911 nodes and 1125465 edges
33225 tip nodes found
830317 tip nodes removed
2219 tip nodes found
8299 tip nodes removed
Found a candidate. k == 45 has n50 == 159706 (NOT computed vs a reference genome but is relative to DeBruijn graph formed from reads as this is a de novo assembler).
Trying k == 50
started loadKmerCoverageMapFromReads
kmer "AATTGCAACTAAACATTAATATTATTTTTTTTCCTACACATTCAATTAAA" shares a hashCode with 0 other kmers
kmer "TAATTTCTGAAAAATGAACAAATAAATCATCTCCTCCGTCGTCAGGAGAT" shares a hashCode with 1 other kmers
Generated kmers from reads for k == 50, processing...
Removed 500000 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 616899)
Removed 443940 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 172959)
Removed 0 extremely low coverage kmers (kmerCoverageMap.keySet().size() == 172959)
estimatedGenomeSizeBp: 172959
estimatedAverageCoverageOfGenome: 14.123566856885157
started loadKmerCoverageMapFromReads
Starting to generate De Bruijn map [(k-1)mer to (k-1)mer edges] from kmers
Processed kmer 10000 out of 1116899
Processed kmer 20000 out of 1116899
Processed kmer 30000 out of 1116899
Processed kmer 40000 out of 1116899
Processed kmer 50000 out of 1116899
Processed kmer 60000 out of 1116899
Processed kmer 70000 out of 1116899
Processed kmer 80000 out of 1116899
Processed kmer 90000 out of 1116899
Processed kmer 100000 out of 1116899
Processed kmer 110000 out of 1116899
Processed kmer 120000 out of 1116899
Processed kmer 130000 out of 1116899
Processed kmer 140000 out of 1116899
Processed kmer 150000 out of 1116899
Processed kmer 160000 out of 1116899
Processed kmer 170000 out of 1116899
Processed kmer 180000 out of 1116899
Processed kmer 190000 out of 1116899
Processed kmer 200000 out of 1116899
Processed kmer 210000 out of 1116899
Processed kmer 220000 out of 1116899
Processed kmer 230000 out of 1116899
Processed kmer 240000 out of 1116899
Processed kmer 250000 out of 1116899
Processed kmer 260000 out of 1116899
Processed kmer 270000 out of 1116899
Processed kmer 280000 out of 1116899
Processed kmer 290000 out of 1116899
Processed kmer 300000 out of 1116899
Processed kmer 310000 out of 1116899
Processed kmer 320000 out of 1116899
Processed kmer 330000 out of 1116899
Processed kmer 340000 out of 1116899
Processed kmer 350000 out of 1116899
Processed kmer 360000 out of 1116899
Processed kmer 370000 out of 1116899
Processed kmer 380000 out of 1116899
Processed kmer 390000 out of 1116899
Processed kmer 400000 out of 1116899
Processed kmer 410000 out of 1116899
Processed kmer 420000 out of 1116899
Processed kmer 430000 out of 1116899
Processed kmer 440000 out of 1116899
Processed kmer 450000 out of 1116899
Processed kmer 460000 out of 1116899
Processed kmer 470000 out of 1116899
Processed kmer 480000 out of 1116899
Processed kmer 490000 out of 1116899
Processed kmer 500000 out of 1116899
Processed kmer 510000 out of 1116899
Processed kmer 520000 out of 1116899
Processed kmer 530000 out of 1116899
Processed kmer 540000 out of 1116899
Processed kmer 550000 out of 1116899
Processed kmer 560000 out of 1116899
Processed kmer 570000 out of 1116899
Processed kmer 580000 out of 1116899
Processed kmer 590000 out of 1116899
Processed kmer 600000 out of 1116899
Processed kmer 610000 out of 1116899
Processed kmer 620000 out of 1116899
Processed kmer 630000 out of 1116899
Processed kmer 640000 out of 1116899
Processed kmer 650000 out of 1116899
Processed kmer 660000 out of 1116899
Processed kmer 670000 out of 1116899
Processed kmer 680000 out of 1116899
Processed kmer 690000 out of 1116899
Processed kmer 700000 out of 1116899
Processed kmer 710000 out of 1116899
Processed kmer 720000 out of 1116899
Processed kmer 730000 out of 1116899
Processed kmer 740000 out of 1116899
Processed kmer 750000 out of 1116899
Processed kmer 760000 out of 1116899
Processed kmer 770000 out of 1116899
Processed kmer 780000 out of 1116899
Processed kmer 790000 out of 1116899
Processed kmer 800000 out of 1116899
Processed kmer 810000 out of 1116899
Processed kmer 820000 out of 1116899
Processed kmer 830000 out of 1116899
Processed kmer 840000 out of 1116899
Processed kmer 850000 out of 1116899
Processed kmer 860000 out of 1116899
Processed kmer 870000 out of 1116899
Processed kmer 880000 out of 1116899
Processed kmer 890000 out of 1116899
Processed kmer 900000 out of 1116899
Processed kmer 910000 out of 1116899
Processed kmer 920000 out of 1116899
Processed kmer 930000 out of 1116899
Processed kmer 940000 out of 1116899
Processed kmer 950000 out of 1116899
Processed kmer 960000 out of 1116899
Processed kmer 970000 out of 1116899
Processed kmer 980000 out of 1116899
Processed kmer 990000 out of 1116899
Processed kmer 1000000 out of 1116899
Processed kmer 1010000 out of 1116899
Processed kmer 1020000 out of 1116899
Processed kmer 1030000 out of 1116899
Processed kmer 1040000 out of 1116899
Processed kmer 1050000 out of 1116899
Processed kmer 1060000 out of 1116899
Processed kmer 1070000 out of 1116899
Processed kmer 1080000 out of 1116899
Processed kmer 1090000 out of 1116899
Processed kmer 1100000 out of 1116899
Processed kmer 1110000 out of 1116899
Setting up the graph...
processed 0 out of 1102027 nodes
processed 10000 out of 1102195 nodes
processed 20000 out of 1102364 nodes
processed 30000 out of 1102537 nodes
processed 40000 out of 1102686 nodes
processed 50000 out of 1102850 nodes
processed 60000 out of 1103014 nodes
processed 70000 out of 1103172 nodes
processed 80000 out of 1103345 nodes
processed 90000 out of 1103512 nodes
processed 100000 out of 1103674 nodes
processed 110000 out of 1103835 nodes
processed 120000 out of 1103994 nodes
processed 130000 out of 1104157 nodes
processed 140000 out of 1104303 nodes
processed 150000 out of 1104480 nodes
processed 160000 out of 1104652 nodes
processed 170000 out of 1104803 nodes
processed 180000 out of 1104984 nodes
processed 190000 out of 1105155 nodes
processed 200000 out of 1105308 nodes
processed 210000 out of 1105477 nodes
processed 220000 out of 1105632 nodes
processed 230000 out of 1105806 nodes
processed 240000 out of 1105963 nodes
processed 250000 out of 1106126 nodes
processed 260000 out of 1106296 nodes
processed 270000 out of 1106463 nodes
processed 280000 out of 1106625 nodes
processed 290000 out of 1106796 nodes
processed 300000 out of 1106947 nodes
processed 310000 out of 1107108 nodes
processed 320000 out of 1107281 nodes
processed 330000 out of 1107441 nodes
processed 340000 out of 1107610 nodes
processed 350000 out of 1107777 nodes
processed 360000 out of 1107952 nodes
processed 370000 out of 1108119 nodes
processed 380000 out of 1108292 nodes
processed 390000 out of 1108461 nodes
processed 400000 out of 1108619 nodes
processed 410000 out of 1108782 nodes
processed 420000 out of 1108963 nodes
processed 430000 out of 1109130 nodes
processed 440000 out of 1109293 nodes
processed 450000 out of 1109450 nodes
processed 460000 out of 1109623 nodes
processed 470000 out of 1109778 nodes
processed 480000 out of 1109942 nodes
processed 490000 out of 1110115 nodes
processed 500000 out of 1110268 nodes
processed 510000 out of 1110430 nodes
processed 520000 out of 1110580 nodes
processed 530000 out of 1110764 nodes
processed 540000 out of 1110920 nodes
processed 550000 out of 1111088 nodes
processed 560000 out of 1111245 nodes
processed 570000 out of 1111405 nodes
processed 580000 out of 1111562 nodes
processed 590000 out of 1111719 nodes
processed 600000 out of 1111871 nodes
processed 610000 out of 1112040 nodes
processed 620000 out of 1112204 nodes
processed 630000 out of 1112369 nodes
processed 640000 out of 1112545 nodes
processed 650000 out of 1112702 nodes
processed 660000 out of 1112872 nodes
processed 670000 out of 1113037 nodes
processed 680000 out of 1113201 nodes
processed 690000 out of 1113344 nodes
processed 700000 out of 1113504 nodes
processed 710000 out of 1113679 nodes
processed 720000 out of 1113849 nodes
processed 730000 out of 1114022 nodes
processed 740000 out of 1114184 nodes
processed 750000 out of 1114350 nodes
processed 760000 out of 1114517 nodes
processed 770000 out of 1114665 nodes
processed 780000 out of 1114823 nodes
processed 790000 out of 1114985 nodes
processed 800000 out of 1115148 nodes
processed 810000 out of 1115303 nodes
processed 820000 out of 1115469 nodes
processed 830000 out of 1115630 nodes
processed 840000 out of 1115805 nodes
processed 850000 out of 1115954 nodes
processed 860000 out of 1116100 nodes
processed 870000 out of 1116277 nodes
processed 880000 out of 1116436 nodes
processed 890000 out of 1116604 nodes
processed 900000 out of 1116781 nodes
processed 910000 out of 1116957 nodes
processed 920000 out of 1117120 nodes
processed 930000 out of 1117277 nodes
processed 940000 out of 1117442 nodes
processed 950000 out of 1117606 nodes
processed 960000 out of 1117780 nodes
processed 970000 out of 1117935 nodes
processed 980000 out of 1118095 nodes
processed 990000 out of 1118276 nodes
processed 1000000 out of 1118434 nodes
processed 1010000 out of 1118597 nodes
processed 1020000 out of 1118761 nodes
processed 1030000 out of 1118915 nodes
processed 1040000 out of 1119090 nodes
processed 1050000 out of 1119283 nodes
processed 1060000 out of 1119442 nodes
processed 1070000 out of 1119605 nodes
processed 1080000 out of 1119762 nodes
processed 1090000 out of 1119938 nodes
processed 1100000 out of 1120106 nodes
processed 1110000 out of 1120135 nodes
processed 1120000 out of 1120135 nodes
average read coverage: 3.529106929095648
Built DeBruijn graph (of k-1 length strings for nodes, kmer strings as edges) with 1120135 nodes and 1116899 edges
36308 tip nodes found
916274 tip nodes removed
1870 tip nodes found
4513 tip nodes removed
N50 dropping by more than 50% from peak, stopping search and using best k so far kWithBestEstimatedN50WithoutReference=45
Detected circularity in longest contig, removed 44 bases
>CONTIG1 (159662bp)
TAGCCGATTCCATAAATTTTTGTTAAACCAAAAACAATATTTTTTTTTTTTGAAATTAAAACTCCACAAATATTAATATTCATAATTGTTTTTGTTTATG
TTTTTTAATAACACAAAAAACAAAAATTTTTTTAAATCTTTTTACAATAACACATTTTTTACAAATTTTTTTTATTGAACTTTTTATTTTCATATAATTT
TATTAAATTTTTTATAATTTTTTTTTTTAATTTTTTTTTTAAAATTAATTTTCGGTTTACATTTAATAAAATTTGTTTGACCTCCTTCAAATAAGGTTGG
TATATTATAACCTGATCTTGTTTTTTGTCCTTTATGTCCTTTTCCACAAGTTTTTCCTTTTCCACTCGATAATCCTCGTCCTATTCTTTTTTTTTTTTTA
TTAAATTTATTTAAAATTAAAATCATTATTTTTTAATTAAAGACAATGCTTTAATAGTACAAATTACAACATTTATTGCATTTGTAGAACCATGTATTTT
AGAAAATACACTAAATACACCTATTATATTTAATACTTGTTTTACGTATTTACTTGCAATTAAACCAATAATATCATTTGCAGGATAAAGAAAAACTTTA
GTTGCACAATGTTTTGCAAAAAGAAAATAAGGAATAATATTTTTTGAAAAATTTAGTTTTAGTATATTTTTTTTAGCTTTTAAATAAGATTTTTTAATAG
ATTCCAATACTTCACGAGATTTACCTCTACCTAATCCTATATTTTTTTTACGATTACCAACAATAGAAGTTGCTGTATAGGAAAAAATTCTTCCTCCTTT
AACTACTTTTGTAACTCTATGTATTTTTATTAACTTTTTTTCCATAATAAGAATCAACAATTTTTTTTATATGACCTTTATATTTATAATTTTTAAAATG
AAATTCTTCATTATAAAAATTAGTAATTTTTTTAAATTCGTTTATTAAAATTCTTATTTTATTTTTGACGTTTTTTATTACAAAAATTATTTTATTATTT
TTTAAAAAACAAATACTTATAATTTTATTTGTTATATTTAGTTTAATTATCATTGTTTTTTTTTAGAAGATTTTTTAGTAAAACTTTTATTTAATAAAAA
AATACCTTTTTTTTTATATGGATCATATTTTTTAATTTTAGTTAATTTTTCAGTAAAAACACCTAATTTATCTTTATAAACTGATAATAAAATTATACTA
TTTTTGTTTATATTTACTTTAATATAATCTGGAATAATAATTTTTATAGAATTTGAATATCCAACAAACAAAATTAAGAAATTTTTTTCTAAGCTAATTT
TATATCCTATTCCAGATATTTCCAATTTTTTAGTCCATAAATTTTTGATTCCTAATATTAAATTATTTAAATTTGAAATAAATGCATTGAACATACATTT
CTGTCTTTTTAAATTTTCAACAATTAAAGTATCATTTTTAATTAAAATACTTATAAATTTTGGTATAAATAATTTACAAAAGTAAACATTTTTTAAATAT
AAAAAATCATTTTTCAGAAAATAATTTATATTATTTAATAATATCTTCATATTATTGAAAAAAAGATTTTGCCTCCTATTTTTAATTTTAAACATTCTCT
AATTGTTAGTAATCCAATATTTGTACTTATAATTAAACCATTATTAAATTCTATTTTTTTTAAATATTTATTTTTAATATATCTTTTATTTGATGGTTTT
GAAAATAATTTTATTAAAATTATAGTTTTAATTATTACAAATATTTTTTTTTTATTTTTAAATAATAAAACAAAAAAACTTTTAATTATTTTTTTTTCAA
TTAAAATTCCAAGTAATAACATGCTAATTTTCGAATAAATACAAATACAAATATTTTTACCAGAATTATAAGAATTTATCAATCTAATAATATTATTAAT
TAATGTTATCATTACCAACTAGATTTTTCTATTCCTGTAATATAACCAAAATTACCTATTTTTCTAATTAAATTTCTATTTAAACTAAATCTGTTATAAA
AACTTCTAGTTCTACCACTAATATAACATCTATTTATAAATTTAGTTTTTAACTTTTTTTTTTTAATTGATTGAATTTTTAACAATAAATTGAATGTAAA
TTTTTTTTTTAGCTTATTTTTAATAATAATAATTTTAAAATTTTTTTTACTATTTTTAAATATTAATTTATTTTTTTGAATTAAACTTTTTTTTGCCATA
TTAATTTATTGGAAAATTAATTAATTTATAAAATTCTAAAAAATTTTTTTTTAAATTTTTTAAAACTATACTAATATTAATACCAAATTTTTTTTCTGAA
AAAAAATCAGGAAAAACAGAACAATCATCAATTCCAAAATGAAAATTTCCAAATCTATCTAAAGATTTTATTTTAAATCCATTAAATTCTCTTATTTTTG
GAGCAACAATATTTATTAATTTATAATAAAAATTCCACATCAAATTTTTTCTTAAAGTAATATTTATACATGCTAAATCTCCTTTTTTTGATTTAAAATT
AGAAATTGGTTTTTTAATTTTAATTAAAATTGGTTTTTGTCCAGTTATAAAATTTAAAGATTCAAAAATATATTTTATATGTTTTTTATCATTGTTAAAT
TTTGCAATTCCTGATCTTATAATAATTTTATCAATATTTGGATTAAAAATATTTTTTTTTTTTGAAAATTCTTTAATTAATTTAAAATAATACATTTTCA
TATAATTTCATTTGATAAAGAAATTAATTTTGCAAATTTTTCATTTTTAAGTTCTCTTAAAATTATTCCAAACACTCTTGTAGAAATGATTTGTTCATTA
TTATTTAATAATATAACAGAATTATCATTAAATTTAATAATAGTTCCATCTAATCGTTTTATTCCAGTTTTGCTTCTAATTATCATTGCTTTTAAAACCT
GACTTTTTTTTATTTTAGATTTAAAATTGGCTGTTTTAACAACAACTTTAATTATTTCACCAATATTAGCATATTTTTTATTTGTTCCGTTTAAAACTTT
TATGCATTTAACAATTTTTCCACCACTATTATCTGCAACTTTTAAATAAGTTTGTTCTTGAATCATATTATTTTAACTTTTTCAATAATTTTAACTAATA
TCCAAAATATTATTTTTTTTATTGGACGAGACTTACGAATTAAAACATAATCACCAAAAAAGCATTCTTTAAACATATCGTACGCATTTATTTTTAAATA
ATAATTTATTCTTTTTTTATATTTCTTAATTAAATAATTTTTTTTAACTATAACAGTAATTATAAAATTATTAATTTTAATTATTCTACCTAATAAAAAA
TTTTGTTTATTCATAACATATTTCTGTTTTAATTGATAATTTTTGAGAAGCAACATAAAATGCTTTTTTGACAATTAAATAATTAAAATTTAATACTTCA
AACAATATTATTCCAGGTTTTATAACTGAAACCCAGTCATAAATTGGACCTTTTCCTTTTCCCATTCTTACTTCTAAAGGTTTTTTTGTTAAAACTTTAT
CAGGAAAAACTCTTATTAAAATTTTGCTACTTTTTCCTAAAAAATGAATTAGCATTTTTCTAGCAGATTCAATTTGTTTAAAATTTAAAAATCCATTTTC
AATAGATTTTAAACCATAATTACCGTGAATTAAATAATTTCTTTTCGTTGAAAATCCTTTATTTCTATTTTTTTGATATTTCGAATACTTTAATTGATCA
GGTTTATTTTTCATTATTTTTTTTAAATAACCATATTTTTACACCTAAAATTCCATATTGTGTTAAAGTTTCGCATTGATAATATTCTAAATTGTATTTT
ATTGTATGTAAAGGTATTCTACCTATTAAACTCCATTCTTTTCTCGCGATATCAACACCTTCTAACCTACCTGAAATTTGAATTTTACAACCAAAATTTT
TTCTATTTTTGAGTAATTCTTCTTTAATAATTTTTTTAATAGAATTTTTGTTTAAAATTTGATTAACAACATTATAAGCAATATTTTTAGCGTTTAAAAC
ATGATTAAATACAAAATTTAAAATTACATTCTTTTTTAAAATTTTAGAAATTTGAAAAACGAAAACATCTAAATAATTTTCAATAATATTAAATTGATCA
ACATTGTTAATATATAAATTTATTGTTAACTTATTTGATATGATAATATCAATATAACTTAAATTAATAAATAAAAAATTTCTTCTAATAATTTCTCTTA
TTAAAATATCACATTTTAAATAATAACAAAAATTTTTTTTAAAAGTATACCATAGACTGTGATATACAGTATTTTTTTTAAGCCTAAATAAAATTGGATT
AATTTTTTTTCCCATAAACTTTTAAAATTATATTAGTATACCTTTTTATAACAAAATCAATTTTTCCTTTTGCTCTGTAATTTATTTTTTTATAATATAT
TCCTTTATTTATTAAAAAATTTTTTATAAAAATTTTTGTTGAATTAATACAATTTATAATAGTAAAAGCTATTTTTTTTAAAATAAAATTTATTTTATTT
GTAAAATTTAAATAAAAAAAATAATCAATAGGAATATTAGATAATAATTTTGAATAACTATACATTTTTTTATATGAAATTGGTAAATTTTTAAAAAAAA
ATTTTTTAATCATTTTTTTTTTCGTTTAATTTTTTTTGAATGAGATTTAAATTTTCTGGTAAAAGAAAATTCACCCAGTTTGCAACCAATCATATCTTCT
ACAATAAAAACTTTTTTAAATAATTTACCATTATGAATATTAATATTTATTCCTATCATATCAGGTATTATAGTTGAATTTCTAGACCAAGTCTTTATAT
TTTCCAGTTTGTTATTTTTTATTTTTAAAAATAAACTTTTGTGTAAAAAAACACCTTTTTTAATTGATCTACTCATTTTTTTTTTGTTTTATAACCTTTA
GTTTTAATACCCCAAAGTGAACAAGGGTGTCTTCCACCACTAGTTTTTCCTTCACCTCCTCCATGTGGATGATCTACAGGATTCATTGCTACTCCTCTTA
CTGTAGGTCTTATTCCTCTCCATCTGTTTTGTCCTGCTTTATATAATTTTTTTTTAACATTAAAAGTACAGATTTTACCTATTGTTGCCATACAATTTAA
TGATATTTTTTTTTTGATTCCAGATGGTAGTTTTATTACTCCAAATTTATTATCTTTAAAAATAATTTCAGATTCTGAACCAGAAGATCTTGAAAAAATT
CCTCCAATTTTTGGAAAAATTTCAACACAATTAATTAATGTACCAATATTAATATTTTTAATTATAGTTGAATTACCAGAATTTAAACTTATTTTATTTA
AACTGGAAACAACAAAATTCGAAATTTCAAGATTGTTTGTTTTAATAATATATCTTTTTTCACCGTCTAAATACTTTATTAATGCTATTCTTGCATTTCT
ATTAGGATCATAATGAATTGATTCAACAATTCCTTTTATATTGTATTTGTCTCTCTTAAAATCTATTTTTCGATAAAATTTTTTATGTTTTTTACCAATA
TGTCTTACTGAAATTTTTCCTCGACTACGTCCACTTTTTTTTTTTTTAATTAAGCATAATGACTTTATTTTATTCTTCATTTTGTAAAATTGTTAATTAA
AAAATAAAATGAGTTTTTACTAAAAAAAATTTTTTCATATTTTAAAAAATCTGTTATTTTAATTTTTTTAATATTTTTAATTAAAAAATTGTTTTTATTT
TCCCAAAAAAAAATATTTTTTTTATCAAAAAAAAAAGTTAAAAAAAAAAAACAACTAATATCAAAAATAAATATTTTGTTACTATAAATTAACAAATTTA
TTAAATTAATTAAAACATTTTTTTGTTTTTTTTTTTTTTTTTTAAATAAATAATTACTAAATGTAACTCCTCCACTTCTCCAAATTGGAGATTTTTTATC
TCCAGCTCTAGCTCTACCTGTTCCTTTTTGAATCCAAGGTTTTTTACCAGATCCAAAAACCAGTCTTCTATTTTTTTTATAAGTTTTATTTTTGTTTTTT
AGATTTATAATAAAATTAAATAATATTATAATATTTATTTTTTTAAAACAATAATAAATTTTACAAGGAAAAAAATTTATAATAGGTAATATCATATTAT
ATCAACATATAAATAATCGTTTTTACTACCTGGTATTATTTTTTTAATAAAAAAAATATTTTTTTTTCTTAGAAAAACGTCTACAAAAAAAGTATTTTTC
TTATTTCCCATTTTTCCTGGCATTTTTTTTCCTTTAAAAACTCTTCCTGGATCTTGACATTGACCAACAGAACCTAATGTTCTGTATGAAAGAGAACATC
CATGGCTTTTATCTTTAGTTTTAAAATTCCATTTTTTAATAACCCCTGCAAATCCTTTTCCTTTACTAATACTAATTATTTTTAAAGTTTTAATTTTAAA
ATTTAAAATTTTTTTGTTAAATTTTGAATATTTAATTACATTTAACAAATAAAACTTATTATTAACATAAAAATTAAAACTACCTTTATTAATAAATATC
ATATAAAATTATAATTTATTTCTAAACCGTTAATAAATTTATTTTTTAATAAAAAAAAAATTAAAAATTTATCAAAATTTCTTATATATAAAATTGTTTT
ATAATATTGTATTTGAATTTGATCTCTTGCATGTTTATCTACATGTGGAGAAATTAAAAAAGTAAATTTTTCAATTTTTTTAGGTAAATTAAATGGACCA
ATAATTGAATACTTTTTTTTTATATTTTTTAAAAAAAAAAATAAATAAGATTTAATCTTAAAAAAAAAAAAAGATTTAATAGTAACTTTAATCATTTAAA
ACTTCAGTTATAATTCCAGCACCAACTGTTTTTCCTCCTTCTCTAATAGCAAATCTTAAACCTTTTTCAATAGCAATTGAAGATAATAATTTAACTATTA
ATTTTACATTATCTCCAGGCATTACCATTTCAATATTTTTCGGTAAATCGCATATTCCTGTAATATCGGTTGTTCTAAAATAAAATTGTGGTTTGTATCC
TTTAAAAAAAGGAGTATGTCTCCCACCTTCTTCTTTTGATAATATATATACTTCACAAATAAAATTAGTATGAGGTTTAATTGTGCCTGATTTAATTAAA
ACTTGTCCTCTTTCAACTTCTTCTCGTTTTATACTTCTTAATAAAATTCCTACATTTTCACCTGCAAATCCTTCATCTAATGTTTTTTTAAACATTTCAA
TACCAATAACAATAGTTTTTATTGTCTCTTTAAATCCCACTATTTCTATTTCTTCACCAGTTTTTATAATTCCTCTTTCTATTTTACCAGTAACTACTGT
CCCTCTTCCAGAAATTGAAAAAACATCTTCTATAGGCATTAAAAATGGTTTATCTATTATTCTATTTGGAACCGGTATATTTTTATCTAAAATTTCTAAT
AATTTTATTATTGAACTAGTTCCTAATTGATTATCATCTTTATTTTCTAAAGCTAATAAAGCAGACCCTATAATAATTTTAGTATTATTACCATCAAAAT
CATATTCTGTTAATAATTCTCTTATTTCCATTTCTACTAATTCTAATAATTCTTTATCTTTAACACAATCTGCTTTGTTTAAATAAACAATAATAGTAGG
AACACCTACTTGTCTTGCTAATAAAATATGTTCTCTAGTTTGAGGCATAGGGCCGTCAACAGCTGAACAAACTAGGATTGCACCATCCATTTGTGCGGCA
CCAGTAATCATATTTTTTATGTAATCTGCATGTCCAGGACAATCAACATGAGCATAATGCTTTGTTTCTGATTCATACTCAACATGTGAAGTTGAAATTG
TTATTCCTCTTTCTCTTTCTTCTGGAGCATTATCTATTGAATCAAATGGTCTACATTCACTCCCATATAAATCTGAAGAAACTTTAGTTAAAGCAGCAGT
CAAAGTAGTTTTTCCGTGATCAACATGTCCTATTGTACCAACATTTAAATGAATTTTTTCTCTATTAAATTTTTTTTTTGCCATAATTATTTTTTTTTTA
TTTTTTCTAATATATAATTAGGAGTTTCAGAATAATTGTGAAATTCCATATTATAATTAGCTCTTCCTTTTGTATTAGATCTTAAATCGGTAGAATATCC
GAATAATTCTCTTAAAGGAATTAAAGAATTGATAATTTTGAGATTATTATTATTATCTACTACGGAAATAATATTTCCACGTTTTTTACTAATATCACTT
ATAACTATGCCTAAATATTCTTTTGGAGATATAATTTCTACTTTCATAATAGGTTCTAGTAAAAAAGAATTTGCTTTTTTAAGTGCTTCTTTTAATGCAA
TAGAAGCAGCATTTTTAAATGCATATTCACTTGAATCGACAGGATGAAATGATCCATTTATTAATGTAATTTTTATTTTTGTTACTGGATAACCTAATAC
TACACCACATTTTATTTGTTCTAAAATTCCTTTTTCTATTGAAAGAAAATATTCTTTTGGAATAACGCCACCTACTACTTCTATTTTAAAAATAAAATCA
TCTTTTTCAATCAAAATAGGTTCTATTTTTAATACAACATGTCCGTATTGTCCACGTCCACCTGTTTGTTTTATATATTTTCCTTCTTGGATAATTGTTT
TTTTTATACTTTCTTTATAAGAAACTTGAGGTTTACTTGTTTTTGTTTTTATATTAAATTCATTATTGATTCTATCTATTATAATTTCTAAATGCAATTC
GCCCATTCCAGATAATATTAATTCTCCAGTATTTTCATTTATTTTAAATAATAAAGAAGGATCTTCTTTGCAAAATTTGTTAATTAAATTTAATAATTTT
TCATAATCGTTTTTTACAATTGGTTCAACTGAAACTGAAATGACTGGTAAAGGAATATTAATTTTTTCTAATAAAACTTTTTCATTATCAAAACTTAAAG
TATCACCAGTAAATGAATTCTTTAATCCGATTAATACAACTATATCTCCTGCTGAAGCAATATTTAAATCTTTTTTTGAATTTGCAAACATTCTAATAAT
TCTAAAAATTTTTTCTTTTATATTTTTTGAATTATTAAATATAATTTGTCCAGGTTCAATTTTACCAGAATAAATTCTTATAAATGACAATAGACCAAGA
TATGGATCATTAAAAACTTTAAACAATAATGCTAAAAATTTACTTTTAATATTTACAGAATAATTAATATTTGAAACGTTTTTAATACCAATATCTATAG
GAGATGGAAGAAAATTTACTATTGAATCAAGTAAAAATTCAATTCCTTTATTTTTTAAAGATGATCCGCATGCAATTGGAATTATTTTATTTAGAATTAC
TAATTTTCTAATTGATTCAATAATATCTTTAATTGAAAAATTTGAATTTATATATTTTTCTAAAAATATATCATCATATTCTGACAATGTTTCTAATAAA
ATATTTCTATACTTATTTGAAATGTCAAAATTTTTATTTGTAATATTTCTAATTTCTAATTGTGAATTATTCCAGATTAATTCTTTCATATTTATTAAGT
CTATAATTCCAGAAAAAGAATTTTCTATTCCTATATTTAAATTAATTATTAGTATATTGCAAAAAAATTTTTTTTTAATATTTTCAATTATACTTAAATA
TTTAGCTCCAATTCTATCTAATTTATTTACAAATAAAATTTTTGGAATATTAAATTTTTCAGATTGATTCCAAACAGTTTCAGTTTGTGGTTGAATACCA
GAAGATGCACATATTAAAATAACTGCTCCATCTAAAACTCTTAAAGATCTTTCTACTTCAATTGTAAAATCAACATGACCTGGAGTATCAATTAAATTAA
TTGATGAATTATAAAAATTAGTTTTCCAAAAAAAAGTTACTGATGCTGATGTAATTGTTATTCCACGTTCTTGTTCTTGTTTCATCCAATCTGTAATTGT
GTTACCAGTATGTACTTCTCCTATTTTGTGAGAAAAACCAGAAAAAAATAGTATTCTCTCTGTAGTAGTAGTTTTTCCTGCATCTACATGTGCTATTATT
CCTATATTTCTAATATTTTTAATATCATTCATATTTTAAAATTAGAGTATGCTTTATTTTGATCTATTATTTTATTTAATTCATCTTTTTGTTTTGTTGA
TAATGAATTATTATAATATGAATCTATTAATTCGCCTACTAATTTGTTTTTATACCCGTTTTCATTTCTTAGTTTAGAATTTTTTACTATCCATTTCATA
GAAAACATTAAACTTCTTTTTAAATTTATTTTAATTGGTATTTTATAAAATGATCCTCCAATTTTTTTTTTTTTAATTTCAAAATTTGGTTGAACATTAT
ACAAAATTTTTTTTATTAATATGAAAGGATTTTTATTTAGTCTAATAGAAATTATACTTATAGAATAATAAAAAATTTTCTGTGCTAAATTTTTTTTTCC
GTTATTCATAATATAATTAATAAATTTTGCAATAATATAGCTTCCAAACTTTGGATCATTTAGAATAACTGTTTTAAAATATCTTTTTTTTCTTGACATT
ATTTTTTTTTTCCATATTTTGATCTAGATGTTTTTCTATTTATAACACCTGATAAATCATATACATTTCTTATAATATGATATTTTACACCTGGTAAATC
TTTAACACGGCCACCACGTACTAATACATTAGAATGTTCTTGTAAATTATGACCTTCTCCTGGAATATATGCAGTTATTTCATTTTTATTTGATAATTTT
ACTCTGCAAACTTTTCTTAAAGCAGAATTTGGTTTCTTGGGAGTTGTTGTATACACTTTAATACAAATTCCTTTTTTTTGTGGAGAAGATAATAAAGCTG
GAGTTTTTTTTTTTTTAACAGATTTTTTTCTCTTAAATTTTAATATTTGATTTAATGTCATTTTAAAAATAATTTATATAAAATTTTTTTTTTCTATTTT
TAAAATCAAGATTATTTTTACTACTTTTAGATTTTTTTTTAAAGTATTTAATTAATCCTGTACCAGCAGGTATTAATCTTCCTACTAAAACATTTTCTTT
TAATCCTAATAAATAGTCTGTTCTATTTTTAATAGCAGAATTAACTAAAATTTTATTTGTTTCTTGAAAAGAAGCAGCAGAAAAAAACGATAAAGATTCT
AAAGATGTTTTAGTAATACCTGTAACTATTCTTTTATAATAAGACAATCTTCTTGAATTTACTAAAGTAGAAAAATTCTCGTTTATAGCATCTTCAAGAA
ACAAAATATCTCCTTGTTGACAAGAACTTTCTCCTGAAAAAATAATTTTTACTTTTTTTGTCATTTGTTTTAAAATTAATTCAATATGTTTACAATTAAC
ATAAATATTTTGCGGAAAATAAATTGAATTTATTTCGTTAACAAAATAAGAAAGTAAATAATTTATACTTATTAAGTTAATAATTTCATTTAAATCTGGT
TTACCATCAGACAAAATATCTCCTATTTTAACATAATCTCCATTATTAATATATAATTTTCTCAAACAGCTTAAAGTATATTCTTTATAAAAACCAAATT
TTGAAATAATAGTGATTATATAATTTAAATGAGAAAATTTAATTTTACAAACTCCATCAATTTCTGAAAGTAAAGCTTTCAATTTTGGTATTCTTGCTTC
GAACAATTCACTTAAACGAGGTAAACCACCAATTATACTTGATTTTAAAATAGTAACACTAACTAATTTTGCTATAATTTCTCCAGGAAAAATATAATCA
TAACAATAAACAACTATGTTAAAGTTTTTTGGAATAAAATAATTTTTTTTTATACTTTTATTAATAATTGTTATTTTAAAATTTTTTAAAGTTTTATAGT
AAAATTTATTGTCTTCTAAACAATAATTTTTAACAATATTATTTTCTTTTTTAAAATACACGTAACCAATATTTTCTGAATATATATAAAAATTATTATC
GTTGAAATTTATTAATTTAGTGTTTTTTTTAATTAAAAATCCATTTCTAAACTTTATTTTTGTACCATAATAAAATTTATATTTTTCCAATATAATTTTG
TTAAAAATCATAAATTCTCCATAAAGAGAAACAATAATTACTTCTCCAAATTTATTAATTACACATTTACATTTTTTAAATTTTACAAATCCAGAATTAC
TAATAATTAAATTATCGCTGTAAAAAAAATAACTTGCAACACCACCTGTATGAAATGTTCTCATTGTTAATTGGGTACCAGGTTCTCCTATTGATTGAGC
AGAAATAATACCAACAGAAACACCAATTAATACTATTTTATTAGTTGATAAGTCTATTCCATAACAAAAACAACAAACACCTCTATTTGAAATACAAAAT
AATACACTTCTAACAAAAATCGTATTTATTTTTTTTTTTATTAAAAGAAACAAAACTTTATTATCTATTAACGTATTTTCTTTTAATAATATTTTATTTT
TAAGTGATACAGAAAAAAGTAATAGTCTCCCGTATAAATTTTTATATAATTCAACTAAAGTAAAATAAGAATTTAAAAATATTTTTATTCCAGTTTTAGT
TTTACAATCAATTTTATAAATTACTAAATCATGAGATACATCAATTAATTTTCTTGTTAAATAACCAGAATTTGCAGTTTTAAGAGAAGTATCAGTTAGA
CCTTTTCTAGCACCAAAAGTAGAAATGAAATAGTTTTTCATATTTAAACCATTTTTTAAATTATCAAAAATCGGATCTTTGATAATATCTCCATTTGATT
TTGAAAAAAATCCTCTAAATGCAATTAATTGTTTAATTTGAAGCATAGAACCTCTAGATCCAGAATCTAACATTATAAATAAATTATTTGTTTTTTTTAA
ACCTGTTTTATTTATTTGTATTTTTTTAATTATTATTTTAACAAAAAAATTGTTTAATAATTCAACAATAGAAAAAATATTATTATTATATTTTAATTTA
AAAACAACTTTATTAATTAAAAAAGTGTTTTTAACTTCTATTAAATCAAAATAACTAATAGTTATTCCTGAATATGTTGATATTAAAAATCCTATTTTTT
TTAATTTTTCTAAAATAATAATTATCTTACTTAGTTGATTAAAATCGAATAAATATTTAATTAAATAAGTTAAAATTTTTTTTTTAAAAACTATATTTAA
AAAATAATCAAAATTTTTAAAAAACAAATAATTAATTATAATTCTACCTAAGGTTGTATTTATTTTTTTTTTATTATTAAATTTAATATTTATTTGAAAA
TTGTAATTATTAAATTCAAAAAAGTTTATAATATCGTTAATATCATTTAAACTAATACTTTTATTTTTAAAATCAAAAGTAAAATAGTATATTCCCATAA
TTATATCTTGAGTTGGAATAATAATTGCATTACCATTTGATGGCGAAATAATATTATTAATTGATAATAATAAATAATTAGATTCTACTTTTGCATTATT
TGTTAGTGGTAAATGAATTGCCATTTGATCTCCATCAAAATCTGCATTATAACATAAACAAACAAGTGGATGTATTTTAATAGTTTTATCTTCAGTTAAT
AAAATTTTAAAAGATTGAAAATTCATTCTATGCAAAGTAGGTGCTCTATTTAGCATTATGAAATAACTTTTAACTTTATTTTTTAATAATAAAATTATTT
GTTTTTTATTTTTTTTATAATAATCATCTATGAAACTTATAGTAGTAATTATTTTTTTTTTTTTTAGTTCATAATATAAAATTGGTTTAAATAATTCTAA
TGCTATATAAATTGGAATTTTACATTCATATAAAAACAAATTTGGTTCAACAGTAATTACTGTTCTTCCTGAAAAATCAACTCTCTTACCTAAAAGATTT
TGTCTAAATCTTCCATACTTTCCTTTTATTAAACTAGAAAAAGATTTTAATACTCTTTTTGAAGAAGTTAAGATAGGATTTATTAATTTTTCATTATCAA
ATAATGCATAAATAGATTCTTGTAATGAAATTCTTTCATTTATTAATAATTGTTTTGGATCTAAGTAAGTAGTCATTTTTTTTAATTTATAATTTCTTTC
TAAAACTATTTTGTATAATTCATTAATATCAGAAGAAGCAAATTTTTTTTCACTTAGAGCAATGAGCGGTCTCATTCTAGGAGGTAAAACTGGTATTTTT
TTTAGTACCATCCAAGATGGTTTATTTCCGGAATAATAAAAAATTAAAATTTTGTTAATATGTTTTAATATAGAAAATAATTTACTAAATGAATTGCATT
TTTTAATTTTAATTTTTAAAATTAAACAATCAATAAATAATTCACGATCAGATAATAATTTTTTTAAAGCTTTTGCACCAGACAATGTAAAAATTATATT
TTTTATTTTTTTATCAAATATAAATAAATTAAATTTTTTAATTTTTATGTCATATGACTTTATTACAACTTTTAATTTGAAATTAATTATTTTTTCTAAT
ATTTTATAACTCAAATTTAAAATATTACTTATATTACTATGCGAAGATTTAATATACCATATATGTACTGCTGGATAAAAAAGTTTTATATGTCCTACTT
TATATCCTGTGTGTGCTTTTTTTAAATTTAAAAATTTATCATTATGACAAAAACATTTTTTTTTTAAACAACATAATTCAAAATCGCTAAATATTTTTAA
ACAAAATAAACCATTAGGTTCTGGTTTTAAATTTTTAAAATTAATAAAATTAGAATGTGTAACTTCTCCAAATGACCAACTTAAAATTTTAGAAGAAGAA
GCCATTTTAATTGAAACTTTTTTAAAAATTGATTTTATCATAAAATTTTAATATCAAAACATAAAGATTGTATTTCTTTCATTAAAACTTGAAATGATTC
TGGAATTCCAGAATTAGCATCATTAATTCCTTTAATTATATTTTTGTATAACTCTATTCTACCAAGTATATCATCTGATTTAATCGTTAACATTTCTTTT
AACAAAAATGCGGCACCATATGCTTCTAAAGCCCAAACTTCCATTTCTCCTAATCGTTGTCCTCCAAGATTTGATTTTCCGCCTAAAGGCTGTTGTGTAA
CTATAGAATACGGACCTATAGATCTTGAATAAATTTTATCAATTACTAAATGGTTTAATTTCATAAAATATATATATCCAACATTTACCAATTGCAAATA
TCTTTTTCCTGTAATTCCATCAAATAATAATAATTCACAATTTTTATTTATTCCTATTGTCTTAATTATGTTATTAACTTTATTAACATTAAAATTATAA
AAATTATGAACGCATACATTTAATTGATTTTTAATATTTTTAAAAATTTTTAATACTAAAGAATTATTAAATATATTTAAATTAATATTTTTATCATATA
TACACTTAAAAATAATTTTTATAAATAATTTCATTTTGAAATAAGAAATTTTGTTTAAGTTTTTAATTTTTAAAAAAAAAGATTTAATTAAATTTAAAGA
TCCAGCTAAAAAAACTTCTAATAATTGTCCTACATTCATTCTAGAAGGAACACCAAGAGGATTTAAAATTAAGTCAATTTTATTTCCAAATTTATCGTAT
GGCATATTTTCATAATCAATAATGTTTGAAACAACTCCTTTATTTCCATGTCTACCAGACATTTTATCACCTATTGCAATTTGTTTTTTTACAGCAATTT
TAATTTTTATAATTCTTATTATTGAATTTTCAAAATCATCGTGTTTAATAAAATTTATTTTTTTATAAACAAATATATTTTTTTTTTTTAATAATTCATT
TGAAATTATGTTTTTAAATATATTTAATTTAAAATTAATTTTTTTATTAAAACATTTAATTTTAAATATATTGTTAATATTTATATTATAAGAATTAATT
CTTTTTTTTTTAATCGTTATCTTATTATTAAATAGTAATTTTTTAATTTTAGTTAAATAATAATTAAAAGTTTCATAAAAAAAGTTATTTATGTTTTTAC
AAGTATAATTAAGTTGTTCGAACTTTAACAATTTAAATATTTTATTTTTAAAATAAAATATTTTAAAATCGTTAACTGCAATTATTGTTCCTTTAATATT
TTTAGGAACTGTTAATGGTTGTTCATAATAATTAAAATTACTTTCTGAGAACACTATTTTAAATAATTTTTCTTCGGGAGAAAACTTTCTTTTTTTTTTT
GGTATCATTTTTCCAACAATTACATCTTTTGAAAAAACAAATTCACCAATTTTTATAATTCCGTTTTTTACTTTATTTTTTATTTTTTCGTTTGAACCAA
AACATTCATTAGAAACTATTTCAAAACCATTTTCGTTATACTTTAAAACAGTTATAAATTCGTATATATGAATAGAATTAAAATTATTTTTATTTAAAAT
TGAACTCGATAATAAAATTGAATCTTCAAAATTATATCCATACCAACTCATAAATGCTACTCTTAAATTTTTTCCTAAACTTATTTCTCCATTTTTAGTA
GAATTTGAATCTGCAATTATATTACCAACTTTAACAAAATCACCTTTTAATACTTTAGTATATTGATTTAATATTGTATTTTGATTTGTTCTAGTATATT
TTTCTAAAAAATATGTTTTTATAAAATTATTATTTTTAATAATTATTTTATAATTATCTGAATATAAAACATAACCATTTAAATCAGACAATATATTATA
GTTTGAATTTAGACCTATTTCTAATTCATTTCCAGTTCCAACAATTGGATTTTCAGAATCTATTAAAGGTACAGCTTGTCTTTGCATATTGCTTCCCATT
AAGCATCTATTCGCATCATTATGCGATAAAAAAGGAATCAAAGACGCACCAACAGAAATAATCTGATCTCCACAAATTTCAATTAAATCAATATTTAAAA
ATTTTTTATAATAATAATAATTAGTTTTTCTTGCTTCAAAGTATGGAGATCTAAATATTTCACCATATATAGTTTCAATTGTAGAATTAAAATTTACAAT
AAATTTGTGTTTATCAAATTTATTATCTAAAAAAGTTATTCCAAGAATTTTTCCTAAAATTGATATTTTATAAATAGTAGATATAAAATTATATTTATTT
ACTTTAGATAAGTAAGCTAAAGAATTTATTAATCCAATATTGTGTCCTTCTGGAGTATCAATTGGACAAACTTTACAATAATGAGAATAATGAATATCTC
TAATATCAAATCCACAATTTTCTTTTTCTATTCCGATTCCGCTTATTAATGAAATTTTTCTATTATGAGATATTTCAGCAAGAGGATTGTTTTGATCTAA
AAATTGAGATAACTCATTATTACAAAAATAATCTTTTAAACCTATAGTAATTATATCACTATTAACAATAAAATCTAAATCTTTATACTTTTTAAAATTA
TTCATTTTATAATTAATAAATTTAATTACTTTTTTAAATAAAAAATCAAATTTTATAGATAATAATTTTCCACAATTTAATATTAATTTATTTTCTAAAT
TATCAAAATTATCGTTTTGAATATTAAATTTTATAAATTTTAAAATTTTTTTAATTATTTCTAAATAAATAGTTAATTTATTATTTTTTTTTAGAAATAA
TCTTTTAAATATTCTCTTTGCTCCTATTAGTGAATAAAAATATTTTTTTTTATAAAATAAATTTCTAAAAAAATTTATAATATATCCTTTTTTATATTGT
AATAATTCTTTATTTATTTTTAATTTAAAATTTTTATATATGTATATCAAATAATTTAAATTTATTATTTTTTTTAAAAAAATAATAAAAAAAAAATAAA
TAGGATAATTTAAAATCAAATTAACTTTATTAATATCAATTATATTAAATTTTAAAAAAGTATTTTTTTTTAAAAAATAATTTTTACCAAATAATTTTCC
TAAAATTTTTTTAATACTAATGCAAATATATTTATTATAAATTTTATAAAAATATTTTTTTTTATTTAATAAGAAAAAAATTCTTTTTTTTCCTCTTATT
ATTTTTATTTTTATCTTAAAAAAAAAAAAATTAAAAAAATATTTTTTATTGTATCCTAAACAAATTAAAAATACATTAATTTCAAAATTTATTTTTTTAT
CAAATACTATTAAAAAATTATTAGTTATAATAAATTCTAACCAACTACCTTTTAAAGGTATAATTATACATTTTTTTTTTTTTTTTTCGGTATAAAAATA
TATTCCATAAGATTTTGTAAATTGAGAAATTAAGATTCTGTCTATTCCATTAATAATAAAATTTCCTTTTTTTGTCATACTAGGAATATTACCTAGAAAT
ATTTTTTTATAAACATTTAAATTTTTATTTATTACGTATATACTTATATACATGTACAAAGTACTAAATAAATGTAAATTTCTAATTTTAGTATATTTTT
CAGTATTACATGGTTCTAATAATTCAATTTTTTTAAGTTTTATAAAAATATTTTTATTGTTACAAATAAGTGGAAAATATTGCTTAAGTACTTTTTTAAT
TGAAAATAAGGAATTATATTTAGATGACAAAAAACTATTATAAGATCCGATTTGATTGTGTAATAAATATGGTAAATTACATTTATGAAAAAATTTTTTT
TTAGAAAAACAAAATCTATTAAATTTTAAATAATTTATCATTTAATTTCTACTATACCTCCTGCTAATGTAATTTTTTTAACAATTTCATCGCATTCTTT
TTTATCAACATTTTTTTTTATTAAAGAATTTTTGTTATCAACTAAATTTTTTGACTCTTTTAAATTTAAATTTGTTATATCTTTAATTGTTTTAATTAAA
TTTAATTTATTATTTCCTATTTCTTTTAAATAAATATCAAATTTGTTTTTTTCAGCAATTTTTTCGTTGCTCACAGTTAAATTAAATTTTTTTTCCAAAT
TAGATATTAAAGTTGTAATTTCTACTAAATTTAATTTAGATATTGAATCTAAAACATTATTTACAATATTGTTTTCCATAATTTTTAATCAAATTAAGTA
ATTTAAGTATTATTATTTTAAAATTTTTTATTAATTTTAAAATTAAATTTTCTTTTGATAATTTTGAATTTTTTTCAATTATTAAATCATAAACTTGATA
ATTTTTATTTATTATTAGTAATGGTTTAATATCTAAAAGATCATTTATGTTGTTTATTAAGCTGTAAACTAAATCATTATTTTTATGAATTAATAAGTAT
AATGATTTTTTATAATTAAACTTAACGATTTCAAATATAATTTTGTTTGAACAAAAAAAAAATTTTGTTTTGTAATTTTTTTTTTTAATATTATCAATTA
AAAAAAATTTTATTTTTTCAAAATTATAAAACAAAATTATATAATTCTTATTAAGAAAATTTTTAATGTTTTTAATAATATTATAAGAGTAATTCATTTT
AATAAAAAACTTTTGCTCATTGTACTTGAAATATAAATTTTCTCAATAAAAATATTATTGGGAAAAAAAAAGTTTTTTTTTAAATTAAATATTAAAAATT
CATAATTTTTAATAATCATAGAATTATAAAAAACTGTTGAAGCAATTTTTATATTTAAAATATTATTTTTGTTTAAAACTAAATTAATTGTTTTATTATA
AAATTTATTTTTTTCAATTTCATAATTTTTATATAAACAATTAACAATTATTGATTTATTAATTTTGCCTATTTTTTCTTTTACTAAAAGTAAATTTTCA
CTATCTGTATATACTTTAGAAAAATTTACTTTTTTAAGTAAAAAATCATCTATAAAATTATTTCCATAATAAACATTATTTTCTATTTTTTGTATATTTG
ATAAAAATAATATTTTTTTCATTTCATTATAAGAATAAAATAAATTAGTGCTTAAATTAAAATATTTTTTTTTTTTATTTAATAATATTAAATTTAAATC
AATAGATTCAATAAAATTTTTTTTTTTTTTGTTTAAATAATTTAAAAAATCATTTATACTTTTTTTCATTAAAATATCCCATGCTAATTAAAGTTCCTAA
CATCATTTTTCTTAATGAATTTAATTTTTTTTTTTCATAAATATTTCTAAAATAAGATATTTCAACAACATCGTTTAAATTTATATTTGTTAATTGATTA
TATCCTGGTTTTTTAGAAAAACTTAAAATATTTTTTTTTTCTTTTATATAATCAGATAAAGACATTGTATTTATTAGTATTTGATAAGAATTATCATTAA
AAATTGCAATTTTAACGTTTATTTTTTCTAAACATATAAACTTAGTTAAATCATTAAATTTATTACAAAAATCAATTAAATTAATTCCATATTGTCCTAA
TATAGGTCCAATTCCACTTGAAGGAGTTGCTTTGTTTGGTTTTAATATTAGTTTTAATTTAAATTTTAAAATTTTCATTTTTTTTTTTTTTTTTAAAAAA
AATTGAAAACAAGCCAAAGAGGAATTGAACCTCTATTATTGATTTTGGAAATCAATACTTTACCGTTAAGTTATTGACTTTTTTTTTAAAAAAGAATTAA
AAAGCGAAAAAAGGGAATCGAACCCTTTACAATAACTTGGAAAGATAATATTTTACCAATAAACTATTTTCGCTATGTGAAGAAAAGAGGATTCGAACCT
CTGAAGCTTTCGCAGCAGATTTACAATCTGTTCCTGTTGACCTCTTAGGTATTTCTTCTAATTTTTATAATTAAAGTGTTTATTTAAAATTCGTTGTTTA
CAATACTTAAAGTATTTTAAAATAACGTTATTTATGTAAGTTATTAGTAATCCATATAATTCCATAAATTAATATATATATAATAAAGAATAAAAAAAAA
ATCGATAATATAAAAAATTTTTTTTTAAATCTATTAAATTTAAAAAAAAAAAATAATTGAAAAAATACTTGTATTATAAAATTAATTAATACTAACTTTG
CATTATAAAATATGTATATTAGATAATTAGTTAAAAAACTAAAGTATAAAAAAATTTTTTTTTCATTAATAAAAATTACTTTTTTATATAAATTATAAAA
GAAAGACATAACCAAATAATATGTATTATATGCCAAAATATTAAAAATATTATATTGATTATTTTGAATTTATATTTTAGTTTATTTAATAATATTAAAT
TTATACATATTATAATTGCAACAAAAACATGAAGTGCATGGAAAAATAAAATTAAATAATAATTTGACAAATAATTATTTAATTTAAAACTTAAATTTAA
AAAAAACAAGTGTTTAATATCTTTTAATTCTATTATTAAAAATAAAATTGAAAAAAAAATATTTAAAAAGTAATATTTTATATTTGATTTGTTTATTATT
TTTATAGTTAGAAAACTACTTAATAATAACAAAATAGTTTCTATTAAAATTATTCTATAATTATAAATTATATTTTTTAAAAAAAAAATATGTTTAGAAA
TTAAAAATGATAAAAAAATAGTTGAAAACATGATACAATCAGTAATAATATATAACCAAAAACCAATTACATTTTTATTAATGTTATTCATAAAAATAAA
TTATTATATAAAAAAATAATAATAAAAAAGTTTTTTTAAAAAAATAAAAATTTTTCCATATTAAAAAAAAAACTGATAAAAAAATAATTAATGAAATTAA
TAATGGGTTAATGCTTTTATTATGTGATGTTATAAAAAAAATATTTTTAATATACTTTTTTTTATAATAAAAGTCATTATAATTATTAATCAATAAATTA
TTAAAATACAAAAATTCAGAACTTCTATATATATATTGTTTATAAGAATGTTTATTTTTTTTAAAAAAAATAAAACAAATTAATTGAAAAATTACAGCGA
TTAAAATAAGTATACTTCCAATAAAACTAACAATTAATAAATTTTTCCATTTTATATTATAATATATATTAATTCGTCTAATCATACCACAATATCCTAA
AAAATAAAATGGAAAAAATGTTAAAATAACACCTATTAACCAAAATAAGATATTAAAATTAAATTTATTGTTTAATTCACAATTATTATATATTAAATTA
TACCATAAACTTATACCAGATAAATAACCAAATAATACACCTCCTATTATTACAGAATGAAAATGAGCTATTAAAAACATTGAATTATGGAAAATAAAAT
CTAAATTAGGAATAGATAGTATTACCCCAGAAAAACCACCAATTGAAAAAATTATTATAAATGATATAAACCAAAAAAAAAGAATATCTTTTTTAAAATT
AGAAAATAAAATTGTAAAAATCCAATTAAATATTTTTACACCTGTTGGAATTGCAATTATCATTGTTGATATACTAAAAAATATATTGGATAAATAACCA
GATCCCATTGTAAAAAAATGGTGTAACCACACTAAAAAAGATAAAACTGTAATTGATATTGTTGCGTATATTAATGAAATATAAGAAAAAACACTTTTTA
CACAAATATAAGATATTATTTCTGAAAATATTCCAAAAGAAGGAAGTATTAAAATATATACTTCCGGGTGTCCCCATGCCCAAATTAAATTAATGTATAA
CATTTGTTGGCCTCCATAAAACGATGAATAAAAATGAGTACTTAAGGTTCTATCTAAAAAAATTTGAAATAACAAAATTGTTAAAATTGGAAAAGAAATT
AATATTAAAATATTAGAACATAAGCAAGTCCACATAAAAATTGAAATTTTATTGAAAAAAAGTTTTTTTTTTCTTAATTTTAAAATAGTTGTTATAAAAT
TTATGGAACTAATAATTGTTCCTATTCCAGATATTTGAATACTCCAATTCCAATAATCTACTCCTACCCAAGGACTAAAAGTTTTTTCAGATAAAGGTGG
ATATCCTAACCATCCAGTTTTTGCAAACTCACCTATTATTAAAGAAATATTTATTAATAAAACTGAAAAAAATGTTAACCAAAAACTTAACAAATTTAAA
GATGGAAACGCAACATCTTTAGAGCCAATTTGTAATGGTATTATTATATTCATAATACCAATTATCATTGGCATAGCAACAAAAAATATCATTATGTCTC
CATGTGCAGTAAATATTTGACAATAATGATGATTAGATAAAAATCCTAAATGTTTATAGCATAAAATTTGTTGTAATCTTATTAGTAATGCATCTGTTAA
ACCTCTAAACATCATTAAAAATGATATTATTAAATAAATTTTACCAATTTTTTTATGATTAAATGTTGTTGAATAATTTTTAACTTTTTTTAAAAATTCA
ATTTTAAACGGTAAAATAAAAAAAAAAATAGTAAAAATTATAATAAATAAAAATGTTAAAACAACTTCATTTATTGGAAATATATTCCAATTAAAATTGT
TCATTTTTTTTTATAATTTAAAATTAAAAAAAATAATTTATTATTGTAAATACTAAAAAATTTTGAATTATCTATATAACTTTGTTTTGTAATAAAATTA
TAATTTTTAATATTGTAAAAGATTTTTTTTTTTTTTATATTAATTAACCAATAAAAAAAATTTTTTTTAATTAATGAAAATAAACATATTTTCATATGAG
AATATCCAATTCCATTGAAATTAGACGATAAACCATGCATTATTCCGTGTTTAATAGCTATTAAATTTAATGATTTTTTATAGTTAGACATACAATACAT
TTGAAGACCAATTTTTGGAATACATATTGTATTCATTAAACTATTAGAAATTATTAACATTTTAATAGGTATTAATAATGGAATACAAATTTCGTTTAAA
AGTAAAATTTTTTGTTTTGGATAAATTATAATCCATTTCCAGTTTGTTGAAACAATTTCTATAATTAATGGTTTTATATTTTTGTAAATATTATTCAACG
GATTTAAATAAATCGAGCTTTTTATAGTATATATTGATAAAAAAAAAATTATACTTGTAGGTAAAAACCAAACTAAAATTTCTAGACAATTTGAATTAAT
TAAATTTGGCAAAAATATTTGTTTTTTAAATAAAATAATTGTTACTATTGCAATTAAAAATATACTAATTATAAGTATTAAATATATTGTATTTAATAAA
ATTTTATTTTCTAATATTCGTTCTAATCCAATAAAATGTTTGTTAATAAAAATTATCATTTTTTATTCTTATTAATAAAATTTAATTCTAAATTAAAATT
TAAAAAAAATAGTTATGGTACCAATAGTAAAAATTTTTTTTTTTTAACAAAAAATGTTTTTTTTTTATCGTTTTTAAAACAAATAAAAAAAATAAAAATT
TTAACTATTAAAATATTTCAAAATAAAGTTTTTTTACCTGTTGCAAAAATAATTATTAATAAAAATAGTTTAGTTTTGTATAAAAAGAATATTGGTTTAA
TTTTATATTATCCTTGTATTTATAACATTTGTAAAATTATTTTTTATTTAAAAAAAATTTTTTATAAAAAAAACAAAAATAAAAAGAAAAATACTTTTTT
TTTTTTAAAAATACATAAAAATTTAATTTCTAAAATTTTATATTTAAAAATAAAAATTAAGAACACAAATTTTATAAATCATAACGATGTATTTAAAGAT
AATATTCTTTTTTATGGTAATACGTTATCAAGTATAATTGATTTTAATAATTATAATTATTTGCATTTTAAAGTAGATTTTAAAATTTTAAAATTAGAAA
TGAATTACAATTATTCAATTAATAATATTTTATTATTAGAAAATTTTAATAATTTTAAAAAAAAGAATAATAATTTATATATATACATGTTTTATTTTTC
TAGAATTAAAAAATTGAAAATTTTTAAAAAGTCAAAAAGTTATAATTCTTATATAAAAAAAGTATTTTATGATTAAATTTTTTGAATATAAAAAAATAAA
ATTTAGATTTATAAATTTTTGGAATAATAATTGGAAGATGGGATTTTTTATAAAATTTTCAAATATTAAAATAAATGAATCTTCTACTTCAATTAATTAT
TCACAACAATGTTTTGAAGGTATAAAATGTTTTTTATATAATAAAAAAAAATATTTAAAGTTAAAAAAAAATTCGTTTAGATTTCAAAAATCTAACAAAA
AAATTTTAAGTCCTATTATATCAATTTTTAATTTTATAGTTACAATTAATAATATTACATATATTAATAAAAATTATATTCCAAACATTAAAGTTGGTTA
TTTATACATAAGACCTTTATTGCTTGGAATAGGAGGAATTCTAGGAGTTAAAAACTCTAAAAAATTTTGTTTTTTAATATATTGTTGTCCTGTTAAAAGT
GTTATTAATAATAAAGTTTTTATTAAAACATGTTTTTTAAAAAAAAATATTAAAAACCTTGGAAGTTTTAAACTAGGATGTAATTACATTACTAGCATTT
TTAATGATTATTATAACAAAATTAAAACTTTTGATGATATTATTTATTATGAAAATGAATATTTTGAAGAAGCAAGCACATCTAATTTAATATTATTTTA
CTATAAAAAATTAGTTACAAATTTAAATAAAAATATTTTACCTGGTACTAATAAAATTAATATAATTTATTTATGTAAAAAAAAAAAATATAATATATTG
TATAAAAAAATTAATTTTAATTTAATAAATAATTCTAAAATTGTAATTACTTGCGGAACTGCCGTTTTTATTAAAAATTTAAAAAGTATTTTGTTTAAAA
ATAAAATAATAAACTATAAAAATAATTTTTTAATTGAATATATTAATTTTTTTAAAAGTTAATTTTTAAAACATGAATATTTTAAGAGGATTTAAAGATT
TGCTATTTTATGAAAAAAGTAAAAAAACAACTATTATTAATTTTATTAGAAATGTGTATCTTAAAAAAAAATATTTCGAAATTGAAACACCAATTTTAGA
AAAATTTGAAATTTTTAACAAAAAAAATAATTTATTTATAAACGAAATTTACAAATTTTATAATTATAGTAAAAAATTAATTTGTTTAAGACCAGAAAAT
ACAACAAATTGTACTAGAATTTTAATAAAGTATAATAAAATTATAAAAGTATTTTATATTGGTTCAATGTTTAGATATGAAAAACCTCAATTAAATAAAA
TTAGACAATTTACACAATTTGGTTTTGAATTGTACAATAATAAAATTTCTGAAGAGCTAAATAGTATTAAATTAACAAATAAAGTTTTATTATTACATAA
TAATTACAAACTAGAATTAAATAATTTTATAAATTATAAGTTTAAAACAGTATATTTAAACATTATAATACATTTTTTAGAAAAAAAATTTTTAATTAAA
ATTGTTAGAAAGATATTAATATATAAAATTTTAGACAATTCAAAGTTTAATTTTAAAAAATTTGTATTTAATTATAAATTTATAAACATTAAAAACAAAG
ATAAATTAAATAAATTTTTTTTTTATTTAAAAACAAAAATTATATTTAATCCTAAATTATTAAGAGGAATAAATTATTATAGTAATTTAATTTTTGAATG
GAAAAATAATAATAATTCTGTTTGCGGGGGAGGAAGTTATTCTTATTATTTAAGCAAAATTTTAAATAAAATTAATTTTTCTTTTGGTTTATCAATTGGA
ATTGAAAGATTATATTTTAAGAAAAAATATAATAAAAAAAAAATTAAAATAAATTCAAATTTAAACATAAATACTAAAAATTTTTTTTTTGTTAATAAAT
TTGAAAATATAAAAGTTTTAAAAATCAATAATTTTATAATAAAATTTATCATAAAAAAAAAACAAATTATAATTAGTAAAAAAAGAATATTTTTTTTTTT
AAGTATTATTAATAATGAAAAACTTATTAAAACTTAAAAAAAATGTTTTTTTAAAACATTTTTTTTGCAAAATTTTTTATACAAATATTAAAACTATTTC
TTTATTATTTTTTAAAAAGAAAATAATTCAAAAAAACAAAAAATATAATTTTTTAAATTATTATTATTTAGTTAAATTGTTTTTTTTTTTTTATATTAAA
AATAATGTAAAAATAATTTTTTATTGGTTAAAAAAATATAATTTTAAGTTTTTTATTAAAAAAAAAAAAAATATAATAATTTATATTAATAAAAAAAACA
TTAATATTAAAAAAAAAAAAATAATATGTTTAAAATTTTAAATTTTAAAGATAAAAGAATTAGATTATTTTTTAAAAATGTTAAAGTAAGTTTTAACTAT
AATATACTTTATATAATAAAACAAATGATAATATTAATGTATAAAAATAACGGTATAGGAATTTCTAGTAATCAAATAAATTGTTTTAAAAATATTATAA
TTTGTGATGTAAATTTTAAAAAAAAAAAACCATTAATTATGATAAATCCAAAAATTTTAATTAATAATAAAAATCATACTTTAGGTATGGAAGGTTGTTT
ATCTATTAAAAATTTTTTAATCTCAGTATTAAGATTTGATAAAGTTTATATAAAATATTTCAATATTTATAATAAAAAAAAAAAAAAAATATTTAATGGT
ATTAAATCTAGATGTATTCAACACGAAATTGATCATTTAAATAGTAAACTAATATTAGATTATTCTAAAATTATATTTCAAAAATTATGAATTTTATTCA
AAAAAAAAAATTTCAAAAGAAAATAATTTATAAAAATGATTTTTTTTTAGTTTATTATAAATTTAAATTTTTTTATGCAAAAATTAAATGTTTATTTTAT
TTAAATTTGATAAATGTTAAACAAAAAAAAAAATTTTTAAATTTAATAACATATTGTAAAAAAAAAAATATAAATATTTTTAAAAAAATATATTTTTTTG
ATAAAAAAATTATTTTTCAGTTATTTTTATTAGATGAAATCGACAAAATTCATATTTTAATAATTTTTATAAGAAAAATTTTAATCAATTTATCAGAATT
AGAATATAAAACAATTTTTTATAAAAATGAAATAAAAAATTTTACAACAATAGGGCATTTATATTTATATTGGAATAATTTATTAAGTATTGACTGTAAT
TTATTGTTAAATTTAAAATTATTAACAAGTATTTTACCTTCTAACATTTTTAATTTTATTGAAATAAACAAAAATAAATTATTAGTTTTAATAAAAAATT
TATTAGATTTTAAAAGTATAAATGAAAATGTGCTAAACAAAAATTATGATTATTATTTATCATTTACATTTTTTTTAATTTCTATAAATACACATTTAAT
TTCTTTTTTTGAAGAAATAAATTATTTATTTTTAGAATATAAAAATGCAATTTTTATTTTTTTAAATTCTAATTTTTTTTTTTTTTATAATTTATTTAAA
AGCGATAGTTTGAAAATCGATAATAATATTAAAAATATTATAAATAGTTTAAAAAATTTTGAAATAATTAAACCTATTATTATAGAAAGTATTTATTTAA
TAAAGTCAAATTTAAATTTATTTAGAAAAATAATTTCATCAATAAAAATTAAAAAAAAAAAATTTTACTTATTAAACTTTAAAAAAGAAAAAATTTTTAA
AAAAATAAAAAAGTATCTTATTTCAAAAAAAATAAAATCTATTGATTGTAATATAATTATAAAAAAAATAAACAATTACATTACAATTAAAAATATTAAT
TTATATAATATTTCTTTTTATGAATTAAACAAAATTAATAAAAATATTAAAAAAGATTTTTTTGAAAAAATATCATTGGAAAACAAATTTAAAAAAAATA
ATTTTTATGGTAGCGATAATTATAAACAAATTTTAAAATCAATACAAAGAGCAAAATTTGCATTGAATAAAACTTTAATAAAATTTAAAATATGATAGAA
ACTAATAAAAATATTTATATAAATAATATAAATTTAAAAAAAATATTTAAAAATTTTTCTTTACCAATATATATATATGATTTTAAAAAAATACTAAAAA
ATTTATATTTAATTAAAAAAACAAAAATTTTTTGTTTTTATTCTATAAAAACAAATAATAATAAATTTTTATTAAAATTAATAAGTAATATAATTAAAAA
ATTCGACATTGTTTCAATTGAAGAATTTTTAAACATAATTTATATAAATAAAAAAAAAACGTATATAATATTTTCTGGTAGTGGAAAAAGTATATCGGAA
ATTTTAATTTCAATAAATTTAAATATTTTTAGTATTAATATAGAATCAATACAAGAATTATTTAAAATTTTTTTTTTTTGTAATAAACATAAAAAAAATA
TAAATTTAATGATACGATTAAACCCTAATATTGATTCATGTTCACATAATAAAATTTCAACTGGTAAAAATACTAATAAATTTGGAATTAATATTTCAAA
TATTAAACATTTTTTATTAATAATTAAAATTACTAAATTAAACTTTATTGGATATGATTTTCATATTGGATCACAAATATGTAGATTATCACCATTAAAA
AAGTTATTTAAAACAATAAAATTAATTAATGATTACAAAAAAGTTAATTTTGTTGATATTGGTGGAGGAAAAAATATAAATTATTATGAAAATAATATTA
ATATTAATTTTAATAATTATTATCGTTCAATAAAAAAATTAATCAAAAAATATAATTTAAGTTTTAAATTAATTATTGAATTAGGAAGGTTTTTTTTTGG
AAATTCTTGTATTATTTTATCAAAAGTTAATTACTTAAAATTTAACAAAATTTTTAAAATAGCAATTTTAAATGTTGGAATGAATGATATTATGAGACCT
GCAATATATAATTCATATCATAAAATTGAAAGTAATAATATTGGATACAAAATTAAAAATGTATTTGGCCCAATATGTGAATCAAGCGATAAATTTTATA
ATTCAAAATATATAAAAATTAATAGCAATAGTTTAATTATTGTTTATTCTTGTGGAAGTTATAGTAAAGTTATGAGTTCTAATTATAATAGTAAAAAAAA
AATTTTTGATATTGTTGTATATAAAAATAAAATTAAAATAGTTTATAAAAAAGAAAAATTTAATAATTTAATTAGCAATTATGTATAAATTTTTAAAACT
TCATTCTTGTGGAAATGACTTTTTAATTTTTTTTAAAAAAGTTAATGATTTTTTAATATCAAAATTTATAAACAAAAAATCTGGAATTACTTGTGATCAA
ATATTAATAATTAAAAAAATTATACTAGACAAAAATATAATAAAAATTGAAATAATAAATAACAATTTATCAAGAGCAAACAATTGTGGAAATGGAATAA
GATCATTATCTTGGTTTTTTTTAAAAAAAAAAAAAACAAATAGTATTATTAAATTTCCTGTAAAAAATAATTTTATTTTATCTTATAAGAACAAAAAAAA
AAATATTATATCTTTTTTCAAAATACCTAATTTTTTTAAAAAAATTATATTAAGAATAGAATTTTATTTTTTAAAATCAGGTTTTGTTGAATTAAATAAT
TTACATTTTATTACAATTATTAAAAATATTAAAACTTATTATTTATTTTTTTTTTATAATAAAATTAATTTTTTTTTTAACAATATTGTAAATATTGAAT
TTATTCAGATTGTTAAAAATAACGAATTTTATATTAGAATTTTTGAAAAAGGAGTAGGAGAAACATATTCTTGTGGAACAGGAATAATAAGTTCATGTTA
TTATATTAATAACATTAATAAAAAAATAAACTTATTTAATGTTTATTCTTTAGGTGGATTTAGTAAAATCAGTTTTTTTGATAAATTTATAATTTTGACT
AATAAAATTAATTTTTGTTGTTTTGGATATTTATGAATTTTATTTTAAATAATATAATTTATTTTTGGAAAATTAAAAAATTTTTATTTTTAAATAAAAA
TTCTAAAAAAATGGGAGCAGCAACTTATAATTTTAAAAATATTTGTTCCGTTTTAAATAAAAAAAAAATTAGTATATTATTTATACAAGAATGTTATAGA
GAATTTGATTCTTTTTTACCTAAAAGTAAAAAACTATTTATTCACAATCAAATTCAAGTTTTATGTAAACCTATTTTTTTTAACTTTGTTAATATTTATA
TTGATTCTATTAATTTTATTAATAAAAAAATTTTTTTAAAAAAAGATAATTGGAATTCTCCAATTTTAGGAGCAAAAGGTATTGGTTATGAATCTTCAAT
TAATAATTTAGAAATTTCACAAATTACTGTTTTTTATTTATTTGGTAATAAAAAATTATTTAAACCAATATTAGAAATAACATATGGATTAGAAAGAATT
TTTTTTTTATTTTACACAAGAGTTTTTTTTGATGAAAGATATTTTTTAATCAATAATTTATTTAAAATTAAAAATATTAAAATATTAATTTATGTTTATA
AAACTTTTTTTTTAAACTATAATAAAAAAAATTTTAATTTTTCTTATAAAATATTATTAAAAATTTCGAGTTTATTTAATATATTTGATAATTTTTATTA
CAACAATAATTATAATAGAATTAAAATTCTTGTTTTAATAAACAAAATATCAGAAAAAATAATAGAAAAAATATGAATAATATTTTAAAAAAAAAAATTT
TTTTTTTAAGTTTTATTTTTACAAAAAAAAAACTAAAAATACTTAAAAAAAGAAAAATAAAATTTTTTAAACTTGAAAAATTCTTTTTTAAAAAATTAAT
ACTAAGAAATATAATTTCACATTATTCTTATTTAAATATTAAAAATATTTTGTTATATTTTAAAAAACCAATAAAAAAAATTTTTAATATATTTTTAATA
AATAAAAAAAATGATATTAATTTAACAAAAAAATATTTTATTAAAATAATAAATTTTAAAAACAAACAATTTTATAATTATAAAAAAAGTATATTTAAAA
CTTATGTTTTAAGTTTATTTTTTTTAATAAATAAATATGATTTTATCAATTATTGGTTAATATTAAAAAATATTAAAAAAATTTTTTTAAACAAAATTAA
AAATAAAAATTATATTATATTTTATAATAGAATGTTAATTTTAAAATTATCAAACCTATAGAACACTATAGATCAGACTTATATTACTTTAACATTTGAA
GCTTGCAATCCTTTCTTTCCTTGAGTTACATCAAATGATACTCTTTGTCCATCTTGTAATGATTTAAAACCATCAACTCTAATTTCTGAAAAATGAACAA
ATAAATCATCTCCTCCGTCGTCAGGAGATATAAAACCAAATCCTTTTGTATCATTAAACCATTTTACTGTTCCGTTTGACATATTTTTTTTTAAAAAAAA
ATTTATTTTAAAAATTAATAACGTATTTTATACTATTCCATATAAATTTAGAATAGATTCTGAAATTATTACATCAGTATTGTAATCAGTTAAAATTATT
TTTTTTTCTTTTTTAGTGTAAAATATTAAAACTTCGTTTTTATAATAATTTATTATACTTAAAAATTTTTTCAAATTAACAAAAACAGTAAAATAACTTA
TATTCTTTATTAATATTTTTATATTAAAACTAAAAGTTTTTAATTCATCTAAAAATTTAAAATTTAAAAAATTCTTTTTAAGAGTAATTTTTAAAATAGT
GTTTTTATTAGATTTTATTTCTTTTATGTTATTAAATTTAAAATTTGTAACTATATAAAAATATTTTACTATAAAAAATTTAAAAATTAAATTTGAACTT
ATACTATTAAAATTAGTTTTTATCTTATTTTTAAAATACCATTCTTCATTTACTGAATTAAACAACAAAAATTTTTTAAAATTAAAAATATAAAAAAAAT
CTTCTGTTAATTTTAATTTTAAAAAAAAGAATGTTTTCTTATTAATAATTATTGATTTTAAAACTAAAATATTTTCGTTATAATTTATTAATAAAATGCT
ATTTTTTAAAAAAAAAAGTTCAATTAAAATTGTTAGTTCTAATAACTTTATATTATTATATTCGTATATTAACAATTTATTTAAAATTAATTTTTTGATA
AAAAAATAAAAATTTTTATAATTAAGTTTTAAAAAAAAATATTTATTATTGTTATAAAAAGTAATCACTTGATTCTTATTTTTCAAAATTAACTTATTCT
TATTAATAAAATAATTTACAAAATTATTAAAGTTAAAAAAAGTTAATATTTTTTTTAAATCGATTATAAACAAACTATTTATTTCAGTAATAAAATTTTG
ATTAAATAAATAAATTATATTTTTAAATTTTGTAAAAAATAAAATTTTTTTTTTATAAATTTCAATATAAATTTTTTTGTTACAGGAGTTTGATAAAAAT
AAAATTTTATTTAAAAATAATAAATTGTTTAAATTAACCATATGAATACTATATTTTCAAGAATAACACCATTAGGAAATGGTACGTTATGTGTTATAAG
AATTTCTGGAAAAAATGTAAAATTTTTAATACAAAAAATTGTAAAAAAAAATATAAAAGAAAAAATAGCTACTTTTTCTAAATTATTTTTAGATAAAGAA
TGTGTAGATTATGCAATGATTATTTTTTTTAAAAAACCAAATACGTTCACTGGAGAAGATATAATCGAATTTCATATTCACAATAATGAAACTATTGTAA
AAAAAATAATTAATTATTTATTATTAAATAAAGCAAGATTTGCAAAAGCTGGCGAATTTTTAGAAAGACGATATTTAAATGGAAAAATTTCTTTAATAGA
ATGCGAATTAATAAATAATAAAATTTTATATGATAATGAAAATATGTTTCAATTAACAAAAAATTCTGAAAAAAAAATATTTTTATGTATAATTAAAAAT
TTAAAATTTAAAATAAATTCTTTAATAATTTGTATTGAAATCGCAAATTTTAATTTTAGTTTTTTTTTTTTTAATGATTTTTTATTTATAAAATATACAT
TTAAAAAACTATTAAAACTTTTAAAAATATTAATTGATAAAATAACTGTTATAAATTATTTAAAAAAGAATTTCACAATAATGATATTAGGTAGAAGAAA
TGTAGGAAAGTCTACTTTATTTAATAAAATATGTGCACAATATGACTCGATTGTAACTAATATTCCTGGTACTACAAAAAATATTATATCAAAAAAAATA
AAAATTTTATCTAAAAAAATAAAAATGATGGATACAGCAGGATTAAAAATTAGAACTAAAAATTTAATTGAAAAAATTGGAATTATTAAAAATATAAATA
AAATTTATCAAGGAAATTTAATTTTGTATATGATTGATAAATTTAATATTAAAAATATATTTTTTAACATTCCAATAGATTTTATTGATAAAATTAAATT
AAATGAATTAATAATTTTAGTTAACAAATCAGATATTTTAGGAAAAGAAGAAGGAGTTTTTAAAATAAAAAATATATTAATAATTTTAATTTCTTCTAAA
AATGGAACTTTTATAAAAAATTTAAAATGTTTTATTAATAAAATCGTTGATAATAAAGATTTTTCTAAAAATAATTATTCTGATGTTAAAATTCTATTTA
ATAAATTTTCTTTTTTTTATAAAGAATTTTCATGTAACTATGATTTAGTGTTATCAAAATTAATTGATTTTCAAAAAAATATATTTAAATTAACAGGAAA
TTTTACTAATAAAAAAATAATAAATTCTTGTTTTAGAAATTTTTGTATTGGTAAATGAATATTTTTAATATAATTATTATTGGAGCAGGACATTCTGGTA
TAGAAGCAGCTATATCTGCATCTAAAATATGTAATAAAATAAAAATAATTACTTCAAATTTAGAAAACTTAGGTATAATGTCTTGTAATCCTTCAATAGG
AGGTATTGGAAAATCACATTTAGTTAAAGAATTAGAATTATTTGGTGGAATAATGCCAGAAGCATCTGATTATAGTAGAATACATTCTAAATTATTAAAT
TATAAAAAAGGAGAATCTGTTCATTCTTTAAGATATCAAATTGATAGAATTTTATATAAAAATTACATATTGAAAATTTTATTTTTAAAAAAAAATATTT
TAATAGAACAAAATGAAATAAATAAAATTATTAGATTTAAAAAAAAAATTTTAATCTTTAACAAATTAAAATTTTTTAATATAGCAAAAATTATTATTGT
TTGTGCTGGTACTTTTATTAATTCTAAAATATATATAGGCAAAAATATTAAAGCTTTGAACAAAGCAGAAAAAAAATCTATTTCTTATTCTTTTAAAAAA
ATAAATTTATTTATTTCAAAATTAAAAACAGGCACACCTCCAAGATTAGATTTAAATTATTTAAATTATAAAAAATTAAGTGTTCAATATAGTGATTATA
CTATTTCATATGGTAAAAATTTCAATTTTAATAATAACGTAAAATGCTTTATAACAAATACTGATAATAAAATTAATAACTTTATTAAAAAAAATATTAA
AAATTCATCTTTATTTAATTTAAAATTTAAATCTATAGGACCCAGATATTGTCCAAGTATTGAAGATAAAATTTTTAAATTTCCAAATAATAAAAATCAT
CAAATTTTTTTAGAGCCAGAAAGTTATTTTAGTAAAGAAATTTACGTTAATGGATTATCTAATTCATTATCTTATAATATTCAAAAAAAATTAATAAAAA
AAATTTTAGGAATTAAAAAAAGTTATATTATAAGATATGCGTATAATATTCAATATGATTATTTTGACCCTAGGTGTTTAAAAATTTCTTTAAATATTAA
ATTTGCTAATAATATATTTTTAGCAGGACAAATTAATGGTACAACTGGTTATGAAGAAGCTTCTTCACAAGGTTTTGTTGCAGGAATAAATTCCGCAAGA
AAAATTTTAAAACTACCTTTATGGAAACCAAAAAAATGGAATTCTTATATAGGAGTTTTATTGTATGACTTAACTAATTTTGGAATTCAAGAACCTTATA
GAATTTTTACTTCAAAATCAGACAATCGCTTATTTTTAAGATTTGATAATGCAATATTTAGATTAATAAATATTTCTTATTATTTAGGATGTTTACCTAT
TGTTAAATTTAAATATTATAATTCTTTAATATACAAATTTTACAAAAATTTAATTAATATTAGAAAAATAAAGTTATTTGATAATTTTTATTTGTTTAAG
TTAATAATTATAATGTCAAAATATTATGGTTATATTAAAAAAAAATATTTTAAATAATTTTCTTAATTTTAAAATAATTGATTTAAATTTAATAATATTA
TTATTATTTATACATTTAATTGTATTTTATTTATTAAAAAATAATAATTTAATGATATTATTATCAATATATTTAAACAATTTTATTAAAAATTCTATCA
ACCTAAATTCAAGAAATATAATTTTTTTTTTTTCACTAGTATTGTTTAATATAATATTATTTTCTAATTTTATTGATTTATTTCCAAATAATTTAATAAA
AAATTTTTTAAATTTAAAACAAATTGAAATTGTTCCAACTTCAAATATAAATATAACTTTTTGTTTTTCAATAATTTCTTTTTTAATAATTATAATGTTA
ACACATAAAAAAATAGGTTTTAAAAAGTATATATATAGTTTTTTTATTTATCCAATAAACACTGAATACTTATATTTATTTAATTTTATTATTGAAAGTA
TTTCTTATATAATGAAACCGATATCTTTATCTTTAAGATTATTTGGAAATATTTTTTCTTCTGAAATTATATTTAATATAATTAATAATATGAATGTATT
TATTAATAGTTTTTTAAATTTAATTTGGGGAATTTTTCATTTTATAATTTTACCTCTTCAATCTTTTATTTTTATTACATTGGTTATAATATATGTTTCA
CAAACTTTAAATCATTAAAAAAAAAAATGAATAATTTATTAATATTATCTTCATCAATAATGATAGGATTATCATCTATTGGAACAGGTATAGGATTTGG
AATTTTAGGAGGAAAACTTTTAGATTCCATATCAAGACAACCAGAATTAGATAATTTATTATTAACTAGAACTTTTTTAATGACAGGATTATTAGATGCT
ATTCCAATGATAAGCGTAGGTATAGGTTTATACTTAATATTTGTTTTATCAAATAAATAATATGAATTTCAATTATACTATTATTAATGAATTTGTATCT
TTTTTAATTTTTTTTTATGTTTCATTTAAAATTATATTTCCAGTTATATTAAAAAAAATAAATAATTTTTTAATAATTGATTATAAAAATTTTGTTTTTA
ACAATCAAGAAAAAATTATTAAAAAAAAATTATTAGATGAAATAGTTAAAAACGAAAATTTAACAAATAAGAAATTTATATCTTTAATAGAAAAAATAAA
AAAAAGTATTTTATTAGAAAAACAAAATTTTATTAATTTTATAAAATTAGAAAAAATAAACGTTCTAAAAATTTTTAAAAAAAAAATATTAAATAATAAT
ATGTTAATTATTAAAAACTTTTTAATTGAGATTAAAAAATTGTTTATAAATAGCTTTAAAAATATTTTTAATGAAATTATTTGTTATAACAATGAATTTA
TAATTAATTATGTTTAAATTTATAAACAGGTTTTTAAATTTAAAAAAAAGATATTTTTATATTTTTTTAATAAATTTTTTTTATTTTTTTAATAAATGTA
ATTTTATTAAAAAAAAAAAAATATATAAAAAAATAATTACTAAAAAATTTGAAAATTATTTATTAAAATTAATTATTCAAAAATATGCTAAATGAAGGAA
TAATAAACAAAATTTATGATAGTGTAGTTGAAGTTCTTGGATTGAAAAATGCTAAATATGGTGAAATGATTTTATTTAGTAAAAATATTAAAGGAATAGT
ATTCAGTTTAAACAAAAAAAATGTAAATATAATTATATTAAATAATTATAACGAGTTAACACAAGGAGAAAAATGTTATTGCACAAACAAAATATTTGAA
GTTCCTGTTGGAAAACAATTAATAGGTAGAATAATAAATTCTAGAGGAGAAACTCTCGATTTGTTACCAGAAATTAAAATAAATGAATTTTCACCTATTG
AAAAAATAGCACCAGGTGTTATGGATAGAGAAACAGTAAATGAGCCATTATTAACTGGAATAAAATCTATTGATTCAATGATTCCTATTGGAAAAGGACA
ACGAGAATTAATTATTGGTGATAGACAAACTGGAAAAACTACAATTTGTATTGATACTATTATTAATCAAAAAAATAAAAATATTATTTGTGTTTATGTT
TGTATAGGTCAAAAAATATCTTCTTTAATAAATATTATTAATAAGCTTAAAAAATTTAATTGCTTAGAATATACAATTATTGTAGCTTCAACTGCCTCAG
ATAGTGCAGCGGAGCAGTATATTGCTCCATATACTGGAAGCACAATAAGTGAATATTTTCGTGATAAAGGACAAGATTGCCTAATTGTTTATGATGATTT
AACAAAACATGCTTGGGCATATAGACAAATTTCTTTACTATTAAGACGTCCACCTGGTCGTGAAGCTTATCCTGGTGATGTATTTTATCTTCATTCAAGA
TTATTAGAAAGATCATCTAAAGTGAACAAATTTTTTGTAAATAAAAAATCTAATATTTTAAAAGCAGGTTCTTTAACTGCATTTCCTATAATTGAAACTT
TAGAAGGAGACGTAACTTCTTTTATTCCAACAAATGTTATTTCTATAACTGATGGTCAAATTTTTTTAGATACAAATTTATTTAATTCAGGAATTAGACC
ATCAATAAACGTTGGATTATCTGTTTCTAGAGTTGGTGGCGCTGCTCAATATAAAATTATTAAAAAATTAAGTGGAGACATTAGAATTATGTTAGCTCAG
TATAGAGAATTAGAAGCATTTTCTAAATTTTCATCCGATCTTGATAGTGAAACTAAAAATCAATTAATAATTGGAGAAAAAATAACAATATTAATGAAAC
AAAATATACATGATGTTTATGATATATTTGAATTAATATTAATATTATTGATAATTAAACATGATTTTTTTAGACTAATTCCAATAAACCAAGTTGAATA
TTTTGAAAATAAAATTATAAATTATTTAAGAAAAATTAAATTTAAAAATCAAATTGAAATTGACAACAAAAATTTAGAAAATTGTTTAAACGAATTAATA
AGTTTTTTTATATCAAACAGTATATTATGATTATTAAAGAAATAAATAGTAAAATAAAAATAACAACAAATATCAATAAATTAACTAATACTTTGAGTAT
GATTTCATTGTCTAAAATGAATAAATATATAAATTTAATTAATAATTTAGATTATATTAACATTGAATTAAAAAAAATTTTAGAATATATTATTATTAAC
ATTAAAAGTAACGTATTTTGTTTAATAATAATTACTTCAAACAAAGGATTGTGTGGAAATTTAAATAATGAAATTATTAAATACTCGCTTAATTATATTA
AAAACAATAAAAATTTAGATTTAATTTTAATAGGAAAAAAAGGAATAGATTTTTTTAATAAAAAAAATTTTTATATTAAAGAAAAAATAATTTTTAAAGA
CAATGAATTAAAAAATTTAGTTTTTAATAATAAAATTTTAAATGATTTAAAAAAATACGAAAATATTTTTTTTATTAGTTCAAAAATTATTAAAAATAAC
GTTAAAATAATAAAAACAGATTTGTATTTAAAAAAAAAATATAATTATTTAATAAAACATAATTTTAATTATGATTGTTTTTTAAAAAATTTTTATAATT
ATAATTTAAAATGTTTGTATTTAAATAACTTGTTTTGTGAATTAAAATCTAGAATGATTACAATGAAGTCTGCTGCTGATAATTCAAAAAAAATAATTAA
AGACATGAAATTAATAAAAAATAAAATTAGACAATTTAAAGTTACTCAAGATATGCTTGAAATAATAAATGGAAGTAATTTATGATAGGAAGAATTGTAC
AAATTTTAGGTTCTATAGTAGACGTTGAATTTAAAAAAAACAATATTCCATATATATATAATGCTTTATTTATTAAAGAATTTAATTTATATTTAGAAGT
TCAACAACAAATTGGAAATAATATTGTAAGAACTATAGCTTTAGGTAGTACCTATGGATTAAAAAGATATCTTTTAGTAATAGATACTAAAAAACCAATT
TTAACTCCTGTTGGAAATTGTACTTTAGGACGTATATTGAATGTTTTAGGTAATCCCATTGATAATAATGGTGAAATTATTTCAAACAAAAAAAAACCAA
TACATTGTTCACCGCCAAAATTTTCAGATCAAGTATTTTCAAATAATATATTAGAAACTGGAATAAAAGTAATAGATTTATTGTGTCCATTTTTAAGAGG
AGGAAAAATTGGTTTATTTGGTGGAGCAGGTGTTGGTAAAACTATAAATATGATGGAATTAATAAGAAATATTGCAATTGAACATAAAGGATGTTCTGTA
TTTATAGGAGTTGGTGAAAGAACTCGTGAAGGAAATGATTTTTATTATGAAATGAAAGAATCAAATGTATTAGACAAAGTTTCTTTAATATATGGTCAAA
TGAATGAACCTTCAGGTAATAGATTAAGAGTTGCATTAACTGGATTAAGTATAGCAGAAGAATTTAGAGAAATGGGTAAAGATGTACTTTTATTTATAGA
TAATATTTACAGATTTACGTTAGCAGGTACTGAAATTTCAGCATTATTGGGAAGAATGCCTTCAGCTGTTGGATATCAGCCTACTTTAGCAGAAGAAATG
GGAAAATTACAAGAAAGAATTTCTTCAACAAAAAATGGAAGTATTACTTCAGTACAAGCTATATACGTACCTGCTGATGATTTAACAGATCCATCTCCAA
GTACTACTTTTACTCATTTAGATTCTACTATTGTTTTGTCTAGACAAATAGCGGAATTAGGAATTTATCCTGCTATTGATCCATTAGAATCTTATTCTAA
ACAATTAGATCCTTATATAGTAGGAATTGAACATTATGAAATTGCTAATTCTGTAAAATTTTATTTACAAAAATATAAAGAATTAAAAGATACAATAGCT
ATTTTAGGAATGGACGAATTATCAGAAAATGATCAAATTATTGTTAAAAGAGCAAGAAAGTTGCAAAGATTTTTTTCTCAACCTTTTTTTGTTGGTGAAA
TATTTACAGGAATAAAAGGAGAATATGTAAATATAAAAGATACAATTCAATGTTTTAAAAATATTTTAAATGGTGAATTTGATAATATTAATGAAAAAAA
TTTTTATATGATAGGAAAAATATGAATTTATTAATTTTAAGTATAAAAAATATTATAGAATATAAAAATGCTTCTATATTAAATGTAAAAACATACTTAA
AACTTTTTTCAATTATGAATAATCATATAAATAATATTTGCGATGTTAATCAAATTAAGTTAATATTTAAAAATAAAATCATAAATATAAGAATTAATAA
TGGTTTTTTATTTCAAAAAAAAAATAATACTAAAATAATATGTAATTTTTATGAATTTTTATAATAAACATATATTAAATGATTTTTCTTTTAAAAAGTA
TGAAATTTTAACTTTATTTGAAATTAGTAAAAAAAAAATAAAAAATTTTTTAAATAATAAAAATATTTGTATTTTAAATGATAAAAAATCATTAAGAACA
ATTAATTCACTAATTAATAGTTTTAATTATTTAAATATTAAATATTTGCAAATTTTAAATAATCATAATATTAAAAAAGAAAGTTTTAAAGATTTTTCAA
GAACAATAGGTTTAAATTTTGATTATTTATATTATAGATGTTTAAATGACAAAATATTAAAAATTATTGCAAAATATTCAAGTTTAATAATTGTAAACTT
ATTAAGTAATGGATATCATCCAATTCAAGCATTAACTGATATTAATAGTTTTTTTTATAATAAAAAAGATGTTTTAATGTATATAGGAAATATAACTTCA
AATGTAATTAGATCAATAATTATATTATTATCAAAGATAAATTATCTTGTTGTTTTAATATCACCTATTAAATATTGGTTTAAATTTTTAATAAAAAAAA
TTTTTCCAAAAAAGAAAATACTTATAAGTGAAAAATTAATTTTATTTAAAAAAAAATATTATGTATATACAGATGTTTGGGAATCAATGAATAATAAAAA
TGTAAAAATAACTGATTTTTTAAACTTACAAATTAATAAAAAATTATTTGATTTAATTAAAATAAAAAAAGTATTACATTGTATGCCAAGATTTAATAAA
AGTTATTTAGATTTTGAAATTTCAAATTTAGTATTTGAATCAGATTACTTTTTAGTTAATAATTCGATAATTAAAAAAAATAAAATATTTAAAAGTTATA
TTTTTATTAGTAATTCATTTTTTTTTAAAATCATTTAGTTCTTTTAAATTAATATTATAAGATAGTTTGTTTATATAATCAAAAATTTCATTTTTTTTAT
ATTCAATAATTTTAATAATTTTTTTCATAAACTTAAAATATAATTTATTGCATGAAAATATCCATTCTCTTTCATACCTGAAATTACAACATTTGAATTA
AGTGAAATTATAGAATTTTTTCTATAATATTCTCTATTGTATATATTTGATATATGTAATTCGATAATTTTACCTTTAAAAATTTTTATACAATCTAATA
AAGCAATTGAATAATGACTATATGCACCTGGATTTATAATAATATAATTAAAGTTTATATTTTTTTGAATAAAATTAATTATTTTTCCTTCGCAATTTGA
ATTATAAAATTTAATATTTATAATATTTTTTGAGTATTTTAAAATTTTTTTTTTTAATTTTTTAAAAGAAATTTTAGAATAAATTTTTTCTCTTTTTTTT
AAAAAATTAATATTTGGTCCATTTATTATTAATACATTTATAATTTTATTACAAACAAACATAGTTTAATTAAAAATTTTTTGTTTAAATTAGTTTTTTT
TTTTAGTTCTAGTTCGTTACTAGAATATCCAAATTTTTTTATGTTTAACACATACGTAAAGTATTTTTTATATTTATACCAAAAATCATCATTTGAAGAT
TCAACAAAAATTATTTTTTTACAATTTAAAATTTTGTTTTTATATATTTTTTTTTGTTTATCAAATAATTTATTGCAGAATAATGAAATTATTTGAATAA
TAAAATATTTTTTTAAGAAAAAATAACATTCAAAACAGATTTCTAAATCTGATCCATTAGAAACAATTATTAATTCTATTTTTTTTTTATAAAAACATGA
ATAAGTACCAGTTATAATATTTTTAATATTATATATTTTAATAAAATTGTTTTTAAAATTTTGTCTTGATAAAATTAGAGACGAACAATTATTTAAAAAT
TTCAAAATTAATATCCAACATAGAATTAATTCTATATAATTATATGGTCTAAATATATAATTTCTTGGTATTATTCTAATTGAATGTAATTGTTCAATTG
GTTGATGTGATGGTCCATCTTCACCAACTAAAATTGAATCATGTGTAAATATAAAAATATTTTTAAGTTTAGATAAACAAAAATTTCTTATTGCACTATA
CATATAATTTGAAAAAACTAAAAAAGTAGAACAATAATTTATTCCTATTTTATCAGAAGATAACCCGTAATTTATTAATCCCATTGTAAATTCTCGTACT
CCATAATTTATATATCTATTTTTAAAATTTTTATATCTAATAGAATTAATAAAATTGTTTTTTGTTAAGTTAGAATTTGTTAAATCTGCGCTTCCTCCAA
ATGTTTCATTTATTGCATATATATTTTTTAATATATTAGAACAAACAAATCTAGTAGACTTATTTAAATTTATTTTATAGTATTTAAAATATAATTTTAA
AAAATTTATTTTTGGTATAATGTTATTAAAAATTCTTATTAACTCAAAAAAATATTTTTTATATTTTTTTTTGTAGTATATTAAATATTTTTTTTTATTA
TCAAAAAACATTTTTTTAACATAATCATATGTTAATGTAAAATTTTTTAAAATTTCTAAAAATTCAAATTTTGTAAAAATATTTCCATGAGAATTTTCAT
TATATGATTTACATGGAGAAATAAATCCTATTATAGTATTGTAAATTATAATTGTTGGAAAATAACTTTTTTTTGCTTTTAATAAAGATTTAATTATTGA
AAAATAGCAATGTCCATTTATTGGTCCAATAACATTCCAATTTAATGAAATAAATTTTAACTTAATATTTTCATTAAAATAATTTTTAACATTTCCATCT
ATTGAAATATTATTACTATCATATAATAATATAATATTGTTAATATTATAGCATCCACAAAAAGAACATGATTCGGAGGACACTCCTTCCATTAAACATC
CATCTCCACAAAATATCCAAACTTTATTATTGAATATATTAAAAAAATTATTAAATTTATTTTTATACTTTTTACTTTTTAAACCAATTCCAATTCCAAT
TCCAATTCCTTGTCCTAATGGACCAGTTGAAGCATCAATAAAATTTCCAATTTCAGGATGACCTGGTGTATTAGAATTAAACCTTCTAAAATTTATTAAA
TCTTTTATTTTATATACATTGTATAAATAAAGTAATACATAATTTATAATTATTCCATGCCCATTTGAAATTATAAGTTTATCTTTATTAATTGATTTTA
AATTGTTAAAATTTATTTTATAAAAATTTAAAAAAAAAATCGTAAATACATCACAAATTCCAAGAGGCATACCGGGATGTCCAGAATTAGCTTTTGAAAT
TGATTTAATACAAATTAATCTAATATTATTTATTATGTTATATAACATTTTAAAATTTAAAATTTTTTTTTTCAAAATTTATTCAATTTGTAATTATAAA
ACAAATACTTTTTCTATTTTAAATAAAAAAATAAAATATAATTTTTTTTTGAACTTTATTAATTATTATATAAATTATTTAAATTATAACAATAAAAAAA
AAATTGGAATTTTAATGTATTTTAAAGTATCAAAAGTAATTTCTTCTTTTAACATAGAAAAAAATGGTATCTTTTTTTTTTCAAACAAGAATGTTTTTTT
ATATAAAATATTAAAAAATTATGATATAAACAATATTTATCACGTAATTAAAATAATTAAAATAAATAAAATAAAGTTTAACTTAAAAATTTTAAAAAAA
ATATTTACAAAAATTTTAAAAAAAAAAAGAAAAGAAGTATATGAAAAATTAGAAGAAAGATATTTAATTACAATACTATTAAATAATTTAAACGAAACAA
AAAATAAGATTATTAATATTTATAAATCATTAATTAATTATAATACTAATAATTTTTTTTTAATTAATAAAGAATTTAACAAAGTATGTTCTTTACTGTA
TTTAAGTAAAAATGAAAGTTTGTCGAAAAAAATTCATTTAGGATTAATAAAAAATAATTTTAAAGAAGAAACTCCTTTTTATTTAAATTACATATTTAAT
TATTTCTTAAAATTTAATGAGCTAAAATTAACAATTTCAATTGAAATTTATAACTTAGATATTTTAAAAATAATTAAAACAATCAAAAAAAATAAAAAAA
TAAAAATTTTCATTAATGTTGGTATAAATGATTTATTTTTTGAAAAAATTTTTAAAAAAAAAAAAATAATTTTATTTAATTCGTTTAAAATAAAAAAAGA
ATATGGTTATTACGTACAAAATTTTTTTGATGAATATGTTGGATATGGATCATTTAGAAAAATGTATTTTAAAATATTTAAAAACAAAAATATTTTTAAG
ATAAAAATTTGTGCTAAATATTTTTTTTTAAAAATTTTAAAAACTAAAAATTTAAAAATTTATTTTTTAGATTCTTTAAACAGAAACAATTTAAATAAAC
ATATTAGTAATTTACTTACTGGATTTTTTCATCCAAAAATATTTGATAAAAATAATTTTTTTAAAAAAAAATATTTTTTTTACAAAAACAATAATATTTT
AATAAATAAAAATAATTCTTTTTATTTAGAAATAAAATTTTTTGTAAATTTTAAAATTTGTAAATATATTAAAAAAAAAATTGTTTTTTTATATAAATTT
TTTAACAAAGAAAGTGAAAATTATATTATAAAAAAAGAAATAAATTTTTGTTTAAATTATCGAATAAAACCAATAACAATTTATTTTCATGTAGTAAATA
AAAAAGTTGAAGAATATATTAATTTTTTAATTTTACAAATTAATTGTAATTTATCAAAGAAAAATAATTCATATTGTTGGTACTTTGGTAGTAATATTTA
TAATAGCAATTTTTTTTATATTAAAAAATATATATCAAAAAAATGGAATTTTATTATTAAGAAAATCATTTTATTTAAAATAAAAAATTCTGTTTATTTA
AATTTTAAAATTAAAAAAACAAATTTAAAACTAATATCATTAGATAATTTTTTATTAAAATTAATAATTAAAAATTGGCAAAAAAAAAATGAAAAATATT
AGTTTTGAAATATTTCCTTGTAATAACATTAAAGACTTATCTGTTTTAATAAATTATTTAAACAAAAATAAACCTAGTTTTGTTTCTGTAACATTTGGAA
AAATCAATAACTTAAAATTTGTTAAAAATATACAAAAACAGATTTCTACAAAAATAATACCACATTTAATATGTGATAATATATTTAATATTATTAATTA
TATAATTTATTTTATTAAAATAAAAATATTTAATTTTTTAATAATTACAGGAGACAAAAACAAAAATAATTCTATAAAATATATTTATTTTATTAGATTT
TTGTTTGGTCATATAATTAAGATAATAACAGGATGTTATTTTGAAAATCACAAATTTTCTAAAAATTTTAAAAACGAAATTTTATTTCATTATAAAAAAA
ATAAAATAGGAACTAATATGTGTATTACACAGTTTTTTTATAATTTTAACACAATAAAGTATTACATTAATATTATTAAAAAAACTGGTATTAGTAAAAA
TTTTATATTAGGAATAATTTCAAAAAAAAATATAAAAGATATTTTAAATTATACTAATTTATGTAAAATAGATATTCCAATTTGGATAATTAAAAATTAT
AAAGAATTTAATATTGAACTTTTTTTTGTTAAAAATTTAAAAAAATACAAAAATTTGCATTTTTATACTTTTAACAATATTAATTTAATTAAAAATTATT
TTAAATAAATTTTATTGTTATAAAATAAGTATACAAAATAATTAATAATAAAAAAAAATTTTTTATTAATAAAAAAAAAAATTTTTTTTATTAAAAAGTT
TCTAACAAAATTTAAAACATTTACTTTAATCATTTAAATTATTTTAAAAAAAAAAAAAATAAACAATTCATTATACTAAAAATAGTTAAAATTTAATTTT
TAAATTACTTTATTAAACTTGATATTTTTAAAAAAAAAAAAGAAAGAAAATGGGATTGTAGCTCAGTTGGTTAGAGCATACCTCTGATAATGGTAAGGTC
AATAGTTCAAATCTATTTAATCCCATTCAAAAAAAAATAAAATATAATAGATAAAATAAAACAATCATCTAAAAAACAATTTTCATTTATTTTATGAATT
GATTTGTTTATTAAACCAAATTCAAATATCTCAGAATTATAATAAGTATAAAAACAATATCTTAAATCAGATGTACCTCCTAAAAAATTTATTTTTGTTT
TAATATTTTGAATTGTATTAATAATATTTAAAATTTTTATTAAATAAAAGTTTTTATAACATAAAAACGGAATTCCACTTAAATTCCAATTTATATAAAC
ATTTTTTGATTTTAAAAAAAATAAAAATTTTTTTTTTAAAACAAAAATGCTATAAAAACTATTAAACCTTAAGTTTAGTTTTAAAAATAATTTATTAGAT
AACAAATTGTCTGTATTATTGCCACAAAAAATATTTGTATACTCTAAAATATTTTTTTTAATATCTTTCTTTAAAAAAAAAAAACTAAACAACATTTTTA
TTAAATTTTTACCTAAATAAGCACAATGTTTTATTTTTCCTAAAAAAATAATTTTTAAATTAAAAGATCCTCTTCTAGATATTTTAAGATAATCACCTAT
TATATTTTCAGAAGTTGGTTCTGTGCCAAAAAAAAAAAAAAATTTTTTTTTTCTTGCACTCAATAAACTAATACTATATTGTATTCCATAAATAGAAACA
GATTCTTCATCATTAGTTAAAATAAAATTTAGAATATTTTTATTAAAAAACAAAACACTAATGAAAGAACCTTTCATATCTATTATACCTCTATTTATTA
TTTTTTTTTTAAATAAATAACTTGAAAATGGATTTTTTAACCATTTAACATTACTTTCATGTACTGTATCTATGTGAGATATATACAAATAACTATAACT
TATATTAAATAAGTGTAAGTTTATTACTTTTCTAATTTTTATTAATTCTATATTTATCAATTTAAAGTATTTTATAATTATTAAAAAAATTTTTTTACAA
CAATATGATAAACTTTTAAATAAAATTTTTTTTAAAGAATTAAAAAATAAAATCATTAATTTTTATGTAATATAATGTTTGTTTTTTTTACAACAAAATT
TGTTTTTAAGCATTCAACTTTTCCATTTACAGAATTTCTAAAAAATATAATATTATTCTTAAAAGATAGTTTATTTGCTTTAATTAAAATTTTTTTATTA
TTAAATAAAAAAAAAATTTTAGTACCAGATGTAATATACAATCCAGCTTCTACAATACAATTATTTCCTAAACTAATTCCAATTCCTGAATTAGCTCCTA
ATAAACATTTATCTCCAATACTTATAATATTATTTCCACCTCCAGATAATGTTCCCATTATCGAAGATGACCCACCAATATCAGTATTATTAAATACTGA
TACTCCTGATGAAACTCTTCCTTCTATCATACAATTTTTTCCAATAAATGTATTAAAATTAACATATCCTTCTGACATTATTGTTGTTCCGGAGCATAGA
TACGCTCCTAATCTTACTCTATTTGTATTTGATATTCTTACATTATTTATCACTAAATATTTGCTTAAAAAAGGTATTTTATCTATACTTTCTATTATTA
ATTGATTTAGTGAAATATTTTTTTTATTAATTAAAATAAAAAAATCATTTAATTCATAAAATCCTATATTTGTCCAAACAATATTTTTAAGATTAATAAA
TAAGTTGTTTAAATTTAAACAATTAGGTTTTAAAATTTTTAAACTTATTATATTTAATTTTAATAATCCATCTATATTATCAACAATTTTATTATTTATC
TTTGTAATTTTTATTTGCAAATTAATAAATTTTATTTTTTTCCAAGTATTGTTAAAAATATTTAAAAATATTGCTTGTTTTATATTTAATAGTATCTTAT
TTTTAACAAATTTTATTTTTAACAATTTTAATTCTACTTTTTTTTTTATAAAATAAAATTTTGCTTCTAAAAAATTATTTTTAAAATTTGTTGAATAATT
ACATAAACTAAAAATTAACATAACTCTATTTTTCTTATTGTTGTTATTATTTTTTTTTTATTTAGTATAAACAAAGTATGTTCCCATTGAAATGTATAAA
AATTTTTTTTTGAAACAAAAATATTTAAAAAAAAAAATCCAGAACAATTTCCATAATTAAACATTGGTTCTATAGTAAAACTATCAAAATTTTTTATTTT
TTTTTTATTATTATTTACATTATGAAAAATAATATTTTTATTGTGTAACTTATCAAAAATACCATGACTACAATATTCTTCTGTTATATAAATTTTGCTA
TTTTTAATTTTATTAATTAAATAACCTATTTTAGAAAAAAAATTGTTTCTTTTTATATTTTTTATTAAATTTAAAAATATTTTTTTTAAAAAATTAAATT
TTTTGTTTTTAAATAAATTGATACAAGAATCACTATGCTTATTCTGATAATTTATTGCTATATCTATTTTAATGTTTAAAAAGTTTATTTTTTCAAATAC
TGGTATTCCATGACATACAATATTTGAAATTGATAAACAAGAACAAAATTTATAATTTTTAAAATTTATTGTGGAAGATTTTAAGTTTAATATTTTAATA
TAATAGAAAACAATATAATCTAATTCTGCTAATGAATAACAAAAATAATAATTACTAATAATAATAATTAATATATTATTCTTATTTCCACATTTTATAA
AATTTTTAATATGAAGATTATTCATTTAATTAATTCAGATATTTTTATAGGAAATAATAAAAAATTTTTAAATTTTAGTATGAAAAAATATATAATATTT
ACTATTAATAATTTAAATATTTTTGATCTTAAAAAAATTATAATTTCTGTTTTATTATTAAAAAAATATATAAAATTTATTTATAAAAAAAAAATGAAAA
TTTTGTTTATAGGAACAAAAAATTTTTTTAGAGAATTAATCTATAAATTTTCAAGAAGTATTAGACAACCTTTTGTATGTAATAAATGGATATCAGGAAG
TTTAACAAATTTACAAAATTATAAAAAAATGATAAATAAATTAAAAATTATAAGAAAAAAAATTAAGTTTAAAAGTTATACTAAAAAAGAAAAAATATCT
TTTTTAAAAAAAGAAAAAAAAATAGAAATTTTGTTTGGTGGGTTTAGAAATTTAAAAAAAACTCCAAAATTAATAATTATTAGTGATATAAATAAAGACA
AAATAATAGTAAATGAAGCAAAACGTTTAAAAATTAAGATTGCGTCTTTTTTAGATTCTAGCGATAATTGTTCAAAAATCGATTTTATATTACCTTGTAA
TAATAATTCAATAAATTCAATAAAAATAATTTTAAATATATTATTTAAAAATTTATGTTAAATTTAATTGTTCAGTTAAGAAAACAAACAGGAATAAGTA
TAAATTTGTGTAAAATTTTTTTACAAAAAAATAATTGGGATATTTTAAGTTGTATAAAATATATTGGAAAAATGAAAATTAGAAATAAAAATAATTATAA
TTATTACAGCATTATAAATTTTATTAATAATAATAAAATTTATGTATTAAAAATAAAATATAATTCAATTATAATTAATAATTCCAATATTATTATCAAT
TTAATTAATGATTTAAAAAAAAAAATAATTAAAAAAAAAATAATATTTAAAATTTTAAATATTTTAAGCGTTAAATTGAATGAAAATGTTTTTTTAAGTT
TTTATTTTTGTTTAGTAAATAAAAATTTAGAATTTTATAATCACAAGAATTTACTTTTTTGTTTGATAAAATTTAATAGCTTAAATTTTCTTAAAAAAAA
AGATTTATGTTTTCAAATTATATCTAAAAAAACGAAATTTATAAAATTTAAAAATTGTATTTTTAGTATTTCTAACCAAAAAAGTGTTGTAACAAATGAT
AAATTTGATAATTTATTAATTAATAAAATTGAATATTTTTTTTTTATAAATAAAAAAAAAATAAAATTTATATATGAGTAATTATCAAATAAAAATTAAA
AATTATTTAAATTTATTTAAAGAAAAAGCATTTTTTATAAATATTAAAAATTTAAACATAGAAATTTTTAATAATTTAAAAATTAAGTATAATAATTCTT
ATATTAATTTACGTGAATTATGTTCTATTAAAAATATTGATAAAAAAAAGTTTTTATTTATTTTTAATGATCAAAAAATCCTTCTTTATTTAATAAAAAG
TAAATATTTTGAAAATTTTGGATTCAATATTTTAAAAAAAAAAACAACAATAGAATTACTTGTTCCAAATATTAGTAGAGAATTCAGAATAAATTTTTTA
AAAATTATAAAGCAAGAATATGAATATTTTATAGAAATTTTAGAAAGTTTAAGAAAAAAAGAATTGTTACATATTAAAATTCAAAATATTTCAAAAGATG
AAATTTTAAGACAAGAAAAAGTAATTAAAAATGATTTTATAAATTATAAAAAGTTATTTAAAAATGAATTAGAAAATATATCAAATAAAATTTTTAATGA
TTAAATTTACTCACTTAAACGTACATACTGAATATTCAATAATAGACGGAATAATAAAAATTGGAAAATTGTTACAAATTTATAATAAATTAGGAATAAA
ATCAATTTCAATAACTGATATTTTAACAATTGCTTCTTTTCCAGAATACTATGAATATAGTTTAAAAAAGAAGATTAAACCAATAATAGGAACAGAATGT
TTTTTATTAATAAATGGAAAAATTTTAAATTTGATTTTAATCGCTAAAAATTATTTAGGATATTTAAAATTAATTAAAATTATTTCAAATGCGTGGAGAT
ACGGAAATATTGAAAATGGTGTTTTTTTAAAATTACACTGGTTGTATAAATTTAAAAATAAATTAATTGTAATAATAAATTTGCGATATTATTTATTAGA
AAGTTTTTTTAACTTTTGCGAATTTTCAATTTTTATAAAAGAGTTAACTATTAATTATAATAATAATTTTTATTTTGAAATTAATAGAATTAATCTGTCG
TTAGAAGAATTAATAAATAATAAAATTATTTATTTTTCTAAATTATTTAATATTAAATTAGTTGCTACTAATAGCGTAAAATTTATTTTTAAAAAAGATT
TTCCAATAAATTTAAGTAAAATTTTAATTTGTCAAAATAAAATTTTTGAAGAAAAATTATTTTTTGATTATTCAAATCAACAATACATTAAAACTTATAA
AGATATGAAAAAGTTATTTTTCGATATTTTAGAATCAATCGAAAATATTAATTTAATAATTTTAAATTGTAACGTATTATTTAATTTTTACAAATTTAAT
TTACCAAAAATAAAAATATCTAATTTTAAAATTAGAAAAAAAATTTTTGATAATCTAGTAAAAGCAGGATTAAAAAAACGTATTAAAAAAAAAAATAAAC
AAATAAAAATATATCTAAATAGAGTAAATAAAGAAGTGTTACTAACTAAAAAGCTAGGATTGATCGATTATTATTTGATAATAACAGAATTTATATTTTG
GACAAGAAAAAAAAATATTATTTCAGGACCAGGAAGAGGATCTGGTTCTTCTTCGTTACTATGTTATTCTTTATATATTACTGATATTGATCCAATTAAT
GAGAATTTATTGTTTGAAAGATTTTTTTCTTCTGAAAGATTAGGAATTCCAGATTTAGATTTGGATTTTTGTGTTTTAGAAAGAGATAACATTATTAGTC
ATTTATATAATTATTATGGATATAATAATATTTCGCAAATTGTTACATTTCATACTGTTTCTGCAAAATCTTCAATAAGAGATTTATCTAGAGCTATTGG
AATGGATTATATTTCGGGAGAAAGATTTTCTAGAAGCGTTCCATTTTCAATAGAATTATCTATGGAACATATTTTTAGAGAAAATATTTCAGTACGTAGT
TATATTTCTAAAAATCACAAATGTTTTGAAATTTGGAAAATTTCTTCAAAACTAGAAGGAATAGCGAGAAATATAAGTAAACATGCAGGAGGTGTAGTGA
TATGTAATACTGGATTAAATAATTTTACACCAATATTATTTGATGAAAATGAATGTATGACACATTACGAAAAAATAATATTGCAAGATATTGGTTTAAT
AAAATTTGATTTTTTAGGTTTAAAAACATTATCTACAATTAGTTTAACGTTAAAAATGATTTCTGAAAAAAATACTGGAGAATTTTTTATTGATGATTAC
CATACTTTTCAAATGATTAATAATTTAGATACTGAACTAATTTTTCAATTAGAATCATATGGGATTAAAAAAATAATTAAAAAATTACCTATTGAAAATA
TATTTGATCTAATAAATTTGTTATCATTGTATAGACCTGGACCTATACAATCAGGATTTATTAATGATTTTATTAATAGAAAAAATAATATTGTAAAAAC
ATATCATCCATATTCAGATTGTGATTTTTTACAATCAAAAATAATATTAGCAAATACTCATGGTATGATCTTATTTCAAGAACAAGTTTTGCAATTGATT
TTATTTTATACTAAATGTAGTAATTATGATTCAGAAAAAATATACGCTTCAATGATTAAAAAATCAAAAATTAAATTAAAAATTTCAAAGTTAATATTTA
TAAATGAATGTAATAAATTAGGAATTGATAAAAAAACATCATCTAAATTTTTTAATATTATTGAAAAATTTGCTTCATATAGTTTTAATAAAACTCATGC
TCATTCTTATTCAAAAATAGTTTATCAAACAGCGTATTTAAAATCTAATTATTTATTAGAATATTGTTTATCTAATATTTATGTTGATCAATTATTAGGA
ATTGATATTAATAATATTATAAATATTATTAAAAGTATTAGTGTATTTTTTTATAAACCAGATATTAACCTATCAGATGAAAATTTTAAAATTTATAAAA
AAGGTATTTTATACGGTTTTGATATAGTAACTTTTATTGATGAAAACTTTATTGATAAAGTAATATATTATAGAAACAAATTATTTTATTATAATAATTT
TGAAATGTTTTGTAAAATATTTAGTGTTTTTAAAATTAAAAAAAAAAAAATTATAGAAAATTTAATATTTTCAGGATTTTTCGATTGTTTTAAAATTAAT
AGAGTTATTTTATTTGTTAATTTTCAATTTATATTTGAAAATATTTTAACTTTAAATAATGAATATAGTAGAACTATTACTTATAAATTTGTTAGATATT
TTAATTATGCAAAAAAATTTTTTTTATTTAAAAAAATTTTAACAATATCTTCTATTAATATATTAAATATTGAAAAAAAAATACTAAAATTTTACACATC
TTTTTATCCTTTAGTTTTTTATTCTTTAAAATTAATAGGACACAAAAACTTTAATTTATTTAAAAGAATAGAATTAAACAATTTCAATATTTTAATTGCG
TATGGCAAATCAAAAAAAAATGAAAAAAAAAAAGCTTATTTTATAGGTTATAAAGAAACAATATTGAAATTTTATAAATTTTCTTTTAGAACCATTCTAC
CATGGAATAATATTAACGTTATTTTTTTTTTAAAATATGATTTTAAAAAAAATAAACATTTTATAATACATTGTTTTAGAATTAAACCTTTTTTAAAAAC
TATCGGGTCTATTTTAGTAATAGAAATAAATTGCAATATATTTTTTTTATATAATTTTGTTAAAAAAATATTTAATTTTTATAGCTTTTATGGCGAAAAG
ATTTATTTTTTAATTAAAATTAAAAAAAAAAAAAAAATACTTGATTTAAAAATAAATATTAATTTAAATGATATTTTTTTTCAGTATTTAAAAACAATAA
ATATAAAAAGAATTTTTTATTTAAATTCTTGTAAATGAAAAATTTATTTTTAATTAAAAATTTGTGCAATTTATTTTCTTTTGGTAAAGATAGTGTTTTT
TTCTTTTTTAAATGTTTTTATAAAAAAAAAAAAAACATATTTATTAATTATAATTTAATTTATTTTAAAAAATATTTTTTGCAATTTTTTAATCCAGAAA
TAAATATTTTTAAAAAATTTACTTTTAATGAATTAATATTTAGAAAAATAAAAATAAGTATTGTTAAAAATACATTTGTATACATGCATATAATGAATCA
TCATTTATTAGATAAAATTGAAACTATTTTTAATAAATTAATAAAAAAAAAAAAACAAATTGTTATAATAAACAACTGGTATTTTAAAAATAAATATGTT
ATGAAACCTTTAATTAATATTTTTTTAAAAAATTACAAGTTGAATATAAAAGATTCATCTAATAAATATATTTTTATTGAAAGAAATTTTTTAAGATGTC
TTATTAGTAAGTTTATTAAATAAGAAAAAAAATAAATTTTTATTTTTTTTTTTTTTTTTAAAAAAAAAATGAAAATTCATTTGTTAGGTTCACCTAATAT
TGGAAAATATAGAGAATTAAAATTTATTACAGAAAAATATTGGAAAATACAAAATAATATTAATTTATTAATTTTAAAAATGGAAATTAAAAAAATAAAA
ATGGAAAAAATATATTATCAAGTAAATAATAATTTTAATTATATTGGTTTTGGAGATTTTACATTATATGATAATATTTTAGATATAAGTTGTTTAATTA
ATACAATTAACAATAATTTAGATATTGTTAATATTAAAAAATACTTTTCAATAGCTAGAGGTATAGATAAGTTAAATATCAGTAAAATGACAAAATTTTT
TAATACTAACTATCATTATATTGTACCAAAAAACTTTAATAATTTAAAGATTGTAAATAATACATTATTCAAAGATATAAAAAATATTATTCAATTAGGT
TTAATTCCAAAATTAATTTTGTTTGGTCCTGTTTCTTTTTTATATTTATCAAATATTAAAATTGATAAATTAAAACAATTACTTGAAATATATTTATATA
TTTTAAAAAAAAACCTAAAATTAATAAATTTAACTATTCAAATTGATGAACCTATTTTAAGTTTAAAATTAAATAATTATTGGAAAAAAATATTTTTGTT
TTTTTACAAAAATATACAGAAATTAAACTTTAATTTGATTTTAACAACTTATTTTGAGTATATAAACAATTTAGAAATTTTAAATGATATTAAGAAATGT
ATTTTACACATTTCTCCTAAATATATTAACTTGATTAATAATCACAATAAATCATTTGGAATTATAAATTCTAATATTTTAAAAACAAATATTTTAGAAA
TTTTAAAAATTAAATATAAAAAAAATATTTTTTTTTCTTTAGTTGATAATAATAAATTATTGCCTTATGATATTTCAGTTGAAAAAAATAATTTAATTAA
AAAATTTTTTTCTTTTTTTTATCAAAAAATAACAGAATTAAAATTAATTAAAAATATTTATTTAAAAAAAATTAATTTTTTAGATATTTTATATTTAAAA
AATTATTCTATATTTAATGAAAAAATTATTGTTAAAAAAATTAATAAAAAAAATATAAAATGTAAGAATACAATTAATAAAAAATTATTAAATTATACAA
CAATTGGATCTTTTCCACAAAATAAAGAAATAAGAATATTAAGAAAATTTTTTAAAAAAAATGTTTTGTTAAAAAACGAGTATAAATTAATAATAAAAGA
ATATATATATATATTAGTTGTTAAACAAATAAGTTTAGAATTAAATTTATTAGTGAATGGAGAATTTGAAAGAACAGATATGGTAGAATATTTTGCTAAT
TCTATAAATGGAATGTATATTACAAATAATGGATGGATACAAAGTTATGGTACTAGATATGTTAAACCTCCAATTATTGTTGATATTAGAAACAGTTTTA
ATATTACTGAAGACTGGTTATACTTTTTTAAATACATAGTATCATTACCTAAAAAAGTAATTTTATCAGGACCAATAACAATTATTAAATGGTCATATTG
CATTAATGAAAAATATAAGTTTATTTTTTGTTATAAATTATCCGAATTGTTAAATTCTGAATTAATTAAATTACAATTATATGGTTTTAAAATTTTTCAA
ATTGATGAACCAACAATTAAAGAATGTTTACCTATAAATATTAAAAAATGGAAATTAGAAATAAATAATTTTTTATATTGTTTTAATAATAGTACAAAAA
ATATAAATAAACGAAATGAAATACATACACATATATGTTATTCAATTTTTGATAACATAATTAATATTATAAAAAAAATGAATATTAACGTTATTACTAT
AGAATCAACAAGAGAAAATATGAATAATTTAAATAATTTCAAAAATATTAATTTAAACATTGGTGGTGGTTTATACGATGTTCATTCTTCTATAATACCA
TACAAAAATGATATAAAAAAAAGAATAATAAAACATACGAAAATAATAAATTTAAATAAAATTTGGTTTAATCCTGATTGTGGATTAAAAACAAGAAATT
GGTATGAAATAATTTTTACATTAAACATTATTAAAAATGTGAAAAAAAAAATTTTAAATTATTATTCGTAATATTTTAAAAAAAAATAATTTTAATTAAA
GTTATATTAAAAAAAATAAATTTTAAAATTTAAAAAAAATAACGTTTTATTATTTATATTTATGATTTTTATACATTTTGTTAAAAATAATTTTTTTTTT
AAAAAAAAAATGAAAAAAATATATTTTTTATGCATTTGTGGCAAAAAATATTCGATACTATCTAATTGTAATAAAAATATAAATATTAATATTTGTAGTG
GATGTCATCCTTTTTTTACAAAAAAAAAAAACAATTTTAATAATTCTGAAAAAACAATTAAATTTAATAAAAAATATGAATTATTTTTTAAAAAATAATT
TTATTTTATTTGGAGATTCTTGTTCTGGTAAAACTTTTATTTATAAAAAAATAGAATTTTCAAAAATAGATATTGATTTTTTTTTATCATACAAAAATAT
TTTTTTTATAAATGAATTTTTTTTTAGATATTTTGAAAAAAAAATTTTAAAAAATTTTATTAAAAACAAAATAATAGTATTAGGAGGTGGATGTATACAA
TATTTAAAAAAAAAAAAAATTAAAAAAAGTATTTTAATTTTTAAAAATATTAATTTTATTAAATTTATAAAAATAAAAAATATTGAAAACAATAGACCAT
TATTAAAAAAAAAAAAATATATTAAAATTAGGTTTGTAATAAGAAAAAAAAAATATTCAAAAATAGCAAATTTAATTTTAAATAAATGTATAATTTGTAA
TCTCAGAAAAATATATGAAAATAATAAAATTTAAAAAAAAATCTTTATGTAGTTTATATATAACAAATGATTATAAAAAATTTATTAAAATAAGCAATAA
GCTAAAAATAATTATTACAGATTATAATGTTTATATTAACTTTTATAATCTAATTAAAGAAATTAAAAATTATAAATTTATAATATTTCCATGCGGTGAA
AATTTTAAAAACATTAATACGTTAAAAATAATTTGGAAATTTTTAATAAAGTTTTTTAACAAAAATATTTCTTTATTATCAATTGGTGGAGGAGTAATTA
ATGATATTGTTGGATTTATTTGTTCTGTATATTTAAGAGGAATTAATTTCATAGAAATTCCAACAACTTTATTATGTCAAATTGATTCTAGTATAGGAGG
AAAAAATGCTATAAATTTTTTTTCGAAAAATACAATTGGAACAATTAAAAATCCAATTTTTATATATTTAAATTATTCTATTATTTTTTATATGAATAAA
AATGATTTAAAAGATGGTTTTGCAGAAATTATAAAATATTTTTTATTAAATAACTTAAAATTTCTTTTTTATTTATATAAGATATTTAATTTTAAAAAAA
TATTGATAAGATCTTGCTATATTAAAATAAAAATAATATCTCAAGATTATTCCGAAAAATCTATAAGATCTGTTTTAAATCTAGGGCATACCTATGCTCA
TTGCATTGAAAATAACAAAATTAATAAAATTTCACATGGTAAATCAGTTTTAATAGGAATTATTTTTTCTTTATTTGTTTCTACATTATATTATAAAATA
GATTTATTTAAAATATTTAAAATTTTAAATTTATTTTTAATTTTTAAATTTAAAATAATTAATAAAATTAAATTTTCAGATGAAATGATTAAAAAAATTA
TATTAGATAAAAAATTTAATAAAAAAATTAATTTTATTTTGTTTAAAAAAATAAGTTGTTGTACTAAAAAAATAATAAAAAAAAAAAACTTATTATTATT
AATTATTTTTTTTTATGAAATCAAAATTAATAACAAAATGGCCAATAAAAGCAGCAAATAGAGCCATGTTGCGTGCAGTTGGATATAACGATTATGATTT
TAATAAATTTCAAGTTGGAATTGCATCAACTTGGAGTAATATAACACCTTGTAATAATCATATCAATGTATTAGCTAAAGCAGTTGAATATGGAGTAAAT
TCAAAATTTTGTAAAGGAACAATATTTAATACAATAACAGTATCTGATGGAATATCAAATGGTAATTTTGGAATGAAATATTCATTACTATCAAGAGATA
TAATTTCAAATTCAATAGAAGTTGTTGGCAAAGCACAAAATTTTGATGGTATTATTTCTATTGGAGGATGTGACAAAAATATTCCTGGTTGTATTATTGG
AATGTGTAATTTAGAAATTCCTTCAATTTTTATATATGGTGGTACAATTTTACCAGGCAAAAATAGAACAGATATTGTCTCTGTTTTTGAATCTTTGGGT
AAATTTTATAATAAAACTATTAATGAAAAAGAATTGTTAAATATTGAAAAAAACTCCATTATTGGATCTGGATCATGTAGTGGAATGTATACAGCAAATT
CAATGGCAATTGTTGCAGAATGCCTAGGAATTAGCTTACCAAACTCTTCTATTCAAAATGCGCAATCAATTAATAAAATTGTTAATTGTATTAATTCTGG
TAAATTAATAAAAATATTGTTAGAAAATAATATTACAATAAAAAAAATTATTAATAAAGAATCTATTTTAAATTCTATAAAAGTAATTTCATTATTAGGT
GGATCAACAAATTGTTTAATTCATTTATTAGCAATTGCTAATTGTTTAAAAATAAATTTATCTTTAAAAGATTTACAATTTGAAACTAATAACTTACCAA
CGTTATCTGATTTAAAACCAAGTGGAAAATTTTTTATTTCAGATTTAATTAATACTGGTGGAATACAAAAATTTTTAAAATATTTAATTGATATTAATTT
AATAAATGGTAATTTATTAACGGTTACAGGTAATTCATTAAAAGAAAATTTAAAATTTATTAAAATAAATTATAAAAATAAAATACTAAAAGATATCAAT
AATCCTGTTAAAAAAACTAATCAAATCAAAATATTATTTGGAAACTTATCAATAAATGGTTGTATATCAAAAATATCGGGTAAAGAAGGAGAAATATTCT
TTGGAAAAGCACTAGTATTTAATTCAGAAGAAGAATCTGTTAATTATATTTATAAAAAAAAAATTCAAAATAACACAATTATAATTATAAGGTATGAAGG
TCCAAAAGGAGGCCCTGGAATGAGAGAAATGTTAACTCCTACTTCTGCTTTAATAGGAGTTGGAAAGAAAAATTCTGTTGCTTTAATAACAGATGGTAGA
TTTTCAGGTGGAAGCCACGGATTTGTTGTAGGACATATTTCGCCAGAAGCATATGATTTTGGTACAATTTGTTTAATTAATAATAATGATCTAATAATAA
TAGATACTATAAACAATTTTATTAGTTTATTTTTAAACAAAAGAAGCTTAAAAAGTAGATTTTATAAAGTTAAAATATTTAACAAAATAGTTTCAGGTAT
TTTAAATTTATATAAAAAATATTCGATTTGCTCATCGAAAGGTGCATTATTAAATTATGAATAAAAATATTTTTATTATATTACCTTTAATATTTAATAT
ATACTTTTATAATAATGTTATAAATTCAAATTTGGTATTAAAATTTTCTAATAAACATACAAAAAATTTTATTAATTTTAAAATAATGTTTTGTCATTAT
ATAGGAAATATGACAGGTAGTTATTTTTTTTCTAAAAGAGTAAACGTTAACTTAGTTCCTATTAGTTTATATTTATTTTTATTTTTTTTAAATACTGTTT
TTTTTTTAGTAAAATGTTCTTATAATTTATATTTTTGTTTATTTTGCAAGATATTTATTGGTTTTTTAATTTCAAGTATAAATAATACTTCTGATTTTTA
TGTATCTAAATATAACAAAAAAATTTCTAATTTTTATAACTCTATAATTTATTTTTCTAGTTTTATAGCGCAAAACATACTAAATTTATATAAAATAAAC
TTTTTAAATAGCGATAATTTATATATATTTTTGTTTTTAATTACACATTTAAATAATAAACTATTATTAATTATTGATGATATTGATATTAAAAATAAAA
AAAAAATAAAATTTAAAAATTTAGAAAAACCAAATTTTTACATTATGATTGTTTTTATGTTTATAATTTCGTTTATATCAATATTAAACAATTATATCAA
TAATATGTTTAAAAAAAGTATTACAAATGAATTATCAACATTTACAAATTTTACTTCATTAGGAGGAGCATTATCTTTTTATATAATTACATTTTTTAAT
GATAAAGATAAAAAAAAATTATTAATTTTATTTTTAAGTTTATTATTAACTTTATCTAGTTTAATTACTTATTTTTTAAATTATAAAATAATTAAATGTT
ATATTTTATTTTTAATTGGTTTATTTACTTATCCTATATATTTTGTATCTTGTAGCGTTTTAAAAAAACGTTTTTCTAAAAAAAACAATATATTTTATTC
AATAGCTAATTCAATATCATATGCTATTTCTCCTTTTTTTTCTATTGTGTATAATAAAAATAATATTTTTATTTTTTTTACGTATATAAAAATTATTACG
ATATTATATTTAATTTTACTAATTTTATGTATAAAAAAATATGAATTTTAAAAATATTTTAAATTTTATAAAAATAGAATCTAAAAAACAAGAAAAAACA
TTAAATCTAATTGCATCTGAAAATTATTCTAGTATAACTAGTATTTTATATTCATCATCTTGTTTAACAAATAAGTATACTGAAGGATATCCAAATCAGA
GATATTATTCTGGATGCAAATTTTTTGATATTATAGAAAATAAAACTATTATTGAAACACAAAATTTATTTAATTCAAATTTTGCTAATGTTCAATCTCA
TTCTGGATCTCAAGCTAATTTTTCAGGTATACAATCTTTAATTAATAAAAATGAAAAAATTTTATCATTAGATTTAAAATCAGGAGGGCATTTAACACAT
GGATTTAGTAAAAATTTTTCTGGTAAATACTTTGATATTGTTAATTATTTGCTAGATAAAAATTTTTCTATAAATAAAGAATATTTATATAAAATAATTA
AAAAAGAAAAACCAAAAATTTTAATTTTAGGATATTCTTCATATCAAAAATATATAGATTGGGATTTTTTTTATTATTTATCTATTAAAAACAATTGTTT
TGTAATTTCTGACATAAGTCATATTTCTGGATTAATTGCCTCGGGCTTATACCCTTCACCTTTGAACTATTCTAGCTTAGTAACTACTACTACGCATAAA
ACATTAAGAGGAATTAAAGGAGGAATTATCTTAACACAAAATAGTAAAATTATTAAAAAAATAAATTTATCTGTTTTTCCTGGACAACAAGGAGGTTGTA
TATCAAACAATGTTTTAGGAAAATTAATTACTTTTAAAGAAGCTAATAATATTAATTTTTTAAATTATACAAAACAGATAATTATTAATTCAAAAATTAT
GTTAAAAACATTTTTATATAGAGGATATAAAACAATTGATTTAAAAACTGAAAATCATATGTTTATTATTAAAGTAAATAATAATAGCTTCTATTTAGAA
AAAAAATTAGAAAAATACGGAATTTTAATAAATAGAAATTTTATACCAAACGATAAAAATAAATCTTTAAATCCAAGTGGTATTAGAATTGGAACATCTT
GTATAACAACAAGAAAAATAAAAAAAAAAGGTTCAGAATTAATTTCAAATTATATTTGTGATTTAATTGAATTTAATAATAATATTATTAAAATTAAAAT
ACGTGTTTTGTGTCTTATTTTTCCTATTTATAAGTAATAAATTAAAAAAATATTATGTTATTTTTTTTTAATTAGTAGTGTTTAATTTAAAAACATAAAA
ATCTTATTAAATTTTTTTTACAATTAATAAAATTTTATTTTTTTTTTAAAAAAAATTATGAATAATTATAATTATTATAGTATAATTAATTCATATCAAA
TTGAAAATAAATATTGTGAAAATTTAATAAATGAAATTTCATATATTAATACAGAACCAAATTGGTTAAAATCTTTTAGAATTAATTCATTTAATATAAT
TAAAAAAATTAATATTCCAAAATGGGGAAATTTTTTTTTAAAAAAAATTATTTTAAGTAAATCTTGTTTTTATAATTTTCATGATAATAATAATAATTAT
TCATTAAAAAAAACTTTTAAACAAATAATAAATAAAAATATAGCAACAGATTTTGTTTACAATTCTATTTCAATAAAAACAACAATGAAGAAAAAATTAT
TAAAACACGGTATAGTTTTTTGTTCAATAAATGAAGCTATTAAAAATTATAGTGAATTATTAAGAAAGTATTTAGGTAGTATTGTAAAACCACAAGATAA
TTTTTTTAGTTGTTTAAATTCATCTATATTTAGCGATGGTACTTTTGTTTTTATTCCAAAAAATACAATTTGTCCAATAGAACTATCTTCATATTTTAGA
ATAAATGATGAAATTGGACAATTTGAAAGAACATTAATTATTTGTGATGATAATTCTCAATTATCTTATTTAGAAGGTTGTACCGCTTCAATAAAAAAAA
AACAACAACTACATTCTGCAGTTGTTGAATTAATAGCAAAAAAAAACTCTACTATAAAATATTCAACAATTCAAAATTGGTATGTTGGTAATAAATTTAA
TTATAATGGAATTTATAATTTTGTAACAAAAAGAGGTTTATGTTTTGGTAATAATTCATTAATTTTATGGATCCAAATAGAATCTGGTTCTTCAATTACA
TGGAAATATCCTTCATGTATTTTAAAAGGTAATTTTTCTAATAGTGAATTTTATTCAATAAGTATAACTAACAATTATCAACAAGTAGATACTGGAACAA
AAATGATTCATTTAGGAAATAAAAATTATAGCATTGTTAATGCAAAAAGTATTGCATTAGACTATTCTTTTCAAACATATAGAGGAATAATTAAAATTTT
AAATAATGCAAATTATTCAAAAAATTATACTTCATGTGATTCGATATTAATAGGTATGAGTAAAATTTATACATTTCCTTTAAATATCGTAAACAATAAA
TTTTGCAAGATTGAACATGAAGCTAGTGTATCACAAATATCATTTGATGAAATTAATTTGTTAAAAAGTAAAGGTATAAAAAAAAAAGATTGTTATAATA
TATTAATTAATAATTTTTGTTATGAGATTTTTAAAAAATTACCTTTAGAATTTAACAACGAAATTGAGAATTTAATATCTTCAATCATTAAATATTCTGT
GATATGATTAAAATTAAAAAATTATTTATTAAGTGTAATAATTTTTATATTTTTAGAAATTTAAATTTTTTTTTTTTAAAAAATAATATTTATGTATTAA
CAGGAAATAACGGAACAGGTAAATCTTCTTTTTTAAAATCGTTTGTTAACGATGAAAATTATTTTTTAACTGGAGAAATATATTTTCAAAAATATAATAT
CAAATTGTATGAATTAGATTTTGTTTCTCGTATTGGAATTTTTATTTCTTATCAAAATTCTATTGAAATTAAAAATATCAAAAATATTTTTTTTTTAAAA
ACTTGCTTTGAAATATTTAATTTTAACAAAAAATTTTTTTTCAAAAAATTAAAATGTTATATTAAATTGCTTTTTTATAAAAAAGACTTATTAAATAGAA
GTTATAATGTTGGATTTTCTGGTGGTGAAAAAAAAAAAAACGAATTTTTATTTTTATTAATTATTAATCCAATTTTAATTTTACTTGATGAAATCGATTC
TGGATTAGATCAAACATCTGTTTTAATAATTTTTAATTATTTAAACCTAATTAAAAAAAATAAATATATAATATTAATTTCACATAATAAAAATATAAAT
AATTTTTTGTTAATAGATTTTTATTTAAAAATTAAAAAAAATAAGATAAATATATTAAAATGTATATAAAAAAAATTTTAAATAGAAAAAAAATATTTTT
TGTTATTAAAAAATTTTTTCACTTAAATAAATATAAATTTTTATTGAATTTTAATTATTTAAAAAGGTTATTTTTTAAAAAAGTTAAATTAATTAATTTA
AAATGTTGTAACTTTAAAAAAAAAAAAAAAATAAATAATATAATTAATTATTTTTTTTTAAATACTTTATTTTTAAAAGTTAATTATGTTAATTATATAA
TAACAATTATAAATTTTAGTAATTTAAATATTTATTTTTTTATTAAAAATTGTTTTTTAAAGATAAAAATTATAAGTAAAAACAAAATTAATTTTTATAT
TTTTTTTAGTAAGATTATTTTTTTTTTTAAATATTCTTTTATTAAAATTATTTGTTTTTTATTTTGTAGCTTGTTTATTAAGTTTTTTTTTTTTTTAAAA
ATTAACAAATTATTTTTTACTTCTTACTTAATAAATCAAAAAAATTTATTATTTATTTATTACCTAAATTATAATAAACATAATATTGTAAAGTTTAAAA
TATTTTTTTTTCATTTTTCAAAATTTAATTTTTCTATTGAAAAATTTTTTTCATACAATAATAACAATATTTTTAATGTTAAAAGTAACATTTATATTAC
AAAAAAAAAAAATTTTATTTTTCAAAATTATGAAAATTATTTAAAAAAAAAAGTTAAAATAAATATTAATTTAAAAATTTTTTTAAACAAAAAAAATGTT
TTTTGTAAACATAAAATTATTTTTTACAAAAATATTAGCAAATATTATTTATTAGAAAAAAAATTCTTCTGTAATTTGTTATTATGAAAAAGTATTTTCC
TATATTTTATTATTATAAAAAAATAATCTATTTTGATAGTTCATCAACTAATCAAAAACCTAAAATATATTTTAAATCAATTTTAAATTGTGTTAAAAAA
AAAAATTTTAATATTAACAGAGGTGATAATTATCTGATTAAAAAAGTTAATTTAATAATTAAAAAATTTAAATTACTGATTAAAAAACTTATTTTAAATA
ATTATATAGAAGAAATAATTATATTTTTTAATTCTACTTATTCTATTAATTTTATTTTAAATAGTTTAATAAATTTAATTAATGAAAAAACAAATATATT
AATTTCTAACATGGAACACAATTCTTTGTTTCTTCCTATTTTAAAAATAATTAAGTTTAAAAATTGTAAATTAATAGTTTTTCCAACTATTAACAATAAC
TTATATTTAAATTTAATTTGTAACTATTTTAAAAATAAAATTTTTTTTTTTATAATTAATCAAATTTCAAATTTAGGATTAATGAATCAGATTAAAAAAA
TATCTACAATAGTTCATTACAATAATAGTATTATTATAGTTGATGCTACTCAATCTATTAGTTTTATAAATGTAAATATTAAGTACAGCAAAATAGATTT
TTTTTTTTTTTCATTACATAAAATTTTTTCTTCTACTGGTATTTCTATATTGTATTATAATTTGTTCTATTTAAATAAACTACTATTACCAAATTATGGT
AGTGGTTCAACTTTAAATATAAATTATAAAAAAATTAAAATTAAAAATTTTAATGAAAAATTTGAAATTGGTACACAAAATTTAATTGCTATTTTTTCCA
GTTATTATACTATAAAATGGTTTTTAAAAAATAAAAAATATTTTTTTTATATTAATAATTTTTTAAAAAATATTTTTTTTATTATTTTTAAAAAAAAAAA
AAAAAATATAATTAACAATTATTTAATTAAAATAAATAAAAATATTTTTTTTATATACTTTATCGAATTAAATAAAATTATAAGTAGATGCGATTTATTA
TGTAATTTTAGTAAACAATTATTTTTAAATAAAAATATTAATTGTAGAATTTCTATAAATTTTTATAATAATATTAAAGAAATAATTAAACTTAAATTTT
TAATTTTTTTTTTTAATTATAATATTACTAATAATTAACTCATTATATTTAATAAGACAATTTTGATTTTTTTCAATAAAATATATAAATTTATAACATA
TTAAATAACATTTTTTTTTGTTAAAAATTATTTTTGCATTTTTAAAAATGTTTTGTGAATTAATTTTTATTTTGCAAATTTTTGAATAATTAAAATTTAA
TTTTAAAATAAGTATTTTTGTTTTGCTCAAATAATTATACTTACTAATAATAATTTGATTATATTTAATTTTAACAATAAAATTTTTTTTGTATAATTGT
CCTAAGTTATAAAATAACTTGCAATAGTTTGTTTCAAAAATATTTTTATTTTCAAATACGTAATAATTACAAAAAAAAAAATTTTTATTAATTTTAAAAC
AAATTCCAGTAGTGCTTTTTTTGCATTTACAAACAAAATTTTTTTTTTTTGTAATATAAGTTATTATCAACTTATTAAAAAAACCAATATTTGAATGTAT
TTTTTTTTTAAAATTTAAAAAATAACATTGATCTTTTTTTTGATCTATTGAAACAAAAAAAAAATTTTTTTGTATTCTAAAATAATGTCCTGTTAATAAA
TTTAATTTAAAAATTATTTGCTTTACTAAGTTTATTTTAATAAATTTATTGCAATAAAAATCTAAATTATAGAAAATAAACTTAGAATTTTTAAAATCAA
AAAAAAAATTATATTCATTTTTTAAATTTATAATCAAAATTTTTTTTTTATTAATAATACAATTTAAAAATAAATATTTTTTGTCTAACAAATTTTCACA
AAAATTGTCAATTTTAATAAAAATATTATTATTATTAACTAAAATTGATGAATAATTTGAATCTTTACCACCAGAATTTAATAAAAAATTTAACATTAAA
ATCTATGAATAATTCTTCCTTTAAACAAGTTATACTGAGAAATTATAACTTTTACTCTATCTCCTATCATAATTTTTACATAATTTTTTCGTATTTTACC
AGAAATTGACGCAGTAATAATGCTTCCGTTATTTACTTTTACTTTAAATATACTATTTTGATTTGCATGGATAACAATACCTTCAATTTCAAAAAAAACA
TTATTATTATTATTCTCCATTTAGTAATAAAAGTCTTGTTCGTCTAATGGTAGGACCTCGCTCTTTCACGGCGGAAAAAAGGGTTCGAATCCCTTACGAG
ATAGCGAAAGTATCTTAATGGTAAAGTATCACCTTGCCATGGTGAAAGTTGCGAGTTCGAATCTCGTCTATCGCTAGAGATATGTCTGAGTGGATTAAAG
AATATGTTTGGAATACATATAAAACTAATGTTTTCATGGGTTCGAATCCCATTGTCTCTTAAAAATTTTGTTTAATTAGAAACATTAATTTTTTAATTGA
TAAAATACCAGAATAAAATTTTTTAACTTTATTATTTATTAGAATTGTTGGTATATTTTTAACGTTTAATTTTACTATTTCTTCTTGAAATAAGTCAACG
TTAATTTCAAAAAAAATTAAACTACTATTTTGTTCAATTTTTTTAATTTTTTTTTTTAAAATTTTACAAGGTTTACACCAATTAGCAGAAAACAATAAAA
TACTTTTTTTTTTAATTAATTTATTAATTAAATTGTTATTATTAACTAATATCATAAAAATTTAAGTAAATGATTACCAAAATAATTTCTTTGTAACTGA
ATATATTTAAAATTATAATTTTTGTTGATAATAATATTTATAAAATTAAAGCATGAAAACATTGTAAAACTTATAATTTCACAATTAATTATTAATAAAA
AAAATTCTTTATAATAATTAATATTTTTTTTTATTAAAATTAAAAACTTATTTATTAAAAAAAAATTTTTTTTAAAATTTAATAAATAATATATAATACT
AGAATTAATTATACAACTATCTAAATAAGACTTAAAAACATTATTATAATTATAATGCCAATTGTATTTTTTTTTAATTTTAATAATTTGATTAAATCCT
TGTAAATAACAAAAAAAAATACAAAAAAAAAATGTTTTTTTAACTTTAAATATTAAAATCGATATATTACTTTTTATAAAAAAATTATTGCTTTTATTTA
TTAAATTTCTTATATAAATATTATTAGATATAATTCTAATAATTAATGCTTCAAATATTGAATTAATATTTATGTAATTTTTAATGGAACTAATTACGCT
CCAACTTCCTGTTCCCTTTTGGTCTATTTTATCTAAAATTTTTATACTTTTTTTTATTAATATATTAATTAATATTTTTATTAAATATGAGTTCAGTTCC
GTTTTATTCCATATATTAATTATATTAATAATTTTTTTTTTTTTTTTTAATATTATTTTTAAAAAAAAATATATTTCAGATATTATTTGAAGAATTCCGT
ATTCAATTGCATTATGAATCATTTTCAAATAATGAGCACTACCAATACCGATTGAAATTGAATAACTACATCTTGAATATATATTTAAAGAAATAATATT
AAAAAAAAATAGTAATCTTTTTATCGTTAAAAAATTACCATCTATCATTAAACATAAACCTCTTAATGCACCTTCTGATCCTCCTGATATTCCTGCACTA
ATAAAACTAAATTTTTTTTTAATATTTAAAAAATTAAAATATGTATTTTTAAAATAAGAGTTACCAAAATCAATTAATATATCAGATATATTTAATTTAT
CTTTTATTAGAAATAAAATATTTTTGACAGGTAGTCCTGGCTTTATTAAAATAATTATAATTTTATAATTTGAAAAAGAATTAATAAATTTTTTCAAGTT
ATTAGTTATAATTTTAATATTAAAAAATTTTTTATTTAAATATATTTTTTCTCTATTATATACAGATAAAAAAATTTTTTTTTTAATAAGATTCAACGAT
ATATTTTTTCCCATTGATCCAAATCCTATGATTCCAATATGATTTATTAACATTTTTAAAAAAAACATCTAAAATATTTTTTATACAAACATAATTTGTA
TTAAAATACTTTTTTAAATATTTTAATGAAAAATAATTTTAAAAATTGAAAAATTTTAAAATTTAAAAAAAAAAATTATTTTTTAAAAAAAAATTACAGT
TAAAATAATATAATTTCAAAATTATCGATTTTTTTTTAAAAAAAAAAAAATATCTTTTTATGATTCAGATATTTCAAAACTCAATAAGTTTAGTTTTCAT
AATATTAATGGTATATCATTGTAATTTTTATTGTTAAATTGATAACAATATATTTTTTTTTTTAAATACATGAAAAAAAAAAATAAAAGATTTAATATTT
AAAAATCCTTGGTAGCTCAGCTGGTAGAGCAAATGACTGTTAATCATTTGGTCACAGGTTCAATCCCTGTCTAAGGAGTAGTTATATCCTGATTTGAACA
GGAGACCCAATCATCATGAATGATTTGCTCTACCAACTGAGCTATATAACCTAGTATTTTTAATTCAAAGTATTTGTTATATATTTTTTAAAAAATATAA
TTTTATTTATTTTAAATTTAAAGTTAAAGTCTTTGTATTTAAATTTTTTTTTTATTTTAAATATGTTTAACTTAATTTTTCTTGGAAAATTATTTTTACA
TAAAAAATCAATTTTACAAGTATAATTATTAATATTATTATATGTTTTAATTTTATATTTCAAAAATAAAACATTTTTTAATTCGAAACTTATAATCGAA
TTTTTTATAAAAGAAAATTTGGTATTTTTAATATATACTTTTATATTTTTTTTATAATACAATACGTTTAAAAAAGAATTTTTAAATACTTTATTTATTA
AGTTATTTTTTAATACATAAATATAATAAATATTGTTATCTATATATAATATTCCTGGAACAGTTTTATTACTTAAATTTTTAGTATTGAAAATTTTTAT
TTTAAATTTAAATTTCATACTTTTTTTTTTTAAAAAATTTTTATTATAATGTAATTATTATTTTTTAAATTTTTTTTTAAAAAGTATATAACTTTTTAAG
TATTTTTAAAAAAAAAAAATAATATATTTCAATTTAAATAAAGTTTTTATAAAAACAATTTTAAAATTTTAAATCAACTTATGAAAAATTATTTAATAAT
AGAAATTAGACAAGGAATAGGAGGAGAAGAATCAATAAATTTTGCGAAAGATATTTACAAAATGTATGTTAAATTTTTTGAAAAACAAAATATTAATTTT
GAAATAATAAGTAATTCTAATTTTAAAGAAATTATATTAAAAGTAGAAAATAATGTTTTTGATAAAAAATTATTAAACGAATCAGGAATACATAGAGTAC
AAAGAATACCAAAAAGTGAAACACAAGGTAGAGTTCATACTTCTACATGTACTGTATTTGTAGCAAATATAAATTTAGATAACAACTTAAAATTAAAAAA
TGAAGATTTAAAAATAGAAATTTGTAAATCTAGTGGATCAGGTGGTCAACATGTAAATAAAACAAATTCTGCTATAAAAATAATTCATTTACCAACTAAA
ATTGCTGTAGAATGCAGTGACGAAAGATCACAAAATTTAAATAAAACAAAAGCTTTAATAATTTTAAATATGAAAATTTTAAAATTTCAAAAAAATAATT
ATAATACACAACTAAATAATGTAAGAAAAAAATTAATTTCAAATTCTGAAAGAGCAAAAAAAATAAGAACATATAATTTTACCAATAATAAAATTACAGA
TCATATAAATAACAAAAATTATTTCCAATTAAACAAAATTTTAAATGGAGAATTTAATTTAATAGATTAATATCTAAAAATTAATTAATTCAAAATATTT
TAAAATTTTAATAATTTTTTTTTTAAATTTATAATATAAAACTTTTTTTTTTTTATTTTTTAAAAAAAAAATAATACTAAAATTATAATAAAAAATATTC
ATAATTTTATTTTTTTTTGAAAAAAAAAAAAAAGAAAAATTTTTTGTATTATATATTTTATTAAAATATATTTTTTTTAAAAAAAAGTATTTTTTTTTTG
GTAAAATATTAAAAATATTTTTTAATCTAAATTTTTTATTGAATTTTTTTATAACTTCATAATTATAATAAGACATTATTCCTATATCATTATCAATATC
AAAATTATAATTTTTTAAATTACTATTTCTAATTTTTATTGTTTTGTTCTGAATATAATTTCTCGTTGAAAAATAAGTTTTTTTTTTAAATTCTAAATTT
AACAATAATCTTTTTACTTCATAAAAAGATGATTTTAAAATGTTTTCATATGAATTTAAATTTTTAATTTCTGTTTTTTTAGATTTTAACAAATTAAATT
TATTTATAATTGATAAATTTACATCAATTCTTAAATTACCATTATACATATTGCAATCTGAAATAAAATAGTTTTTAATAATTTTATTAATTTTTTTTAT
AAATAATATTAAACAATTTAAATTATTAAAATTTGGTTCTGTTACAATTTCAATTAGAGAATTACCTGATCTATTATAATTTATATTTTTCAAAAAAAAA
TTTTTTGTTGATGCTGCATCTTCTTCTAAGTGTATTTTTTTAATTAAAATTTTTTTTTCACAATAATTATTTAACAATAAATTAAATATATTTTGAAAAA
CGTTTTGTGTTATTTGATAATTTTTAGGTAAATCGTAGTAAAAATATAATTTTCTTTCAAAAACAGAAAAATTAAAAAAACAAGAATTTAAGTTAAAACA
AAAAAAAACAATTAAATATTTAATTAAATAATTAAAATAAGGTAATATTCCTGGAATAGCAATATCATAATATGAATTTTTTTTGTTATTTGTAAAAATT
TTATTAGTTTTTAAATGAACGTGAATTTCTATTCCTATTTTTAAAACTAAAAAATTTTTAAATTGAAACATAATTATTTAAAAAAAATTTTTCATATCTT
ATAGCTAATTCTAAAATTAAATTATCTGACTGTAATTCTCCTATTATGTTAAAACCAAATGGTTTATGATTTATCATACCAATTGGTATTGTAATTGATG
GATATCCTATCAAATTTGAAAATACTGTAATATAATCACAATATTCTGAATAATTTATATTATCTATTTTAAAGTTATTTAAAGATGGTATGATTAAAAA
ATCTGCAATAATAAATAATTCATTAAAAAAATTCAATATTTTTTTATTAATTTTTATTTTTTGATTTTTTAAGTAAAAATTTTGTACTCCTAATAAAATT
TTATTTTTACATGTTTTATAGAAACATCTATTTAATTTTGTAAAATCATTTATATTTTTAAATATATTTTTTTTATATCCAAATTTGATTCCATCATACC
TGCAAGAATTAGTATAAAATTCTTTAGATGATAAAATTGTATATTCATAAAATAGTACTGACAAATCAATTTTTTTAAATATAATTGTATAATTTAATAT
TTTAAAATTTAAAATAACTTTTTCAAATTTTTTTTTTGATTCATAATCATAATATAAATTATCAAATAATAATGCAATTATTTTTATTTTATTATATTTA
ACATATAAATATTTTTTAAAAGATAATAAATCTAAATTTTTAGTAGATAAAATATTGTATAAAAACTTACAGTCAGATGAATAATTTGTTATAATAGAAC
AACAATCTAAAGAACTAGAATAAGGAACCATTCCATTTCTTGAAATTTTTCCATATGTTGGTTTAAAACCAATTAAATTAGAAAAAATTGCTGGTGTGCG
TATAGAACCACCTGTATCGCTTCCTATAGAACCAATAACACATCCACTACTAACTGAAATTGCAGATCCAGATGATGACCCGCCAGGTAAATATAAATTA
CTATAAATATTTTTTATGTTAATACAATTGTTTTTTCCACTTTCACCAATACAAAATTCTTCTAATTTATCCAACGATATTACTATTAAATTATATTTTT
TTATTATTTTAATAATTGAAGAATTATAAGTAGAAAAATAATTTTTTAATATTTTAGAATTACAACTTAATAATTTGTTTTTAATGCTATAAATATTTTT
TATAGATATTGGAATACTTAAAATTTTTGGATTTGTATTATCTAGTTTTTTTGCAATTTTTAAAGATTCTTTTTCTAATGTTTCTAGAACAAAATAATTT
TTATCATTTAATTTTGATAAATTATTTATACAAACTTTTATTATTTCAAAGTATGATATTTTTTTATTTTTAATTAAAAAAAAAATTTTTTTTATACCTA
ATTTAAACAAATTATTCATATAAAAATATTTCTTATTTTTACTGTTTTTGTTTTTTTATTTTTTTTTTTTATAATATTTAATTTTAAAAAAAAGCTTTTT
TTAAAATTTTTAAATTTAATTATATAATTTTTTAAAAACTTATTAAATTTAAAATTTAACAATTTACATTTTATTGATAAGTATTTTACGTACAAAAAAA
AAAAAATATTAATCATAAATTTAAAAGAAATTATTTTATTTTTAAAAAAAAAAAATGCTAAATGTTTTTTAAATTTATTTTTTTTCTAAAATTTTGTATT
TTAAAATATGATAATAATTAATTGTTTAAATAAAAAACTAGGAAGAATACTTTCAAGTATTTCTAAAATAATAATTTATTTTAATTTTTTTAAAAAAAAA
ATAAAAATTTTTTTAATTAATATAAAAAGTTTTGTTTTAAAAAAAAAAACAACAATAACTCATTCTGGTTATATAGGTAATTTAAAAGTTAAATGTAAAG
AAAAAAATAAAAAAAAATATATTAAAAAAAGTATATATAATATGTTACCAAAAAATTATAATAGATTAATTTTATTAAAAAATTTATTTTTCATATTATG
AAATATATTTCATTTTCTAAAAAAAAAAAAACAATAACTAAAGTATTAATAAAAATAGGAACAGGAATTATAAAAATTAATAGTTTAAGCATAAAAAAAT
ATTTTGGTAATATTTACAATAAAAATTACTTATTTATTCCTTTAATTTTAGTTAATTTAAAAAACAAAAATATTTTAATAAATACAAAAGGAGGCGGAAA
AAACTCTCAAATTACTTCTATTAAAATTGCAATTTGTAAATGTATATTATTGTTTAATTTTAATTATTACAAATTATTTAGAAAATTAAATTTAATAACA
ATTGATGATAGAAAAATTGAAAGAAAAAAATACAGTTATAAAAAATCTAGAAAAAGAGAACAGTATTCTAAAAGATAAATTATTAATTAATGTTTATTAT
ATGTTAAACATAATAAAGGAAAATCCAATTTTTAAAAAATAAAAAATACTAAATAAGGTTAGTTGCTTAAGATATATTAAAAATACAAAAATTGGTTTAT
GAAAAATAAACATACTACTATATTACCAAATAGACTTTTATATAAAAAAGAATATTTATTAGATTTTACTTTTGGTTGTGGTGGATATATTAAAAATTTA
ATATTAAAACAAAAATTTAAAATTATTTATACTATTGATGTTAGTCAAATATCTTATTTAATTTCTAATAAAATATATAATAAATTTTTTTTTTTTTTTA
GATTAAAAATAAAAAATATAAATAAAATTTTCAAAAGATTTAATTTAATAAATGTTGATTTTATAATATATGATCAAGGAATAAATAGTTATGAAATTAA
AAATTTTTACTATAAATTAAATAAAAAAAAATATTTGGAAAATAAAGTAAAATTAAATGTTTTAAATTTAAAATTTATTTTTTTTAAAATTATAAAATAC
TTTAAAAAAAAATTTAAATTATTAATATTAACGTTTAGTTTGTATGAACATTACAAAATTATATTATTTTTAAAAAAAATAAAAAAAATAAAAATAAAAA
TATTTAAACCAAATAAATTTGAAATCAGTTTAAATAATAGTATTAAAAACGTATTAATTCATTTAATTTATGTTAATTAACATTTTAAAAAAAAAAATAC
TTTTCTTTTATTTTAAAATAATTAGTTTTAAAAAAAAAAAAATTAAAACATTAGGACCATACGGAAATTATTGTTATAATATTTTATTAAAAAAAATAAA
TAAAAAATATTATTTTTATCCTATAATAAGTATAAAAAAAATTTTAAATATTAAAAAATATTTTTTTCCAATTGAAAATAATAACGGTGGTTTAGTAAAC
GATTCTATTAATTTGTTATTTAATAACAATTTTTTTTTTAATTGTATTCTAATCATAAATATTAATCATAAAATATTTTTATATAAAAATAAAAAAAAAA
TATTTTTGCATAATCAATCATTAAAACAAATAAATTATAATTTGATGTTTAAATTTTTTAAATTAAAAATTTTAAAAACATTTTCAAATACAATTATAAA
TTCAGGAATTAACATTTGTAATTCTTTAACTAAGACTATTTTGTTAATAAGTATTAAAAATGTTTTTTTAAAAAATAATTTTATTAATAAAACAAAATTT
ATACTATTTAATAATTTTATAAATAAAAAAGTATTAATTTCATTTTTTATTAATAAAAATTTTTTTTTTTTATTTAAAATTATTAAAAATATAACAAACA
TTTATATAAAAAATAAAATATTTTATATTGAAATTTTTTTTTTCAGTTTAAGAATACTTCTTTTTATTATGAAATTATTTAAAAATAAAATTAAAATTAA
ATTTAAAAGTTTTCATTCTATATTATGAAATCATATTACAAAATAAATTGTTTAAATAACAAAAACAGAATAATTAGAACAGGAGACAAATCAATTTCTC
ATCGTAGTATAATATGTATTATTGTTAAAAAACAAATTATTGAAGTATTAAATTTATTAGAGTCATCTGATATATTATCTACAATTAACTTGTTTAGAAA
TCTGAAACTACAAATTTATGGACCAATAAATAATTATTTATTAATATCTAATTTTAAAAAAAAAAATCAAAAAAATAATATTAATTTTATTGGTAATTCT
GGTACAACTATTAGAATAAGTTTAAGTATTTTGTTTAATAATAATATAATAATTGGAGATAAATCATTAAACAATAGAACTATGTATAGAATAATAAAAC
CATTATCTTTAATTGGTTTTATTATTCAATGTAAAAAAAATTTTTTTACACCTTTAATTATTATAAAAAAAAATAATTTTGGTTTAAAATACAATTTAGT
TAATATTAGTTCTCAAGTTAAATCTTGTTTATTATTATATTCTTTATTTTCTTTTGTTAAAATTTATTTAATTGAAAAAAAAACAACTAGAGATCATACA
GAAAGATTTTTTCACTTAATAAATAAAAAAAAAAACGTTTTAGTTAGAATACCTAATGATTTTTCATCATTAACGTTTTTAATGTGTTATTTTATAAAAC
AAAATAAAAAATTTATTTTAATTTTAAATTTTAATAAATTTAGAATTGGATTTTTTGATTTTTTAGTTATTAATAATATTAATTTTTTTTTTATTTATAA
AAAAATAATTAATAATGAACACATTGTTAAAATATTGTTTTTAAATTTTAAACTTAGAATAAAAACTATTTATTCAAATAATATTAGTAAATTAATTGAT
GAAATTCCTTGCTTATTAATTTTTTTATTAAATTTTAATGCTAAAATAAAAATTTATGGATTAGAAGAATTAAAATTTAAAGAATCAAATAGATTTTTAA
ATATTTATAAAAATCTATTATTATTAGGAATAAGAATTATTAAAAAAAAAAATTATTTAATATTTAAATGTAAAAAATTTCATTTAAATTTTTTTAAAAC
ATTTAATGATCATAGATTATTTATGTCAATATTTATAAATAATAGTTTTAATAAAATTTCAAATGCTGAAAATATTATTTCATCGTTTCCGTTTTTTTTA
AAACTATTTAATAATAAAAAAAATAAATTTTATGTTAAAACTAAATAGTACTATAAAAGAAAATATATTATCATTTTACGATGAATATATTATTACAAAT
TTTTTAAAAACAATTTATTTAAAAAAAATAAATTTTAATTTTAAAAAAATTGTTAATAAGTTTAAAAAAAAAAAAATTGTATATAATATTTTTTTTTTTA
GTAATTTAAAAAAAAAAATAAAAACAGGATATTATTTAGATAATGAATGTCTATTATATTTGTCTAATAATTTAATCGGATTAAAAAAAAAAATTTTACA
AAAAAAAATTTATTTGTTCGAAATTTTATTTAAAATAATTATTTTTAACAAAAGTAAATTTATAACTAGAAAAATTAATTTTAAAAAAAAAAAAAAAATA
AAAATAAAAATAAATAAAAAATACATAGGAATAATAAAAAATGTAGTATCATATGGAATATTCATAGATATAGGAAGCCTAGATGGCTTATTACATATTT
CTGATATTCCAAAATACAAAAAAGTATACAAAAATTTATTTACTAAAAATAAAATTATAGTTAAAATAACAAAATTTGATAGAAAATTAAAAAAAATATC
ATTAAATTTAAAAAAAACATATAAAAAATACAATTTGATATTTGAAGAATATATAATTTGTAAAATTAAAAAAATAGAAAAAAATTTTTTTTTATGTTTA
AGTAAGTACAATAAAATTTTTATAAAAAAAAATTCAAATTTTTATAAAAAAAATGATATAATTAAATTTTACTTAATTAAAAAAAATGAAAATTATTTAT
TTTTAAGCAAGTACTATAAACTATACAAAAACAGAAACAAATATATTAAATATAATTTAAAACTAAAATTTAATGAATTTTACTTATTTTCTTACAAAAA
TAAAAAATTTATTAGCAATAAAAATAATAGAATTACTAAAATTAAAAAAAATGAATTTTATAAATATTTTATATCAAAAATTATAAATAATGATTTTTAT
AAAATCATTTATATAAATAAAAAAAATTATATTAAATATGGTAATTTTTTTCTATTAACAAATATAAATTTTATTTTTACTATTAATAATTTTAAAATTA
TTAATAATAATATATTTTTTTTCTTAAAAAAATACTAAACAGATTTTTCATCCGAATAATCAGCTATCATTGCTTCTGTAGTAATTAATAATCCTCCAAT
TGAACCTGCTGATTGTAGCGCACTTCTTGTAACTTTAACGGGATCAATTATCCCCATTTTAAACATATCACCGTATTTTCCAGTTGAAGCATCATATCCA
AAATTATTTGAAGAAGATTTAATATTGTTTAATACAATTGATGGTTCACCACCAGAATTCTTAACAATTTGTCTTAATGGAGCTTCCAATGCTTTAAGAG
CAATTTGTATTCCATAATTTTGGTCTTCGTTGTCTCCTTGTAAATTTTTTAATTTGTTTAAAATTCTAATGAGCGATACTCCACCACCTATAACAACTCC
TTCTTCTACTGCAGCTCTAGTAGAATGTAGAGCATCTTCTATTCGTGCTTTTTTTTCTTTCATTTCTATTTCTGTTGCAGAACCTACTCTTATTACGGCA
ACACCTCCTGCTAATTTTGCCATTCGTTCTTGTAATTTTTCTTTATCGTAATCTGATAACGATTCTGATATTTGTTTTTTAATAGTATTAATTCTTTTTT
TAATAGAATTCTCATCTCCTCCTCCTCCTAATATTATTGTATTTTCTTTTGTTGAAGTTATTTTTTTTGCAAATCCTAATAAACTTAAATCAATATTTTC
TAATTTAATTCCTAATTCTTCACTTATTAAAGTTGAACCAGTTAAAATAGAAATATCTTTTAATATTTCTTTTCTTCTATCTCCAAATCCAGGAGCTTTA
ACTGCTAAAATTTTTAAAACACCTCTAATATTGTTTATAACTAAAGTAGCAAGTGCTTCTCCTTCAATATCCTCAGCTATTATAAAAAGCGATTTGTTTT
TTTTTGAAATTAATTCTAAAATATTAACAATATCTCTAACATTAGATATTTTTTTATCAGTAATTAAAATTAAACAATTTTCAAGTATTGAAGACATATT
TTCTTGGTTAGAAATAAAATATGGAGAAATATATCCTCTATCAAATTGCATTCCTTCTACTACTTCTAATTCATCTTCAAATCCTCTACCTTCATCAACA
GTAATAACTCCATTTTTTCCAACTCTATTCATAGCGTCTGAAATAATTTTTCCAATAATTGTTTCTCCATTTGCTGAAATAGTACCAACTTGAGAAATTG
ATAACGTGTCAACACATGGTATTGAAATTTTTTTTAATTCTAAAACAGCTTGATATATTGTTTTATCAATACCTCTTTTTAAATCCATTGGATTTATTCC
AGAAATAACAGCTTTTATTCCTTCATTTACTATAGTTTGAGCTAAAACAGTAGCTGTTGTTGTACCATCACCTGCAACATCTGATGTTTTTGATGCTACT
TCTTTTACCATTTGTGCTCCCATATTTTCAAATTTATCTTTTAATTCTATTTCTTTAGCCACTGAAACTCCATCTTTTGTTACCAGTGGTGAATTAAAAG
ATTTATCTAAAATAACATTTCTACCTTTTGGTCCTAAAGTAGTTTTAACAGCATCTGCAAGAATATTAACTCCTATTGCTAGACTTTTTCTTGCATCGTC
GCCAAATTTTATTTTCTTATAACCCATTTTCTATAATAGCGATTACATCTTCTTCTTTTAGAAAATAATAATTAATGTTATCACTTTTATACTTTTCTAT
ATTATAATTGTCTTTAAATAATATTATATCTTTTTTTTTAACAATTAATTTTTTAATTTCACCATTTTGTAATAACTTTCCACATCCAATTTCTATTATT
TCACCTTTTATTAAATTATTATCATTAAATGGTAAAAAAATACTTCCAATTTTATTTTCTAATTCAATTTTTTTAACTATAATTTTGTCATACAATGGTA
AAAATTTCATTTTTTAAGTTTTTATTAATAAAAAAAATATTTTTAAATAACACAATGTATTTATTTTTAATTTAAAAATATGTTTTAATTTTAATAACGT
TATTTTAATAAATTATTTGTTTGAAAAAAAATAATATTAATTAATAACATCAATTAGGTAATTAACTCAATAGGTAGAGTATCAGTTTTACATACTGAAA
GTTATAAGTTCAAATCTTATATTACCTATAATAATTGGAGCGATAATTTAGTAGGTTAAAATGTTGGCTTGTCACGTCAAATATCGCGGGTTCGATTCCC
GTTCGTTCCGTAATTTTATTAAAAATTATTATTTTAATATTTAACTTTAGAAAGTTTCTTAAAAATATTAAAAACATAGTTAAATCAAATTTTGAATTTT
TAATAATAAACAATTTTTATAAAAATAATATTTTATACTTTAACATTAATAATTTTTATTTAATATATAATTTTTATAAAAATTTAAATTTTAAAAAAAA
AAATTTTTTAAATTTGTTTTTAGTCAAAAATTTTGGAAATATAAATTCTTGTATTAGAACTTGTTATATGTTTAATGTTTATGTTATAATTAAAAACATA
TATATTAATAGCTTTGTAATTAATTATAATAATATAATTTTTATTAAAAATAATATTTTATTGTTAAAATTCATAAAAAAAACAAATTGTATTGTTTCTT
TATCTATAAGTAGTTATTTTGTTTTAAATAATTTTAAATTATTAAAAAAATTTGTTATTGTTATTGGTAATGAAAAATATGGAATTAATAATTCTATTAT
TTTGCATAGCGATTTTATTCTTAAAATTAAATCTTATAGAAAAAAAAGTTTAAATTTGTCAATTGTAAGTGGTATAACATTGTACCATTTTATATATAAT
AATGAAATATAATATTTTTATTTTTTTTAAAAAATATAAAAATAAATTTTTTTCTATTTTTTTTAAAAATTTAGTTTTTTTCTTTTATAAATTTAATATT
AAAATTTTAAAAATTATAGATTTTGGAAATATTATATCATTTAAAAAAAAAAATAAAAGATTATTTTTAATTGAAATAGAATGTATTAAGAATAAAGTAA
TTACAATTTTTAAATTATTTAAATCAAAAAACGATGTTTTATCTTTCTTTATAATATTAAAATTAAGTATAATTAAATTAATAACAAATAACAATTATAA
AAAATATTTGTCTAAAAATTTTTTAATAATACCGTCTTTACTATTTAAAATAAAATTTATTAAATATAAGAAGATAGTAAAATTAATTAAAATTTTAAGA
ATATTAGGAATTCTTCCTTATTCAATAAGAAATTATATTTCAAAAATAAATATATTATGCTAGAAAATATATCATCGATTAGTTCAGAATATATTATTAT
AAATTATATTTTTAACAATTTTAAATACTCAAACAATATTTTTAATATAATAAACGAAAATGATTTTTTCTATGAAAAAACTAAAAAATATTTTATAGAA
AAAAAATTTTTTATAAAATTAAATTCGATAGATGAAACATATATTTTTTCAAATTTAAAAAATTTAGTTGAAAAAACAAAAAAAAGAAAAATGATAAAAA
TTTTATTTAATATAGTTTTTACTCTTTTAAATGAAAAAGATATACCTGTAATAATTTATGAAAAAGTTAAAATGCTATTAAATTTTAAAGAAACAAATAA
AATTATGAAGTTAAATTATTTTGAAATGTTAAAAGATTATTTGTTCAATAAAGAAACGCCGATTTATACAGGATATAAATCTTTAGACAATATTTTAAAT
GGATTACAAAAAGGTGATCTTATTATTCTTGCTGGTAGACCTTCAATTGGAAAAACTTCTTTTTTATTAAATTTGGTAAAAAATTTACTATTTTTGAATA
ATAAAATTATAATATTTTCATTAGAAATGACTGTTTTACAAATTTTCATAAGATTAATTTCTATAATATCAGAAGTTAATCAAAACAAATTTAAAAATAA
TGATTTTTCTGAATTTGATATTAAAAAATTATCTTGTTTATTTAAAAATTTTAATTTTAATAATTTAATTATTGTAGATTGTTCTTCATTATCTCCAAAT
GATATTGAAATTCAACTTAATTTTTATAAAAAAAATAATTTAATAATTAATATCATTTGTATTGATTATATTCAATTAATGAAATCGGAAATATCGAATA
ATAATAGAGTTTTAGAAATTTCTGATATTTCTAGATCATTAAAATTAATAGCAAAAAATTTTAATTGTGTGATAATATCGTTATCGCAATTAAATCGTTT
AATCGAACATAGAATTGAAAAAACACCAATTTTATCAGATTTGAGAGATTCTGGATCTCTAGAACAAGATGCTGATATTGTTATATTTTTAAATAATAAA
TATAAATTAGTTGAATTATCAATTCTTAAAAATAGAAACGGTCCGCTAGGTAATATTCTTTTCAGTTTTATAAATGAATATACAAAATTTAATCAAGTTT
AAAAAAGAATTTGGACAAAATTACGTTTTATTTCAAAAAAAAAAAAAAATAGATTTTTGTTCTGGATTTAATTTTTTTAATTGTTTAACTAATGATTTAG
ATTATTTAAAAAGTTGTTTTTGTAAAAATAAAATTTATATTTTTAATAACAAAAATTTTTATAAAAATTTTTTAATTTTTTTTTTTCCAATAAATTATAT
TGTTAACATAAATTTACCATTTAATATTATAAATTTTTTTATTAAAAAACAATTTAAAATTAAAAATTTTTATTTTAAATTAATAATTTATTTTTTTTTT
ATTACAAAATTTAAAATAACAGGTATTGTTAACAAAATTAAATTTTTTAAAATTAATATTTTTTATCCTTTTACAAAAGTGAAAATAGTTAATATTTTAT
ACAAAAAAAATGTATATTTAAAGTTTATTTATAAAAAAATAAAATTAAACAATTTCAAATTTATCAAATCTGTTAAAATGCAATTTTTTATATTTAAAAT
AATAAATTCATTTAATAAAAATGTCTTTTGCAATAAGATATAAAAAAGAAAGATTTATTTCTTTTTTAAAAAGATTTAAAAAATCAATAAATAGAAAAGT
TATATTTATAAATAAAAAGAAAAATGTTAAAAACAAAATATTAAAAAAAATAAAATGAATTTAAAATGTCCATTTCATAATGATACAAATGCATCATTAT
CAATTAAAAATAATTTCTATATATGTTATGGTTGTAAAAAAAGAGGAAACACAACTTTTAAATTTTTAAGTATTAATAATTTTGATAATTTATTTAATAA
AAATACTTTATTTATTGCAAGAAATAATTTATTATTTAAAAAAAATTGTTGGTTAAATTATTTAATCAAAAGAAATATTAGTTTTGAAACATTGTTAAAA
TATAATTTAGGTTATGCAAATATTGATTTTAAATACAATAATAAAAAAAAAATGTTATATAATAGATTAATTTTTCCAATACTAAATGAAACAGGTATTT
TAATAGGTATAGGTTTAAAATCTAAACAAAATAAACCAAAATACATTAATTTATTAAAATTTAATTTTAATAAAAATGAATTAATTTATGGAATTTATGA
AAAAAAAAACAATAAGTTTGTAATAATTGTAGAAGGATATTTTGATTTACTAACTTTGTACGAAAATAATATTTATAATGCAATTTCGTTACTTGGTTCA
AATATAAATGAATATAAATTATTATTTATTTTAAAAAAATTTAAAGAAGTATTTTTTTGTTTTGATGGTGATCATTCTGGATATTTAGGAGTATTAAAAT
TATCTTTTTTTAAAAAAAAACACAAAATAAAAAATTTATTATTTAAGTGTTTACCATCAAAATATGATCCCGATTTATACATTAATAAATTTGGTATAAA
AGAATTTTTAAATTATTTAAAAATATGATTTTAAAAATTAAACTTAAAAATAAAAATTATAAAAATATTAAAAATATTTTTAATAAATTTTGTTTTAAAA
AAAAAAAAAAAAAATTTTTTTTAAAAATAATAAATACAAAAAATTTATTCAATAAAATAACTATTTTAAGAAAAATACTTTTAAAATACATTTTTAGAAT
AAAAATTTTTAGATTTTTTTTTAAAAAAATCAATTGTTTTTTCTTTTTTTTAAATTTTATTTATAAAAATACGATTTTTTTTAAAAAAAGATTTAAAATA
AATTTTAAAATACTTAAAATAATAAAAAAAAATTTTTTAAAAAAAAAAATTTTTTTCAATTGTTTTAGTATTTTTAATAATAAAATTAATTATTTTTTAA
ACAAAAATTATTTTATTAAAATGTTAATATTTATTATTTACAAAATAAAAATTTCAATAATTAAAAAAAATTTAATTGCTTTTAAAAAAAAAATAAAAAA
ATATAAAAATTTAATAACAAAAGAAGAACAACGAGAAAATTACTTAATTTTTAAAAATTTAATTTTTAATTATAATTTTTGTGAAAATAAATTTCATATT
AATTTTAATTGGAAAATTAAAAACAAAATTATAAAAAATAAACTAAAAAAAAATATAATAGTTTTAAAAAATAATTGTTATGAAAATATTAAAATATACA
ACGAACCCAAAAAAACATTATTTATTAGAGAATATGATAATAGAAATCTAAATAATTTAATAAGATTAATAATTTTAGAATTACCTAAAAAAGAACAAAT
AATAATAAGATTAAGATTTGGAATAGGATTTCCTAAAAGTTATACTTTAGAAGAAATAGGGTTAATGTATTACTTAACAAAAGAAAGAATAAGGCAAATT
GAATTAAACGTACTTTTTAAATTAAGACATCCTACTAGATCAGAAGTTTTAAAACCATATATAAAATTATTAAATACTGAAGAATAAAATTTTAAATTGG
GGTTATAATTTAATTGGTAAAATATTTGGTTTGCAACCAAAAAAAAAGAGTTCAATTCTCTTTAACTCCATTTTTTTTATTTTTTTTTTAAAAATAAAAA
AAAACGTTTAAATTTATTTTAAATATTTTAAAACCAATTTTTTTTTTTTTAAAAGAATATATTTTTTTTAAGAGCTTTAAAAAAATTAAGTTATGGCTAA
ATAGCTCAGTAGGTAGAGCAAAGGACTGAAAATCCTTGTGTCGGTGGTTCGATTCCACCTTTAGCTATGAGATATAGCCAAGTGGTAAGGTATTGGTTTT
TGATACCAAATATCCTAGGTTCGATTCCTAGTATCTCAGAAAGCTTTTGTAGCTTAGTAGGTAAAGCTGTTGATTTGTAATCAACTGTCTCGGGTTCGAT
TCCTGACAAAAGTATTTTTTTGGAGATGGTGAGAGTTGAACTCACATTCATTATAATTTTATATATTTATATTACATGATTTATTAGTTAATATTTTGTA
TTTTTTAAATTATAATTATTTAATATTGTATATTCTTGTTGAATTTACTAAATTATTTAATTTCATGTCTTTTTTTTAAAAATTAAAATTTTGTATTTAT
AATGTCAAAACCCGTCATCCCCTTTTTAACAATAGATATTTTTAACTTGAATAAAAAATTTAATTTAAAAAGCTCTATTGGAATACATTTATAATTTTTA
TTATAATTATAATATAAAATTATGTTTTGTTCCCTTTTATTGAGTAATAATATTATTTTTTTTTTTTTATTTATTAAATTAATATTATTATCAGATAATT
TTATTATATAATCATTGATATTTAAAATATTTTTTTTTATAATTAAAATATCATTTGATTTTAAAATCATTCCTGCTATAAATGTTTTTATGATAAAGTA
TTTTGTTAATATTTTTTTATTTTTTTTTAACATAATTTATGGAAAAAGAATTTTCATTAATAATTTTTTTAGAAAGATTAAAACATTTAGAAGAAAAATT
TAATAAATTTTTAACAGATACAAAAAATATATGTATAAAAACAATATCAGAATTAAATGTTATCAAAGAAAATTTAATTAATGATATATTATTAGAATTA
ATACCATTAAATGATTCAATGGAAATGTTTTCCAAATCTTTTAAAATAAATCAAACTGGAGAAATGGAAATTTTAGTTTTAATTTTTAAATTAATAAATA
AGTTTTTTTACAAATTTGAAGTTAAACAAATATCAAAAATTGGAATTAGTTTCAATCCAGAAATACATGAAGCAATAGGAATGTATCCTACAAATTTAAT
TAATAAAAAAAACACAATAAAACATGTTTTACAAACTGGTTATAAAAGAAAAATTAAATTATTAAGACCTGCATTAGTAATTGTTTACAATTAAACAATG
AAAATAATTTTTTTTTTCTTATATAAAAATATAATAATATAAAAAAATATGAGTAAAATAATAGGTATAGATTTAGGTACTACTAATTCTTGTATAGCAG
TATTAAGTAATGGTAAACCTCAAGTTATTGAAAATTCTGAAGGTGGCAGAACAACACCTTCAGTAGTTGGTTATACTGAAGATAATAGAATTATCGTAGG
TTTACCTGCTAAAAGACAAGCAATAACAAATCCTAAAAATACTTTATATGCAATTAAAAGATTAATAGGAAGAAAATTTAAAGATGATATTGTACAAAAA
GATATTAAAATGGTTCCTTATAAAATAATAAGTTCAGAAAACGGAGATGCTTGGGTTGAAGTAAAAGATAAAAAATTAGCGCCTCCGCAAATTAGTGCTG
AAATTTTAAAAAAAATGAAAATAACTGCAGAAAATTTTTTAAATGAAAAAGTAACTAAAGCTGTAATTACTGTACCTGCTTATTTTAATGATTCGCAAAG
ACAAGCAACAAAAGACGCTGGTAAAATTGCTGGATTGGAAGTATTAAGAATAATAAATGAACCTACTGCTGCAGCGCTTGCTTATGGTTTAGATAAAAAA
AAAAATGATAGAATCATTGCTGTTTATGATTTAGGTGGCGGTACATTTGATATCTCAATAATTGAAATAGCAAATGTTGATGGAGAAACACAATTTGAAG
TTTTATCAACAAATGGAGATACTTTTTTAGGTGGAGAAGATTTCGATATTAGAATTATTAATAATTTAATATATGAATTTAAAATAGAAAATGGAATAAA
TTTAAGTGGTGATTCTTTAGCAATGCAAAGATTAAAAGAAGCAGCAGAAAAAGCTAAAATTGAATTATCTAGTGTAGAACAAACTGATATTAATTTACCA
TATATTACTGCTGATAAAAATGGACCAAAACACTTAAATATTAAAATTACAAGATCAAAGCTTGAATCATTAGTTGAAGATTTAATTTTAAAATCTTTAA
AACCATGTGAAATTGCGTTAAATGATGCTAAAATAAGTAAAAATAAAATTGATGAAATTATCTTAGTAGGAGGGCAAACAAGAATGCCTTTAGTTCAAAA
AATGGTTAGTGATTTTTTTGAAAAAGTTGTTAAAAAAGATATTAATCCTGATGAAGCAGTTGCAATTGGTGCTTCAGTTCAAGCTGGAGTATTAAGCGGA
GTAGTAAAAGATGTTCTTTTATTGGATGTAACACCATTAACATTAGGAATTGAAACAATGGGAGGAATAATGACTCCATTAATTGAAAAAAATACTACAA
TTCCAACAAAAAAAACTCAAGTTTTTTCAACTGCAGAAGATAATCAAACTTCTGTTACTATACATACTTTACAAGGAGAAAGAAAAAAAGCTTTACAAAA
TAAATCTTTGGGTAAGTTTGATTTAAACAATATTTCTCCTGCACCAAGAGGTGTACCTCAAATTGAAGTATCTTTCGATTTAGACGCAAACGGTATATTA
AATGTAACAGCAAAAGATAAAAAAACAGGAGTAGAACAATCTATTGTGATTAAATCATCAGGTGGTCTTTCCGAATTAGAAATAGAAAATATGATAAAAG
ATGCTGAAGCAAATTTAGAAATTGATAAAAAATTTGAAGAGTTAGTAAAATGTAGAAATGAAGCAGATAGTACTATCTCAATAGTTAAAAAAAAATTAAA
AGATGAAAATTTAAAAATTTTAGATGAAGAACGTGTTTCAATAGAAAAATCTATTTCAAACTTAGAATTATTAATAAAAGGAGATGATATAGATTCGATT
AAAAAAGAAAACGAAGAATTATTAAAATTAAGTGATAATATTATAAAAAAAAAATAACAAAGCATTACTCTTTATTTTATTAAAGAGTAGTGCCTTTTTT
TTAAAAAAAAAAAATGATATTTAATATTAAAAAATTATTAATTAAAAAAAAAAATTTTATTTTAAAAGTTAGAAAAATATGTTTGTTTTGTTATAATTAT
TTTTTTTTAAAATGTAAAAATTGTTTTTTCTTAAATCTTATTAAAATTAAGTATTTTTTTTGTGAAAAATGTGAAATAAACTCTAATTTGTTATGTAATA
AAAATAAAGTTTTTTATATTTTAATATTTGTTAAAAAAAAAATAAAACGTAGAGATGAGTTATTTTTAAATTATAAAGAAATAAATTCTAAAATTTTAAA
CAAAGTTATAATTAAATATAATGATTAAAATTTTTATATTTGGAATAACAGGAAAAATAGGAAAAACAATTTTAAATTTTATTAAATTAAATAAAAATTT
TATTTTATTAGGAGGTATTAATAAAAAAAATTATAAAAAATTTATAAATAATAAATATAACTTTATTTTTAAATTAATGACAAAAAATAGTGTTATTATT
GATTTTTCTAATCATTATATGATTAAGAAAATTCTATTTGTTTCTTTATATTATAAAATTTCATTAATTATAGGAACAACAGGTTTTAATTTTAAAGAAT
TAAAATGTATAAAATATTGTTCTAAGTATATAGCTTTAATACTATCATACAATATGAGTATTGGTATCAATATATTAAACTTATTTTTTTTAAATTTGAA
TTATTTTTTTTTAAAATTTAATTTTAATTCTATAATAATAGACATTCATCATAATAAGAAAAAAGATAAACCATCAGGAACAGCATTGATTTTATATTCA
AAACTTAAAAATATTAATTTAAATATATTTAGTTCAAGAATTAAAAATATAATTGGGAATCATATTATTTATTTAATTTCAAACTTTGAAATTTTAAAAT
TTGAACATTATGTTATTAATAGAAATATATTTATAATAGGAATTTTTTATTCAATTATATGGTTAATGAATAAAAAAATGGGATTTTTTTCAATGTATAA
TGTTTTTTTTTCATGTTAATATTAGAAAATGGTTATTTAATAAATTGTTTAAGAATAAATAAAAAAAATATTTTTGGTGAATTATCATTTAGTATATCTA
ATTACGGATATATAGAATCTATTTCTGATCCTTCTTATAAAGGTCAAATATTGATTTTAACAAATTCTTATATAGGAAACGTTGGTTATATAAATCAAGA
TATACAATCAAATAAAATTTATATAAATACAATAATTTCAAACAATTATTCTATATCTAGCAATTTTAGATCTAATTTTAAATTATTTGATTTTTGTAAA
AAAAAAAAAATTCAAATATTAACAAATTTAAATACTAGATTATTGATTTATTTAATAAGAAACACAGGTTCACAAATAGGATGTACAATAATTACAAAAA
ATAAAAAAAATATAATTTTATATATTAGATCATTAACGCTAAAAAAAATTTTATATGATTAATTTATTAATTATAAATTTAGGAACAAAATTAAGTTTGC
TTAAAAAATTAGTTAGTAATAATTATTTTATAATAGAATTAAAAGATGAATTTTATTTAAAAAATATAGATGGTTTATTTATTTCTAATGGTCCTGGTTA
TCCTAAAAATTTTTTAAAATATAAAAACATTATTTTATATTTTTTATATTATAATATACCTATTTTAAGCGTATGTTTAGGACATCAAATTATAAGTATT
TTAAATAAGTTTAAAATTTTTAAATTAAAAATAGGTCATCATAGTTGCAATCATACTTTATATAATGATTTTGAAAATAAAATATATATTACTTTGCAAA
ATCATAATTTTAATATTAAAAAAAAAAATCAAAAAAATTTAATAAATAATTTTAAGTCTTTATTTGATAAAACATTTCAAAATATTTCATCTTTAATATT
TCCTATATTGAGTTTTCAAAATCATCCTGAAGGATGTTCAGGTCCAAATGATTTAATTTTAGTTTTTAAATTTTACAAAATTAATAATGAATAAAAAAAT
ATTAGTAATAGGAGCAGGTCCTATTTTAGTTGGTCAAGCTTGTGAATTTGACTATTCAGGTAGTCAATCATGTAAAACTTTAAAAGAAGAAAATTTTAAT
GTAATTTTATTAAATTCTAATCCTGCAACAATAATGACTGATTATGAAATTGCAAATACTGTTTATATTGAAAAAATAAATAAAAATAGTTTAATAAAAA
TCGTTCAAATTGAAAAACCAGATTTTATCTTACCGACAATGGGAGGCCAAACAGCTTTAAATTGTATTTTAGATTTTATTAATTCAGAGTATAATTTTCC
AGAAAATAAAATTTTAGGAATAAACAAAAAAATTTTAATGAATGCAGAAAGTAGACTTATATTTTATAAATTGATTAATAATTTAAATTTAAAATGTCCA
GACTCTGTTATAATTAAAATTTCAAATTTAGAAGAAATTAATAATGTAAATTTTCCATGTATTATAAGACCTTCTTTTACTTTAGGAGGTTTAGGAAGTG
GAATTGCTTATAATAATAAAAGTTTAAATTTAATATTAAAAAATGCATTTTTATTTTCAAATGAAGTTTCATTAGATAAATCAATTATAGGTTGGAAAGA
ATTTGAATTAGAATTATTAATAGATAATTATAATAATATTATAGTTATATGTTGTATTGAAAATATTGATCCAGTTGGAATTCATACAGGTGATTCAATC
ACTATTACTCCAGCACAAACTATATCAGACAAAGAATATCAGAATATGAGAGATTCTTCTTTTATTATACTAAAATCAATAGGATTAAAAGGAGGTGGAG
CAAATATTCAATTTGCTATTAACCCAATTAATGGTGATTTAATCGTTATTGAAATGAATCCTAGAATTTCAAGATCTGCTGCTTTATCTTCTAAAGCTAC
TGGTTATCCAATTGCAAAAATTTCAACTAAATTATCGATTGGATATAGTTTATTAAAAATATTTAATAATTCTAAGTATGGAAAATTTATTTCAGGATAT
GAACCATGTATAGATTATATCGCTATAAAAATACCAAAATTTAATTTTGAAAAATATTTAGAATTAAACTTTTTAAATACTGTTATGAAATCAATAGGAG
AAGTATTAGGTATTGGTTTTTCTTTTCAAGAAGCTTTTTTAAAAGCAATATATTCAGTATTTGAAAATAATAAGATTCCATCTTTTTTAAAAAATAAGTT
TTATAATTTATATAAAAATAAAGCTATAAAAAAAATTATAAATTCTAATAGCACAAAAATATTTAATATAATAGATTTATTTAGATTAAATGTAAATATT
AAATTTATTTTTTGTGTATCAAAAATTGATCCATGGTTTTTATTTAATATTAAAAAAATTATTGAAGAAGAAAAATTTTTTTTTAATAAAATTAAAAATA
ATAATTTAAATTACAAAAAAATAGATTTATGTTCTAATGAATTTGAAAAACCTTCTTTGTATTATTATTCTAGCAAAAATTTATATTTTAATGAATTAAG
ATATTCTATAAAAAAAAAAATTATAATAATTGGAAGCGGAACTAATAGAATAGGACAAAGCATTGAGTTTGATTATTGTTGTGCAAAATTTTCCAAATTA
ATTAAAAAAATAGGTTTAATTTCTACAATGATAAATTGTAATCCAGAAACAGTTTCTACTGATTATGATACATCTAATTACTTATTTTTTGAACCGATAA
CTATTTTTTTTATAAATAATATAATTAGTTTTATTAAACCAATATTAATTATTTGTCAAATAGGTGGACAATCTCCTATAAATTATTTAATTAAAAATAA
TAAAATTAATAAAATTTTTTTTGGTTTTGATAAAAATAATATTTTTTTTGGAAAAAAAAAAAAATTTAATAAAGTATTATACATTTTAAAATTAAATAGA
ATTAAAAATTATAATTGTATAAATATAAACGATTTAATACTTTATTTTAATATTTTTAACAAAACAATAATAACTAGATTACCAAATATTATTGGTGGTG
CATTAATGAAAATTTTGTGCAATAATAAAGAATTTTTAGATTTTATTAATTTTAATAATATAGTTAATATTTATATAGAAAATTTCTTACTTGATTTTAA
AGAATTAGATTTAGATGTATTAGTAGAAAATGGAAAAATAATTGTTTTAAGTATAGTAGAACATATAGAAAGTACTGGAATTCATTCTGGTGATTCTAGT
ATGATTTTTCCGTCTTATTCGTTATCTTCTAAAATAATTAAAAAAATTTTTTTTATAATTAATATTTTTTGTTATAAACTTAAATTGAATGGAATTATAA
ATTTTCAAATTTGTTTTAAAAAAAAAATTTATATTATTGAATGTAATCCTAGAGCGTCTAGGACTGTTCCTTTTATTTCTAAATCAAATAAATATACACT
AATATATAATTATGTTTTAATTTTGTTAGGATATAATATATGTTTTATTAAAAAAAAAAGCAATTTTTACTTTATTAAATCTTCAATATTTCCGCTTAAT
AAGTTTAAAATGTTTAAATTATCTCCTGAAATGAAATCTACAGGTGAAGAAATGAATTCTGGTTTAAGTATACAAGAATGTTTTTCTAAATTATTTTTGC
TATCAGACAAAACAGTATTTTTAAATACAAAAAATACATATATATATATAAAATATTTAAAAAAAATAAATTTTAAATTAAAAACTAATAGCTTTAAAAA
AAAAATAAATTTAGATTTGCACATAAATTCATATGTTTTAAAAAGAAAAAACAATATGTATATAATTTCCGATAAATTAAATATTTTAAGTAATGATATG
TTTTCTAATAAAAATACAATAATTTTTTTTTTATCTTCTTTTGATAAAAAATTTTTATTAGTTAGAAAATTAAATAAAATTTAAACCGAAATGGTGAAAT
GGTAAACACTCTATTTTGAGGTAGTAGATTTTACGGGTTCAAATCCCGTTTTCGGTATTTTTAAATTTTTATATAATGGATAATATTAAAAACATATATA
TAGATAATAAAGTTTTAATTAAAGATTTATCATTTTTAATTAATAAACCTATTGAAAAATTAATAAAAGACTTGTTCTTAAATGGTATTTTTGTTAAATT
AAATGATTATTTATTTTTTGAAAATGTTAAAAAAATATGTAAAAAAATTTATAACATTGATGTTTTTAAAAAAGAAAACATTAAAATATTAGAAAATAAT
AACAAATTAAATAAAAGAACAGATGTGATCTTTATTACAGTAACTGGAAATGTAAACAATGGAAAATCTTCTCTGATAGATTTTATTTTAAAGAAAAATA
ATGTAAAATTTGAAGTTGGAGAAATTACGCAAAATATTTCTGTATTTAATTTTTTTTTTTTTGAAAAAAAAATTTATTTATTTGATTTACCTGGTCACTC
ATTATTTAGTAAATTAATAAATATTAATTTATCAATTTCAGATATAATATTTTATATAATATCGTACGAAGATAATATTGATTATAAAAAAATTAACGAT
ACAATTATAAAATTTGAAAAGTTATCAATTTCCATAATTTTATGTATTAATAAATATGATAAATTTAAGTTTGGTAAAAAATATTTAAATTTTAAAAATG
AAAAAATTTATATTTCAGCAAAAACAGGATTTAATATTAAAAAATTAATAAATACATCAATTTTAATTTTTAATAAAAAAAACAAATATATAGATTTAAA
TAATCCAGGAAAAGGAATAATAGTTAATAGTTGTTTTAAAAATGATGTATTAATTACAACATTATTTATTTTTAAAGGAACACTAGAAACGGGCAATTTT
TTAAATTTTAAACATATTAGTATAAAAATTTTAGAATTTTTTGTTAATGAAAAAATAACAAGTAAAATTGAATCACCTAATATTATATTAATTAAAAATA
TTCAATTTCCTATTGAAGTATTGTTTGAAATAAGTAATACAAGAAAAAATTATTTTATTGGAATAGATTATAAAGATAAATATTATAATTTTACTAATTT
TTATATTAAGGTTAGTACTCATAATATGGGTTTTTCAATTATGAATTTTTATAATGATTTAAAACTAAACGAGTCTATAAATATTGTTAAATTAAGCATA
GGCGTTTTAAATGATAATGATATTAATTATTGTTTAAATTTTAATTGTATTATTATTACAGTTGGAATTCTTATAAGTAGTGTTTTAAAACAAAAAATTT
TAATTAACAAAATCAAATTTAAAGAATTTGATTTGGTTAACGATTTAATAGATTATTTTAAAAATTTCTACAAATGTAATAAAATTGAAAAAACAATTGG
AAAACTAAAAATACAAGAAATTTTTCCATCTGGAAAATTAAACAAAATAGCAGGTTGCAAAGTCATTTTTGGTGAAGTTAATTTAAAAAATAATATTAAA
ATTTATAAAAATTTAAAATTAATATTTAAAGGCAAAATAAAATCAATTAAAATTAAAAATCAAACTAAAGAAATTGTTATTTTAAACGAAGAATGTGGAA
TTCTTATTAAAAATTTTAATAATTATGAAGTAGGTGATATTATCGAATCATATATTTATGAGTATGATAAAAAATATATTGAATAATAAAGGATTTGGCG
ATCCTAAAATTCAAAATTTTTTTTTAATTAAAAGATTAAAAAAAATTAAAAATCATTTTTTAATAAATAAAAAAGATTTAAAATGTAAAATAGTAATTTC
AAAATTATTATGTAAGATAAAAAAAAATATTAATTATATGAAAAATAAATTATGATATATTTAAACAATATATTTGAAAAAATATTTATTAAAAATAATA
ATTCCAAAATAATAAAATTTGGAATAGATCCAACGTTTTTTTCTATTCATTTAGGACATTTATTTATAATTAATTATTTATTTTTTTTAATTTATAAAAA
ATTTATAATTATAATTATTATAGGTGACTATACTACTAAGTTTAAAAAAAAAATAAATTTAAAAAATTTGATAATAAATTCAATTTGTTTAAAATCACAA
ATTAAAAATATTTTAGGAGAAATAGATGTTGTATTTAATTCTATATGGTATAATAAATTTAATCTATGCTATTTTATAAATTTAATAAATTTAGTTAGTA
TAAAAAATTATATAAATAAAACTTTAAAAAATAATTTAAGCAAAAAAATAAGTAATTATATTTATCCTACAATACAATCTTATGATTCTGTTTTCCTAAA
ATCAGGTTTTGAAATAGGAGGATTAGATCAATTATTAAATATTATCTGTGGAAGAATATTTCAAAGTAAATTTAATATTAAAAAACAAAATATTATTACT
TTAAAAATATTAAGTATAAATAATATTAAAATTTCGAAATCAAAAAATAAACATTTATTAAATGCTTATAAACAAATTTATAATTTTAAACTTTTAAAAT
TAGTTTTTTTTAATTTTAAAAATTACAATAAATGTGTTTATAAAAAAACATTCTTTTTAAATTTAGTTTTACTTAAAAAATGCAAAATTTTTTTTTTTAA
AAAAAAAATTTATTTTTTAAATAGAAATTTTATTAGTTTTTATTTCAAGAAAATATTTAATTTAAAAAAATTTAATTTTTATAAATTAATTTATAATAAA
AATATAATTATTAATAAAAAAATTTTAATCAAAAATGTTTATTTTAGAAAACACATTTTTAATATTAAAAATTATAAAATTTTAGTTTATGATAAATTTT
ATAAGTTTGTACTTAAAAAAAAAAAATAATTTTAGATTTCCACCTGATCCAAATGGAAATTTACATTTTGGTCATACTTTTAGTATTTTTATAAATAAAA
ATTTATCAAAAATTAAAAAAGGCAATTTTTTTTTAAGATTTGATAATACAAATTTAATAAATAATTTTAGTTTTTTTTATAAAAATATTAAAAATGATAT
TTTATGGTTAAACTTAAAATGGAATGGAAAAATTTTATTTTTTCAAAATAAAATTAATATTTTTTATAAATATTTAATTATTTTTTTTAAAAAAAAAAAA
TGTTATTACAAAAAAAAAAAAATAATTAACAGATTTTTTTTAAATTATATTAAAAAATTAAATGTTTTTGAATGTTTTATATTTAAAAACAATTTTTATG
AAAAATATAACTTTGTTATTTTAATTAAAAAAAAAATTATTTACAGAAAAGTAAAAAAACAAAGACATTGGATTATAAATTCTACTTATGATTTTTCTCA
ACCTATAAACGATTATTTAAATTTTATTTCAATATCAATTTGTACAAACGAATTTAAAAATAATTCTAAATTTTATCATTATATTTTTAAAAAAAAAATA
TTACCAATTCAAATTGAATTTAAAAAAAAAAACTTTAAAAACACTAAAATTTCAAAAAGAAAATTAAAATTTAATAAAATTTATAATTTTTTTTTTTTAC
GAAAAATTGGAATTACACCAAAGATTTTAAAAATTTATACTAATATAATAGGTATTTCTAACAAAAACATTTATTTTAAAAAAAAAGATTTAAAAAATTC
TATTTTTTTTGAATTAAACTATTTATTTAAAAACTGTTGTTATTTTAATAATTTTATTAAAGTGAAAAATATTAAAAAAAATATCGTAATTTCTAGTTTA
TTTAATTTTAAAAATATTAATTTTTATTTTTTTAATATTAGATTTAACTTTTTTTTTTTTATAAAAAAAACATATTATAAAAAACTTAAAAAGATTTTTT
TAAATAAAAAAAGATTTTTTATAAATAAAAAAAAAAAGATATTAATAGTTTTAATAAAACTGATTAAAAACAAATATATTAATAATAATATTGATAAAAA
TCATTATTTTAATAATAAAAAAAAACTATTATTTAAAAAACAACTAATTAGTAATGTTAACTAATTTTTTTCGTTTTTTTTTTAAAATTTTAAAAAAAAA
TATCGGAATATAGCGTAGTTTGGTAACGTACTTGCTTTGGGAGTAAGTGATCAAAGGTTCAAATCCTTTTATTCCGAGCGTTTATAGCTCAAATGGATAG
AGCAGTGACCTTCTAAGTCAAAGGTTGTAAGTTCAATTCTTACTAAGCGTATGGTAAATATAGCTCAGCTGGCAGAGCAATAGTTTGTGATACTATTGGT
CGCGGGTTCAAATCCCGTTGTTTACCTTTTTATTATTTTAAAAATATTTAAAATGATTTTTTTACCAACTGTAATAGATAAAAATATAAAAGGAGACAGA
TTTTTAGATCTATATTCTAAAATGTTAAAAGAACGAGTAATTTATCTCAATGGTGCTATTGAAGATACGATGGCATCACTAATTGTATCACAATTAATGT
TTTTAGATTCTGAAAATTCTAAAGATATAATATTATATATTAATTCTCCTGGTGGAGTTGTATCTTCTGGTTTAAGTATTTATGATACTATACAATTTTT
GAAATCAGATGTTTCAACAATTTGTATTGGTCAAGCAGCCAGTATGGCTGCAGTTTTGTTAGCTGCAGGTAAAAAAGGAAAAAGATTTTGCTTTCCTAAT
TCAAGAATAATGATTCATCAACCATTAGGGTATGCTCAAGGGCAAGCAAGTGATGTGGAAATTCATGCGCGTGAAATGATAAATATTAAACAAATATTAT
GTGAAATATTATCTTCTCACACAAATAATAGTATTTTTCAAATTTTTAAAGATACCGATAGAGACAATTTTATGAATTGTAAACAAACTTTAAAATACGG
AATTATAGATAATGTTTTGTATAAAAAATGGATATAAAAGTATTAATTTTTCAAATTATTTAAATCCTGAAAAAATAAAATTAGAATTAGATAGATATAT
TATTGGACAAAATGAAACAAAAAAAATTATTTCTGTTGCAGTTTATAATCATTATAAAAGATTGTTTCTTATTAAAAGTAAAAAAATTTTACTTGAAAAA
AGTAATATAATTTTAGTTGGACCAACTGGATGTGGTAAAACATTAATGGTTAAAACTTTAGCAAAAATAGTAAATGTTCCAATTATTTGTGTAGATGCTA
CTTCTTTTACTGAAGCTGGATATGTAGGTGATGATGTTGAATCAATTATTCAAAAATTATTACATGAATGCAATTATAATGTTGAATTAACAGAAAAAAG
TATTATTTATATTGATGAAATTGATAAAATTTCTAAAAAAACTGATTTTTTTTCTGGAAAAGATGTATCAGGAGAAGGAGTACAACAATCAATGCTTAAA
TTAATAGAAGGTATTACTTTAAGTATTCCTTCATTAGTTGAAAAAAAAAATTCTCAACAAATATTTAATATAGATACAACAAATATTTTATTTATTGTTG
GCGGTGCTTTTTCTGGTATTGAAAGTATAATTAATTTTAGAATAAATCAGGAATTAAATTTTATAAAAAAAAATTTTGAATTAACAGATATAATAAATTA
TACTAGTTCGGAAGATTTAATTAATTTTGGTATTATTCCAGAATTTTTAGGTAGATTGCCTATAATAGCTAAATTTAAAGAATTAAGTGAATCAGAATAC
ATTTATATTCTTATTAAACCTAGAAATTCATTAATAAAACAATTTTGTTATTTATTTTTAGTTGAAGGAGTAAACATAAAATTTACTTTTAATGCAATTA
AAGAAATTGCAAGAATTGCAGTAAAACGTAAAATTGGTGCAAGAGGTTTGAAATCAATTTTAGAATTTGTTTTATTAAAAGCAATGTTTATTTTTCCTTC
TAAAAATAATTTAAAATTAATTTTAATTTACAAAGATGTAATTGTTTTAAATAAAACACCATTGTTTATTTATAAATAAAAGTAAATAAATGTTTATTTT
ATTTTCTCATTTTTTTTTTTTTTTAAAAAAAAACAAACAAAATTTTTAAATAATTATTCAATAAAAAAAAATACAATAAAAAAAAATAGTTTTTTTTTTT
ATTTAAATAATTATTTTTTTTTTAAAGACTTATTTGCTTCAATAATAACGTCAAAAATATATGTTTTTGTACAATTTTTTATTATTAAAAATGACATAAT
TGGTAAAAAATTTTTAAATTTATTAAAAAATAAATCTTATGAATGTAAAATAATTATAATTGTTGATAGAATTGGAACTTTATTTTTAAGTAAAAAATTA
AATATTTTATATTTTAATACAAATAAATTTTTATTAAATTATAGAAATCACAAAAAACTTGTTTTAATTGATAAAAATTTATTATGGTTAACAAGTTCTA
ATATAGGAAAAGAATATTTTAATTTAGATTTAAATGTTTATGATTGGAAAGATTTTTATTTCAGAATTAAAAATTTTTATAGCTATTGTTTCTTATTAGA
TTATTATTATGTTAAAAATATAATACATAAATCGTTTAACAAAGTTTTTGTAACAAGAAATTTTTTTTTAATTAATAATTTAACATTTAATTATGTTGTT
TTATTATTATTTTTAATTAAAAAATTTTTTATAATAATTTCTCCGTACATTGTTATTGATAATTTTTTTATTAAACTTATTAAAGTTTTAATATTAAAAA
ATATAAAATTAATAATTGTACTATCAAAAATTGCAGAAAATTTGTATACTCATGTTTCTTCAATTATTTTTTTAAAAATTTTAAAATTAATTGGTATTTC
TTTTTTTTTTACAAAAAATGGTTTTAATCACAGAAAAATATATATAATTGATAAAAGTTTAATTTTTTTTGGTTCAATGAATTTTGATAATAGATCTATT
TATTTAAATTTTGAATCTTTATTTTTAATTACCAATAAAAATTTTATTAAAATTTTTTTAAAATCTTTGTTTATTAAAATTAATTGTAATTTTTATAATT
ATAAAAAAAAAAAGATAGTTTATAAAATACTATACATTGTTTCATTTTTAAATTACTTAAACATATGAATAATAATATTTTTAAAAAAGAAAATTTCTTT
ATTTTATACGGTATAATTAAAAAAATTAAAAAACTAGGTAAAATAATATTTATAGAATTTATTTCTTTTTCAAAAAAAATAAATTTTTTAATTAAAAATA
TAAAAATTAAAATTTTTAATACAATAACTGGAATTTATTTTTTTAATAATAATAATTTAAATATATTTGAATTATGCTATTTTAAAAAAAAAATTAATAA
AATTAATTTAAAATTTTTAAAACTAAAATCTAAAATAATATATTTTATTAGATTGTTTTTTAGTATAAATAATTACCTAGAATTAGATATACCAATTATT
GAAAAATATACTAGTTCAGGTTCAAAACAATTTTTAATTATTGATAAAAATAAAAAAAAATATTTTTTAAGTTTAACACAATCACCGCAAAAAATTAAAC
AATATTATATGTTTAATGCGATAAATAAATATTTTCAAATTGCAAAATGTTTTCGTGATGAAGATTCTAGATCTTCTAGAATAAAAGAGTTTCAACAAAT
AGATATTGAAAATTCAAATACAATATTTTTAAATTTTAAAAAAAAAATTAATTTATTTTTAAAATCACTTATTTTTTTTATTTTAAAAAAAAAAACATTA
ATATTAAAAATAAAATATAAATTCATTAAAAAATATTTGTTTGAAAAAAAAAATTTAAATTTGCCTTATTTGTATAAAAAAAAAATAATTAAAAATTCAT
ACATTTATATACTTAAAACTAAATTAAAAAAAATTGAAATAAATAAAAGTTTTTATTTTAAATTGAGTAAATATTATATAATATTGACGTTAAAAAAACA
AGATTATAATTTTTGTTTAAGTTTAAAATTATCGCATAAATATAATAATTTAATTAATTTAAATATTATTTTATTGTGGATAATTGATTTTTATTATTTT
AAAAATAAAAAAATTAAACATCATCAATTTACCGCATTTAAGAATAATTTTAAAAATTTTTATAATTCTAAATCCTTAGCATACGATGTATTTTTAAATG
GAATAGAAATAGGTGGTGGATCAATTAGAAATATTAATTTTTTAATACAAAATAAAATTTTTTTAAATTCTAAAAAAAAATCAAAATTTATAAATTTTTA
TAAAAGAGCTTTACCACATCATTGTGGTATTGCTTTCGGATTAGAAAGAATAATAAGTTTATTAATAAAAAAAAATATTAAAAAAACAATTACATATTAT
AATTATTCAAAATTAATAAAATCAAAAAAAATAAATGAATAATTTTCAAAACATTGGAAAGATGGTAGAGTGGTTTAATACATCGGTCTTGAAAACCGAT
AAAGTTTATTCTTTCCAGGGTTCGAATCCCTGTCTTTCCGATAAAATTTATGTTTCATTATTATCTATAAATAAAATTTATATAAAAAAAAAATTGATAC
AAACATTATTATCTAAAATTAAATGTTATCATTTTGATTTAATGGAATTTTCATATGTAAAAAATAATTCATTTTCAATTAATGAAATTAATTCAATTTT
GTTAATATTAAGTAAAATAATGAATATAAAATATGAAGTACATGTTATGAGTAAATATTTACTATTAAGTAAAATAGACAAAAATAAAAGTATTAATCAT
TTAGAAAACAAAATGTATAAAATTAATAATATTGCTTTATCAACTAATTTTTGTTGGAATTATATAAAATACTTTAACTATAATAATATATTAATAATGT
CAGTTATTCCTGGTTTTGGAAATCAAAAATTTTTAATAAGTACTTTAAATAAAGTTAAAAAAAAAATAAACATAGATGGTGGTGTTAATTGTTATATTTT
TAAGAATATAAAAAATTATTTTAATAAAATTATTATTGGTAGTAATATTATTAACATAAAAAATAAATTAAGTTTTTATAAGATAAATTTTATCTTAAAT
GAATTTAAAATATAATTTATTCATTTAGCAAAATTACATTTCCTCCTTGAGCAGTTGTATTTATAGTAATTGTTTTTTCATTAACTAATGAAAGAAGATA
ATTACATCCTCCTGCTTTTGGACCAGTTCCAGATAAATTACATCCACCAAATGGTTGCATTCCAACAATTGCACCTATTATATTTCTATTAATGTAAATG
TTTCCTATTTTTAAATTATTAGATAAATATTTGCAAAATGCTTCATTTCTACTATGAATACCTAGTGTTAATCCAAATTTAGAATTATTAATATCAGAAA
TTATTTGATCAATTTGAAAATTTTTAAATATTGATACATGTAAAATTGGACCAAATTGTTCATTTTTTAAATCATTAATATTATTAATTCTTATTAATGT
TGGATTTAGAAAATTACCTAATATTGGTTTTTCTTCTATAGAATAAACACTCTTTTTGTTATAAGTATTAATATATTTATTTAAATTAAAAAAATTTTTT
TTTGTTATTATAGGACCTATATCATAAAATAAATTTAAGGGATTACCAACATTAATATTTAATAAAACATTCATTAATAATTTTAATGTATCAAAATAAA
TGTTTTCATTTATATATATTAATCTTAAAGCTGAACATCTTTGTCCGCAACTTTTGAAAGCGGATTCAACTACATCATATACTACTTGTTCAATTAAAGC
AGTTGAATCTGCTATAAGAGTATTAATTCCTCCTGTTTCTGCTACTAATTTATATAGAGGAGCATTTTTTCGCATTATTAAATTTTTGCTAATAGTATTA
GCAACTTCATTAGAACCTGTAAAAATAATTCCACAAATTTCTTTATGAAACGAAATTTCATTTCCAATAGAAAAACCTGAACCAATTATTAATTGACATA
CACTAATTGGGATTCCTGCTTTAAATAATAATTTTATTAATTTACAAGCAATTAAAGATGTACTTTCTGCTGGTTTAACTAAAACTACATTTCCAGATAA
TAATGCTGAAATTAATTGTCCACAAAAAATTGCTACTGGAAAATTCCAAGGACTAATTGCAGAAAATATTCCTTTTCCTTCTAGCATATATATGTTATTC
TCTCCTGTTGTACAAGGTAAAAAAACTTTTTTATTTAACAAAATTGATTGATTACAATAATAATTGCAAAAATCTATTGCTTCTTTAATATCTGAGATAC
AATCTATAATTGTTTTACCTGCTTCTATAGAACATAATAATATTAATTCTATAAAATTATGTTTTATTAATAAAACAAAATTTTTTATTATTTTATTTTT
TTTAAAAATACTAATATTTTTCCAAAAATTAAATGATTTTTTTAATATTTTAATTGAATTATTTATACTTAAATCATGTTTTAATAATCCAATAACTGAA
TAATTATTAAAAGGAGGATAAACAATAGTTTTTCTTTTAATATTACTTATAAATGAACTTGCTATCCATTTTTTATTAGTAAATAAATTTATTTTTTTAT
AAAAATTATAATTATGAAATGAAATATTCAAATTATAAAAATATATTGAACTTCGAATTCCTCCAAATAAATCTGTTGGTAATGGTATTTTATTATTATA
TTTTTTATTATTAATTTTAAAAGGATTTTTGCTTAAAATTTCTAAATTAACATTTTTATCAATAATTTTATTAACAAAAGAAGAATTTGCACCATTTTCT
AATAATCTTCTTACTAAATATGGAAGTAATTCTTTATATTTTCCAATTGGAGCATATTCTCTATATGTAATATTATATATTTTTTTTAACGAATGATAAA
CATCATTTCCCATTCCATATAATTTTTGAAACTCATAATTTTTATCATTAGATAATGATAAAATAAAAGAAATTGTTTGTATATTATGTGTTGCAAATTG
TGAAAAAATATTTTTTTTACATAAATTTGATAACATATACATTGAACATAACAAATAAGATAAGTCTGTACAAAATTTATTTATATAAACAGGATACATA
GGTAAATTTAAAGTTTGAGAATATTTAATTTCATAATCCCAGTATGCTCCTTTAACTAATCTTACAGGTATTATTTTTTTTTGTTCTATAGAAATATAAT
TTAACCAGTAAAGAATAGGAATAGCTCTTTTAGAATAAGCTTGAACAACTATTCCAAAACCTTCCCAATTTTTACATATTTTAGAATAAAAAATATTATT
AAAAAGAATTAGAGATAATTCTAATCGATCGACTTCTTCTGCGTCAATAGTTATTGATACAAATGCTTCTTTAGCATTATAAATTAATATTTTAATTAAT
GGTATCATGTCTCTAGTAATTTGTTCAATATTATAAAAAGAATATTTAGGATTTAATGCTGATAATTTTATTGAAATTGATGGTAATCTTTCATTTAAAT
TATGTTCAGTATAACATTTTTTTATTTCGTTAATTGCTAGTTTATATTGAAAAAAAAATTTTCTAGCATCATAATATGTAAGTGCAGCTTCTCCTAACAT
ATCAAAAGAATATTTATTTTTATCATTTATTGATTTTTGTATTGCATTAATTATATTACTGGAATATACAAATTTTTTACCAATATGTTTCATAACATAA
TTTGAAAAAAAAACAACTGTTCTTTTAAAAATTTTTATAAATATTTCTTTATTTTTAGTGTAATGTTTATAATATGCAACAGAACAAAAATCAATAATCA
TATTATAAAATATAATTTTCCAATAATCAGATTCATAATAGTAAGACCAATCTTGAAACGAGATTTTATCTTTTATAAAAGAATCTGCAGAATAAAAATC
TGGAATTCTAAGCAAAGATTCTGCTAAACACATTAATTGTATTCCTTCTTTTGTACTTAAATTATATTCTCCTAATAAATTATCTAAATTATCTAAATGT
ATATTTTTTCTAGTTTCATTTACTAAATTTAAAGATATTTTTTTTACTTTGTTAAAAAAATCATTTGTAAAATTACAATTTTTTAATAATTCTAAAAGAT
AGATATTTTCTTCTATTAAATAATATTTACTTATTATATTTAATAATTGTAAATTATTCATTTTAATTAAATTAATATGTTTAGTTGGACTAGTAATTAT
TTGTATAAAAAACCAATATGGGTAAGTGAAGATTTAAGAGATGGAAATCAATCGTTGATAAATGGTTATAATTTAAAAAATAAAATTAATATTTGGGAAT
TTTTAATTGAAGTAGGATTTAAGCAAATTGTTTTAGGTTTTCCATCTTCTAATAAACATGACTTTAACTTTATAAATTATTTAAAAAAAAATAGACTAAT
ACCAAACAATGTTTTTGTTTCTGTTTTAACTCCTGCGAAAACTAATTCAATAAATTTAACAATAGATTCATTAAAAGGAATTGAAAATTCTATAATTCAT
TTATATAATTCAATATCAAAGATTCAAAGAAAATTAGTTTTTAAAATGAATAAAAATGAAATAAAAAATTTTACAATAAATTTTTTTTTATATACAATAA
GTAAAATAAAAAAAAAAAATATTATATTTCAGTATTCTCCAGAAAGTTTTTCTGATTGTGAATTAATATATTCAAAAAAAATTTGCTATATATTTTCATA
TTTGTGTTATATTAACAACGTAAAATCTATTATTAATTTACCTATTACTGTTGAAAATATTTTAAGTAATAAATTTGTAAATTCAGTTTTATATATTAAA
AAAAAAAAGTTTAATTCTTTATTATCAGTCCATACTCATAATGATATGGGAGGAGCTATTACTTCTTCGATTTTATCTTTATTATCTGGAATTGATCGTA
TTGAAGGAACTTTACTTGGAAATGGTGAAAGATCTGGTAATTCAGCTATTATGATTTTAGCTTCAAATTATTATAATTTAGGTATTGATCCTGGTATAAA
TATATTTAATAATAAAATTTTTTATTTTTTAAAAAAAAATAACAAAAGTAGAATTCCATGGTATTCAAATTTAAATTACGTTGCTTTTTCAGGAAGTCAT
CAAGATGCAATTAATAAATCTTATTTTAAAAAAAAAAAATTCAATTGGAATATTATTTATGTACCAATAAATCCAAAAATTTTTAATTTTAAACATAAAA
ATATGATTAAAATTAATATTCAATCTGGAAAAGGAGGATTAAAATTTGTTTTTAATTATAATTATAAAATAAAATTAAATAAATTAATTTTAATAAAATT
ATATTTTATAATTCAAGATATATCAGAATATCTAATGACTGAAATTTATAAAGAAATGATATTTTCAATTTTAATAATAAGATCTAATTTATTATTTATT
AAAGATTTTAAAATTATTTTTTTAGATATTTCTTTCTTATATAACTTTAAAATAATAATTATCGTTTTAAAGAAAAATAAAAAAAAAACTATAAAAATTA
TAAATTATTATAATGAATAATTCATACGGTGAAATTATTAAAATTTCAACTTTTGGAGAAAGTCATGGTTTAATTATTGGTGCTTTAATTGATGGTTTTT
TTTCAAATTTATATATTAGTGAAAAATTTATTCAAAAAAATTTAAACTTAAGAAAACCATTTACTTCATTATTTTCAACACAAAGAAGAGAACAAGACAA
AGTTAAAATTTTCACCGGAATTTTTAAAAATAAAACAACAGGCGCACCTGTATTAATGTTAATAAAAAATAATGATAAACAAAGTTCAGATTATAATAAT
ATAAGTTTAAATTTTAGACCTGGACATGCAGACTATACTTATTTTTTAAAGTATAAATTTAGAGATTATAGAGGTGGAGGTAGATCTAGTGCTAGAGAAA
CAGCTTGCAGAGTTGCAAGTGGATGTGTGTTTAAAAATTTGATTTATAATAAAGGAGTTATTGTTCGTTCATATATTAAAAAAATTGGTTTTTTAAAAAT
AAATTTTAAATATTGGAATTATACATTAAATAGATTTTTTTCAAATTTATTATTTATAAATGAGATTAAAGATATAATTAATAATTGTAAAAATTCATGC
AATTCGTTAAGTTCAGAAATTGTAATTATTATCAACGGTCTTGAACCAAGTTTGGGAGATCCTCTTTATAAAAAAATTAATTCTACTATTTCTAATTATT
TGTTAAGTATTAATGCAACTAAAAGTATTTGCTTTGGTTTTAACTTTAAAAATAAAAACTCATTTCAAGTAAAAGATGAAATTAAAAATTCTGGATTTAC
TTCAAACAATAATGGAGGAATATTAGCTGGAATAACTAATGGACAACCTTTAGTAATCAAAATATTATTTAAACCTACATCTAGTACTTCTAGAAAAATA
AAAACAATAAACGAAAAATTAAAAAATATTACAAATAAAACTTATGGAAGACATGATCCTTGTGTTGGTTTAAGAGCTGTACCAGTAATTGAATCTATGT
TATATACAATATTAATAAATAAAATTTTAAAAAAAAAAATTTATGAATAAAACCATCTATGATAAAATTTTTGAAAGTCATGTAATTAAAAAATATAATA
ACTTATATATTTTATATATAGATAAAATTTTATTACACGAAGTTACGTCTCCGCAAGCTTTTATGTCAATAGGAAAAAAAATATTGTGGAACAAAAGTAG
TATTTTTTCTACTTCAGACCACAATGTTTCAACTAATTTTAAACATAGATTTTTTTATAATAAAAATTTAAAACAGTTAAAATGCTTAAAAAAAAATTTT
AAAAAATTTTTTTTTAAATATTATGATATTAATAGTCCAAAACAAGGAATAATTCATATTATAGCATCTGAATCTAAAATTTTATTACCTGGAATGATAG
CTATATGTGGAGATTCTCATACAACAACAAATGGAGCATTGTCACTAATAGCAAATGGTATTGGCACAACAGACATAGAAATTGGAATATCAACTCAATG
TATTATTCAAAAAAAATTAAAAAATATGAAAATAGTAATTAATAACTTTTTGAATAAAAATGTAACTTCGAAAGATTTAATTTTATTTATAATAAAAAAA
ATAACTTCAAAAGGAGGAACAGGATATTCAATTGAATTTAAAGGAGATTGTATTAAATCTTTATCAATTTCTGAAAAAATGACTTTATGCAATATGTCTA
TTGAAGCTGGTTCAAAAATAAGTATAATTTCTCCTGATGTTAAAACTATTAATTTTTACAAAAAAAAAATAAAAAATATAAAAAAATTTATTAATTATTT
AAAACAAATTAAATCAAATAAAAAATCTTTTTATGATAAAACATTTTATTATAATGCAAAAAATATTTATCCTCATATTACTTGGGGTAGTAATCTTGAT
ACTATCATAGAACTAGATGAATTAGTTCATTCTGACAATTTTAAAATGTTAAAATATATGAATTTAAAAAGTAACAATTCTTTATATAAAATAAAAATCG
ATAAAGTTTTTATCGGATCATGTACAAATTCTAGATTTGAAGATTTATTAGTTTGTTCAAAATTATTATTAAAATTAAACAAAAAAAAACACAAAAATGT
TATTGCGTATGTAGTATCTGGATCAGAAAATATAAGATTAAAATGTGAATTTTATGGTATAGATAAAATTTTTAAAAAATATGATTTTATATGGAAAAAT
TCTGGATGTTCAATGTGTTTAGCAATGAACGAAGATAAATTAAAACCAGGAGAAAGATGTGTATCAACATCTAATAGAAATTTTGTAGGAAGACAAGGAT
ATAAAAGTATTACTCATTTGTCTAGTCCTATATTTGCTGTTATATCTGCTATTTATGGAGAATTTATAAATTTTAAACTTTATAATTTAATAACAAATGA
TTTTGATTTCTAAATTTTTAATACTAAATATTAATAATATTGATACTGATTTAATTATTCCAAAACAATTTTTGAAAACAATAAAAAAAACAGGTTTTTA
TTATTGTTTATTTTATGATTTAAGATATTTAATAAATCAAAATAATATTTTTTTAAATTATGATTTTCCATTTAATATTAAAAAAAATAAAAATGCTAAA
ATATTAATTTCTAGAAAAAATTTTGGATGTGGAAGTTCTAGAGAACATGCTGTTTGGGCAATAAAAGATTTTGGAATTAAAATAATTATTGCTGAAAGTT
TTAGTGATATTTTTTATGATAATTCTTTTAAAAATAATTTATTTTTAATAAAATTAAAAAATTTTGAAATTAATTTTATAATTAACAATTATGAAATAAA
TATAATTTATATTAATATAAAAAATCAATTTTTAAAATTTAATAATAAATTATTTTATTTTAATATAAATAATTTATATAAAAACATATTATTAAGCAAT
TTTTCAATAATAGACTTTTTATTAGAAAAAAAAGATATAATATTTTTGTTTTATAAAAGATGATTAATTTATTAATTTTACCAGGAGATGGTATAGGACC
AGAAATAATTAAACAAGTAATTAAAATAGTTAAATCATGTATTTATACTGGTTATAAAATAAATATAATTTACAATTATATTGGTGGTATTTCAATCGAT
AAATTTAATACTCCTATTACAAATAATTTAATAAGTATTATAAAATATATAGATACAATTTTTTTAGGATGTGTAGGAGGATATAAATGGAATCATTCTA
TATTTAAACCAGAATATGGTTTATTAAAGTTAAGAAAAAAATTTAATTTTTTTACTAATATAAGACCAATAAAGTGTCCTTTTAAAAATATAGATATAAT
AATTGTAAGAGAATTAAACGGTGGTATTTATTATGGTAAACCTAAAGGTTTTTCTAAACAAATTATAAATCAAATACCAACATGGTATGCTTATAATACT
AAAATATATAATGAACAAGAAATAATAAGATTAGCTAGAATTAGTTTTAATTTAGCTTTAAATAGAAAAAAAAAACTGTGTTCAATTGATAAATCAAATG
TATTAGAAACATTTAAACTATGGAAAAAAACAATAAATTATGTTCATAAATTTTATAATAAAGTAAAACTTTCTCATATTTATATTGATTATGCAACAAT
TGATTTAATAAAAAATTTTAACAAATTTGATGTAATTATAACTTCAAATTTATTTGGAGATATAATATCTGATTTATGCTCGTTATTAACGGGATCATTA
GGTATGTTACCATCAATTTCAATAAATAATAAATCACTAAGTTTGTTTGAACCATGTCATGGTAGTGCTCCGGATATTGCTAACAAAAATATTGCAAATC
CTGTAGGCGCTTTACTATCATTAGTAATGATGTTTGAATATGTTTTAAATGATTTTAAATTATCTAATAATTTATATTATTCTATTTACAAAGTTTTGTC
ATATGGTTTTTGCACAATTGATATGAAAAAATATATTAGAAATTTTAAAATAGTAAGTACCGAAGAATTTGGTGATTTAGTAAATCATTTTTTTATAATA
AATGTTTAAATTAGGTATTATTGGATGGAGAGGATTAGTTGGATCAGTTTTTATAAATAGAATATTTACTTCGAATATAATCAAGTATTTAGAAATATAT
TTATTTTCTACTAACAAAATTTTAAATTTTAATTTAAATAATGCATTTAATTTAAACAATTTAATAAATATGAAATTTATTGTTTGTTGTCAAGGTAGTA
ATTTTACAAAAAAAGTTTTAAAATTATTATTATTAAAAAAATGGAGTGGTTATTGGATAGACGCATCTAGTTATTTAAGAATGAATAAATTTTGTACGTT
AATATTTGATCCAATTAATAAGATTAATATTTTGAAAAATATTAAAAATCAAAAAATATATTCAGGAAGTAATTGTACTGTTAGTTTGTGCTTATTAACT
TTTAGTAATCTATTAAAATTAAACTTAATAGATTGGATTATTGCTACAAGTTATCAAGCTATTTCTGGTGCGGGAAGTAAACTGATTAACGAATTAGTTA
ATAATATTAATAAATCACATAATTTATCAAAAAATTTATTAACATTAGAAAAACAAACAAAACAATCATTTAAAAAAGAAAATCCTATATTATTTAATTT
AATACCTTGGATAGATAAAAAAGTAAAATTTAGTCAAACAAAAGAAGAATGGAAATCCTCATCTGAAGCAAGTAAAATTTTAAATAGAAAAATTTTAATT
GATTCAAATTGTGTAAGAGTTTCATCTTTAAGATGTCACTCTCAACTTTTTACTTTTAAAGTTAATAAAAATATAAGTATTAATGATTTATATTATATAA
TTAATAACAAATTCATTAAAATAATTAAAAATAATGAAATTGATAGTACTAAAAAATTAAATCCTTTTAATGTAAGCGGTAATTTAAGTTTATTTGTAGG
TAGAATAAAAAAAAGTTTAATAGATAATAGAATTTTTAGTTTATTTAGTATAGGAGATCAACTTCTTTGGGGAGCTGCTGAACCTTTAAAAAGATTTTTA
GAAATTTTAATAGAAGAATTATTATAATTTTTTAATCATGAAGTACAAAATTATACTTATTAATTCAATAATAAAAAAAGAAAATTTTTCAATTTTTAAT
TTTAATAAAAAAAATAAAAATTTTTTTTTAAAAAAAAAATTAAATTTTTTTTTTAAAATTAATTTTACTAATTATAATTATATTTATAACACTTTTAATA
ATTATAATAAAAATATTATTTATTCTAATAAGAAAAAAAATTTTTTTTTAAAAATCCCAAATTTGCAAAGAAATTGTTTAATTTCTAAAATCACTAGTTC
AATAAAGTATAATTTATCTAATTTTTTTTTTTGTTGTTTAAATAAAATTATTGATTTTTATTTAACAAAAATTTTAAGTTTAAAATTTAGATATAATGAA
ATTAATATTCCAGTTTTTATTAATTATTCTAATTTATTGTTTTCAGGTCAACTACCAAAATTTTATAATTTTTTATTTAAAATTGAAAATAAAAAATGTT
TTTTAATACCTACATCAGAAGTAATATTAAATTCTCTTTCGTTTTTTTTAAAAAAAAAAATAAATCAAATTAAAATATTTTGCAATAGTTTATGTTTTAG
AAAAGAATCATATAATTTACAAAATAGTTCTGGTTTTAAAAAACAAAATCAATTTAAAAAAATTGAAATTTATCAATTTATAAATAAAAATATTTCATTA
ATTGTATTTTATAATATGTGTAGTACTATATTTTATATTTTAAAATCTTTAAATATAAAATTTAAAATTATAAAAATTAACAATTTTGAACTTAATCCTA
ATACATTTTATTCATTTGATTTTGAAATTTACATTAATAATTGGTTAGAGATATCATCTTTATCATTATGTCTTGATAAACCTTTTTTTTTTTATTTAAA
AAAAAAAAATATGCATATAATTAATGGATCATGTTTTCCTATAGGAAGATTAGTATTAGCAATATTGCATTATTATCGATTAAATAATAGAATTTTTAAA
GTTCCTAAAAAACTGAACAAATATTTAACCGAACTGTTGAAATGGTAAACAATCAAGATTTAGAATCTTGTGCTAACGCTTAGGAGTTCAAATCTCCTGT
TCGGTATATTTTTTTTTATAAAAAAAAAAATATATATATATGATAAGTACTATCGAAAAAATTAAAAATAAGATTATTACAATTAATATTTATAATTTTT
TAATTAAATTAGGAACACAAAAAAAATCTTTTTTAAAAAATAAATACTTATTAAAAAATTGTTTAATAAATACTTGGATTAAAATTATTATTAAAAAAAA
AATTTTTTTATTTGGATATTCTAATTCTTTGATATTATCTGGATTAATAAAAATTATTAGTAAGGTTATTAACAATAACATAAAATTAAACGTTAATATC
TTTTTGAAATACAATTTATTAAAAATTATTAAAATTAAAAATATAATTACAAATGTAAAACAAAATAATTTTAATAGTATTGTTAATCATATTAAATTTA
AAATAAAGTAAAAATGTTAATGTTTAAATTTATTTATTAACTTTTAAAAATTTTAAAATTAATAAATTTTTATATTAATATTTTTTTGATAAATTTTATT
TTTTTAAACATCTTTAAAATAAAATTATTAAAAGTGTTTTTCGATATAGAACAAGTAGAACAATTACCTATAAATGATATTAATAACGTTTTTTTTTTAA
AATTAAATATATCAATTTTAATAGCTCCATTATGTAAAATTAATTTTTTATTAATAATACATTTTATTAAAAAAAAAATTTTGTTTAATATTTTATTTTC
TAAAAACATGTTATTACTTAAAAATGGAGCATTTATTAAAATTATTTTTTTATAATGAATAAATAAATTTTTTAAAAAACAGATACTGTTTATATCAATA
ATCAAATTAATTTTATTAATTTTTATTAATTTGCTAAATAAATTTTTTTTATAATAACTATAGTATAAAAATACTTTTGAATTTTTTTTTCCTTTATTTT
TAATAAAAATTTTATAATAAAAATTATTTTTAATATTGTTAAATAAATAATAGTAACATGAATCAGATATATAAAATCTATAATTAAATAAATTTATAAT
ATTCATATTTCAATACAAAAATCTTTTTCTATAAAACCTGTTTTTATAATTAAAACAGCTAAAAGAAAAGCATCATTTAAAGCTCCATGAATTCTTTTTT
TGTTTTTAATCAAATTAAACTTTAAACATAAATCATTTAAATTATTTTTTTTTCTTGGATATAATTTTCTAAATAACAATAATGAATCTAATATATTAGC
ATAGTTTTGTATTTTTTTAATTTTAAAATTAGTTAAATAAATTTCTTTATTAATAAAATTAATATCAAATTTTGCATTATGAGCAATAATAGTAGAATTA
TTAATAAAACCTATAAATTCATTTATTTTTTCGTAAAACATTGGTTTTAATAATAAAAAATCATCTTTTATACCATGAACATTAAACGCTCCTTGTGTTA
TTTTTACTTCTGGATTAAAATAAGAATGAAATACTCTGCCAGTTAAGTGTCCATTTATTACTTCTACACATCCAATTTCTATTATTCTATCACCGTGTTC
TACAAATAATCCTGTTGTTTCTACATCTAAAAATATTATTCTTTTCATAGATTTATATTTATGTTAAATATGTAATACATTTTTTTTAAAAAAAAAATTT
TAATTTCGTTAATTTTAACAATTAAGTTTATTATTCCTGATTTCAAATATTTGTTATATTTTTTTTCTTTTATTATTATTAGTTTATTATTATTAAACAA
ATAATTTAAAAAATTTAAAAAAACAATTTTGATTAAAATTGAAAATCTAATTTTTTTTTTAAATTTACTTTTTTTTTTAATACTAAATTCAAATATTTTA
GTTTTACAATATTTTAAAAATTTAAAAAATCCAACGCAATTTATATTTTTATCTATTACAAAAAAAAATTTTTTTTTATTATAATTAATAAATTCGCAAT
TAATATTATTAAATATAAAAAAAAAAAATATATTTTTTAAGTTTTTAGAAAAACAAATGTTTTTCGAAATTAAATTTATATAAATTTTGTTTTTATAAAT
ATTATAGTTAATATTATAACAATTTTTTTTATTTTTTTTTATAGAATTTAAAAAATTTAGCTTATTTACAAAAATGTATTTTTTTAAAAAAAAATTTTTA
ATTAAAAAGTTGTTTTTTAAAAAAAATAAATTAGTATAATTGATTTTTTTATAATTTTTAAAGAACAAGTTTAAATTTAAATATCTTTTTAAGTATTTTT
TTTTTTTTTTAAGAACAAAAAATTTTAAAAAAAATTTTATTTTTTTTTTACAAAAATATATCCAATATGGAGATTCATATTGAATCCAATTTTTTTTTTT
TTTTATTATTATAAAACCATAATTTATAAAAATTTTTATTTTCGTAAAAAAAAATCTTTTTATTATTGATTTTTTATAAAAAATTTTTTTTTGGAAAATA
TTTTTAATAACATTTTTTTTTTTTAAATTAACAAATTTAATTAAATTTTTTTTTTTGTAATAACATATTAACTGTTTTCTTTTAAGTTTTATAAAAATAT
TTTTATATAAAATTTTTTTTTTTTTTTTTAGAATAAAATTTTTAAAATATATTGTTATTAATATAAAAAAATATTTTTTTTTAATTTTAAAATATATTAT
TTTAAAAAAAAAATTTAATATTAAATAACTGTAAAAAATAATTTTGCTTTTTTTTAAATTATTATTTTTTTTTTTAACATAAATTGAAAATTCATATAAA
TAAGAAAATATTTTTATAAATTTTTTTAACTTTTTATTTTTAAGAATAGATAGATTAATCATATATATATTTAATATGTTTTTTACTTTTTTTAATCATG
ATTAATCTATCTATACCTATACCTCCTGCAAATCCCAAAGTTTTTAATTTATTATTTATTAAAATATTATTGTTAATAATTCCTATTCCTATTATTTCTA
TCCATAAATTATTATAAAATATATCTATTTCATAAGAATTTATTGTAAAAGGAAATCTTGTTTTTCTTATTTTATAATAAATTTTTTTTTTTAAAAAAAA
AGATAAAAATTTTAACAAAAAATATAAAACATTTTTTAAAGAATAATTATTTTTAAGTATAAAAAAATCAAGCTGGAAAAACTGAAATAAATGAAATTTA
CTAAAATCATTACGATATACTTTTCCAATATTAAAAATTTTTATTTTGTTAGAATAATTTTTAAAATATCTATTTTGAGAACATGAAGTATGTGTTCTTA
AAATTTTATTTTTAAGATAAAACGTTTCTTCAATATTGTTATTTCTAAATATGTTTAATAGTTTGTAATTAAAAAAATAGCTTTCTATTTCTGGTGTTAT
AATTTGACAATAATTATTTAAAATAAAAAATTTTTTAATTTTTAAGATATATAAAGTATTTTTCTTATTTAATACTTTTTGTATTTTTTTATAATATACA
TTATTGTATTTTATGTTTATTAAACAATTGTTTTTAATAATAAAAATAAAATTATGTTTTTATTTAATTTACAATTGTTAATTTTTAAAATATAAACAAA
CTTACTATATGATAAACCAAAGTAAATTCTGCTAAAAAAATTTATTAAAGCTATTTTTTTTTTTAAAATTAGACGCTTTTTATTTTTTTTATCTATATAT
TTTTTAAATAATGATTTTATATAATATTGCTTAGCTAATTTTAAACAATTTTTTTTCCTTCCGTAAAATCCTTTATTATGTTTTAAATATTTTTTTGTGT
TTTTATTAGGTTTAGAACGTGTCATAGATATTTTTTTAATCTTAAAAAATTTGTTTTATTAATAATTTTTTTATTTAAACTTTTTTTTTTGCATTTATTC
AATAATAAATGAGTTTTATTAGACTTGATACATTTAATTTTATTATTTATTAATAACCTTTTTTTTATTGATTTTTTTTTTTTTTCTTTTATCATTTTTT
GGTATAAAATATGAATATAAATGTTTACCGTCTAACTCAATATCAGAAAATTTAAAATATACACCTTTTATCTCGTTTTGCAATCTTAAAATTAATTCAA
TACCTTTTTCCTTATAAATTATTTCTCGTCCTTTAAAAATAATTGTAATTTTTATGCTATATCCTTCTAGTAAAAAACAGTTTGTTTTTTTAATTTTTAA
TTTAAAATCTTGCGTATCAATGTTTATTTTTATACGAATTTCTTTTATTTTTCCAAATCTAGATTTTTTTTTTATTTTTTTTTTTTTTTTTTTAAAATTA
AATATGTAAATAAAAGTTTTTTTTTTTTTGCTAAAGAAAAAAAAAAGTTTTTTTTTAATTTTTTTTAACAAAAAAATAATATTTGTTATAAATAAAAAGT
TTTTAAAAGAGAAAAATTCGTTTTTTTTAACTAAACTAACATTATATAAAACAAGTTTTTTAATCATTAATAATTTAAATTATAATTTTGTAATCTATTT
TCAATTTCATACCAATTTAAATATTCGAAAAAATTTTTAATGTATTGTTTTTTATTATTATTATAATCAATATAATAAGAATGTTCGTGTAAATCTATTC
CAATAAGAGGTATTGAATTAAAACCACCTAAAGATATTGAAAACATTGGATTATTATTATCTTTCGTTTTTAATAAACAAATTTTTTTTTTATTAATAAC
TAACCAACCCCAACCTATATGATCAACAAATTTATTTATAAAATTAAATTTAAAATTTTCATAATTTATATAATTTAATTCGATTATTTTTTTAATATTA
CCAAATATAAAATTTTTATTAAATGTTATATTTTTAAAATAATAATCATGATTTAAATAACCACCGAGCAAACTTATATATTTATATTTATCAATATCAG
ATAAATCATGAATTAAATCAATTAGTTCATATTTGTGATTATAAATAATTTTTTTTTTTAATAAAAAATTGTTTAAATTAATTTTAAAATCATTATAAAA
ACTAAAATGTAAATTAAATTGATTTAACGAATATAAATTATTTTTTTTGTATTTTAAATTAAATATTTCTAACATAATATATGTATGTTTTTATATATTC
TAAAAAAACTATTACAATAAAATAGTTTTAAATAAATTTCGATAAGAATTGTTTTCTGCTTTTTTTCTATAAAAACTTAATCCAAAATGTTTAACATATT
TTTTTTTAAACTTACAACATAAAATAAACGTTTTATTATATTTGAAATTTCTTATTTCAGATAAATAAATTATTATTTTTTTGTCTAAATTTACTTTTTT
TAATTCATTTTTTAAATTTATTTCAAAAATAATTTTTATATATTTATTATAACTATTTATTATATAATTATAAACAAATTTGTTATTAAAGTAAATGTAT
GTATTTGCTAAAAAATAAAAAAAATAATATTTATTCATTTTAATAAAATTAAATAAAAATAATTTTTTATTTTTATAAAAACTTAATTTATTAAAATAAT
AAAATTTGTGTATTAAATAAAATAAATTAAAATTATATAAATAATTACCAATTTCTAAATTTGTTTTATAAATTTTAAATTTCATAATTTAACATTAATA
TTATTTAATAATAATATTAGTATTTTAGTAAATTTTAAATTTATACATGAATCAGTAATTGAAACTCCTATAAAATTATTATTTAAAGTTTGTTTTCCAA
AATTTAAATTTGATTCTAACATAAATCCAATTATATTTGTTTTTAAATATAAAAATTGATTTAAAACATTTTCATAAACGTATAATTGATTTATTGCATA
CTTATTTGAATTCATATGAGAACAATCAACTATTAATGGTTTAATATCTAACTTATTTAAATTAACAAAATTGGGTAAATTACCTCCTCTTAAAACAAAA
TGACAATTTAAATTTCCTAAAGAACTAGTAAATTTTCTATTATTAGAAAAATCAATATAATAATGTTTATTTGAAATAGCTATAAAAGTATCTTTTAAAT
AAGATATTTTACCTGATAACTCATTTTTAATTGCAACTATACATTTTAAATTAGAACAATATTCTCTATGAATTTGCGATAACATTGTTCTTGCACCTAA
ACAAACCCAAAAAATTAAATCTATAAAATAATTTGTTAAATAAAAATTTAAACACTCAACACCTATTAACATATTTTTTTTTATAATATTTATCATTAAA
TTTCTCAAAATATAAATTGAATCTATTATAGAATAACTATTATCTAGATAAGGATCGTAAATATATCCTTTCCATCCAATACTTGATCTTGGTTTTTCAT
AATAAATTCTAATTATTACAACATAATTTTTTGTAAAAATTTTTTTTAACTTATTAATATAAAAATAAAAATCTTTTATATTACTAACAGAACATGGTCC
TATTATTATATATAATTTACTTTTTTTTTTTGAAAAAAAAACTTTTGGTAAATCAATTATAACTCTTTTAACACTTATTAAACTAAACATCTATGTTTTA
TAAAATAATTGTTATATTGATTTTTTGCAATAATTTTATTTTTAAGTAAAAAATTTTTTTTTTTAAATCTAGATTTACTTTTTGATTTCTTTTGTTTTGG
AACTGCCATATGTTTAATTTTATAAATTTATATATAAAAACTTTAGGATGTAATATTAATACATACATTAGTAGCAAAATAATATATAATATAAAATATT
TTAAAATTAAAATTTTAAAAAATTTTTTAAAATCAAATTTATTAATATTAAATAGTTGTGTTGTTAGAAAAAATCCTCAAATAAAAATTTTAAAAGAATT
AAAAAAATGGTTTTTTATTAAAAAATATAAAAAAATAATAATAATATTAACTGGATGTTTAACAGAATTTGAAAAAATTAATAGTTTAATATCTTTAAAA
ATAGATATAGTAATAAACTCATTATCATACATTTTTATAAAAAAAATATTAAATTTGTATTTAAAAACTAAAAAAAAAATTCTTTTGATTAAAAAAAAAA
ATAATTTTAATATTAAAAAAAATATTTTAAATTATATTTCTATAATGAAAGGATGTAATCATAGTTGCACATATTGCATTATACCACAAACAAAAGGTAA
AGAATTTTACTATTCTTTTAGTTATATATTTAATTATATTATTAATAATATTAAAAAAAAAACCACTGAAATAACATTATTAGGACAAAATGTTAATTCT
TATTATAATAAAAACGTAAACTTTAATTCTTTAATTTTTAATATTTCAAAAATTAAAAATATTAAAAGGATAAATTTTCTATCATCTAACATTATTGATT
TTAACAAAAATTTTTATAATTTATACAAAAACGTTAAAAAAATATCAAATCACATTCATTTACCTATTCAAAGTGGATCAAATTTAATATTAAAAAAAAT
GAATAGAAAATATAATTTAAATCATTATATTTGTTTTATTAAAAAAATTCAAAAAATTAAATTTACAACTTTTTCTACTGATATTATTGTTTCATTTCCA
AATGAAAATTTTTTTGATTTTGATCAAACTTTAAAAGTTTTAAAAAAAATAAAATTTTTAGATATATATTATTTTTTATATTCAAAACTTAGAAACACTA
TTTCTTTTAATTTTAAAGAAAATAGTTTTTTTGTTAAAAAATTTAAATTATTTATTTTTCAAAAAAGTATAATTAAAAATTATTATTTATTAAATAATAG
AGTCGTCAGAATTTTAGTTATTGGATATATTAGTAAAAATATTTTTATAGGCAAAATGGATAATTTAAAATTAGTTTTTTTTGAATATTATAAGTATAAT
ATAATTGGAAAATTTATAAATGTTAAAATTATAAAAATTAAAAAAAATATTTTTTTAGGATTATATGAAAATATATATCCTTGTATATAGTTTATTTAAA
ATTAAAATAAACATCAAAAATAAAAATAAAAAAAAAAAAAAAATTATAATTATAATTAAAAAAAATAAATTTATTAAATATTTTAAAATTAAAAATTCTT
TTTTTATATCAAACCTTATAAAAAAAAAAAAAGTTATTAATGAAATTATTATTTTAAAACATTATTGTAAAAATAAATTATTATTACATTGTTTTAAACA
TTTTTTTGAAAAACATCATAATTCATATATTGAAGTTAAAATAAATGATTTATTTTAAAAATCCAGCTTTTATAGAAAATAAAATTACTAAAAATATTAA
AAATATTAAAAATGAAAATTTTTTTTGTATACCAATGTTTCCATATCCAAGTGGAAAACTTCACGTTGGACATGCAAGAAGTTATATTATATCCGATGTT
ATTTCTAGATATAAAAAATTAAAAAAAAACAATGTTTTACAATCTATAGCATGGGATGCCTTTGGATTACCAGCAGAAAACGCTGCTATTAAATACAATA
TAAATCCAGAAAAATGGACAATTTCTAATATTAAATTTATGAAAAAACAACTAAAATACTTTTCTTTAGATTATTCTAATCTAGAATTTTCAACATGCGA
TATTAAATTTTATAAATGGGAATTTTTTTTTTTTTTACTTTTATTTAAAAATAATTTATTATACAAAAAAAAAGAATATGTAAATTGGGATAATGTTGAA
AATTGTATTCTTTCTAATGAACAAGTTAATAATAATAAAGGATGGAGGTCAAATTTTCCAATAAAAAAAGTTAAAATAAAAACTTGGTTTTTAAAAATTA
AAAAATATTCATCAAGATTATTATATGATTTAAATTATAATAATTGGTCAAAAAAAGTTAAAAAAGTACAAAAACAATGGATTAAAATATTTTTTTTTTT
TTTTTTAAAAAAAAAAAAAATTTATTTAAATATAAATAATAAAATTGTTAATTATAATGAAAAAATTTATTTTAAAAGTAAAAAAATAATTTTTTTAATA
ATTAAAAGTTTTAAAACAAATATTATAAGCTCAATAAAAATTTTTTATTATGATAAATTTATTAAAAAAAATAGTAAAATTATAATTTGTAACTATTTTA
TAAAAAAAAAAAAAATAAAAATAGAGTATTTAAATATTCTTATTAATAATAAAAAATGTTTTTTTTTAAAATTAAGTAATTTAAAAAATTGGTCTTTTTT
AAGAGAAAGACGATGGGGATCTCCTTTTTTTTATAAAAAAATAAAAAATAATAATTTTAAAAATTATAAAACAGTGGATACTTTTATCCAATCATCATGG
TATTATTTATTTTATATAAAAACTAAAAATATTAATACAAAAAAAAAAAGTTATTTTTTACCAATAAATTCATACATAGGAGGAATAGAACATATAAACT
TACATTTAATATATTTAAGATTTTTTAATAAAGTTTTTTTTGATTTTAAAATTATTAATGTAAAAGAAGTAATACTGAATCTTATAAACAATGGATTAAT
AAATAATAATGTTTATTATAAAATTAAAAAAAATAAAATTTTATTTTGTAAGTATAACAAAAAAGCTATTTTATTCGGAATAGAAAAAATGTCAAAATCA
AAAAAAAATGGAATCAATCCAATAAAAATTATAAAAAAATACGGATCTGATATTTTAAGATTATATTTTATAACAAATAAACCTATAAATAAAAATATAA
TATGGAATAATTGCAATTTTATTGATATTAAAAACTTTATATTAAACTTAAATAAAAATATAATGTTATTAGATAAAAAAAAAAGTAATATATTTTATCT
AAATAATGTTTTAAATATAAAAAAAATACATACAATAATATCTACAATAAAAAAAATATTATTAAAAAACAATAGTATTATTCAACTTAAGATTATCATA
TATTGTTTATATCCTATTATACCTAATTTATCTAAAATTTTTTGGTTTAAAAATGGTTGTAAACACCCTATTGAAAAATTTAAATTATCTTTAAATTATA
ATAAATTATATAAATTATATTATAAAAATAATTTTATTAAAAAAATTAAAAATTTAAATTTTTTTCTAAACATTAAAAATATGTTTCATAAGATATCAAA
AATAATTATTTCAATGGATGAAATATCTATTATTATTTTATAAATAAAATTTCTTCAAGAGTATTTAAAACAATATGTTTTTTTTTAATTTTTTTAATTT
TAATTAAAACATCATTATTTACTCTATTTATTGCACCAGCGTGACCCATTTTTTTTTCTAATGGAGAAAAAATTCCTACAATATAAAAAAAGATTTTTTT
ATTTACTTTACAAGTAGTAATAACATTTTCAAAATTTCCTCCTATTTCACCAATAATTAATATTTTTTTTGTAAATTTATTTAAAAAAACTAACTTAAAA
ATATTTTTAATATTACATCCTGAGATTATATCTCCACCTATACCAATACAAATTGATTGTCCAATTATTTTTGAAGAAATTTTAATAGCTTCATACGTTA
AAGTTCCAGATTTTGAAATAATAGCTAATTTTCCTTTTTTAATAATACTAATAGGTATTATACCTAGCCTAATTTTTAAAAAAGGTAAAATTAAACCAGG
AGAGTTAGGACCAATAAAAATTATTTTGTATTTTTCACAATAATATTTAATTTTTAAAATATCAAAAACAGATATATTTTCTGTAATGCAAATAATGATT
TTTATCCCACTATAAATATTTTCTAATATTGTTTTTTTACAAACAAAATAAGGTATATATAAAACACTAATTTTGCAATTATTTGACTTAACAGCTTTTA
TAGCAGAACTATAAATAGGAATATGAAAAACAATTCCCCCATTTTTACTATTATTAATTCCGCATTTTATTTTAGTACCAAAATTAATAGAAACTTTAGT
GTGCTGAAACCCGAATTTTCCAGTTATTCCATAAGATAAAATATTATTATTTATCATCTCGTAATTTTACATTTATAAAATAACTTTTTGTAATTTTTAA
TTATTATTATTTTTTTAAGAATAATTATTTTTTTTAATGCTAATTTTGATTTTAAACCATTTAATTTAATAATAATTTTAAAATTAAAATCATAATAAAA
TATGTTTAAAAAACCAATAAGTATTTTTTCACAAGATACAATTCCACCAAATATATTAATTAATAAAAAATTTAATTTCTTATCGATTAAAATTACTTTA
AGTAAATTATTTATATTAACTTCTGTTATCGAACCGCTTAAATCAATAAAGTTTGCACATTTTAAATTATTAAAAGATAAAATATCTAAAGTTTTAAGTG
CTAATCCAGCACCATTTACAATACAACAAATATTACCTTTTAATTTTATATAATTTATTTTTAAAATGTTAGAAATATTGTCCGAGAAATTTTCTTTTTT
AATTTTTTCTGTTATAATTTTGCAATCTAAAATATAAATTTTTTTTTTATAGACTATTAAAGGATTAATTTCTACTAATAATAAATTATTACATATTATA
ATTTTGTAAAGTTTAAAAACTAAATTTAAAATACTTAAAATATTTTTTTTAAAAATTTTGCAATTTAATAAATAATCATATATTGAATAAAATACTATAA
AAGTTTTTATGTTTAATTTTAAAAAATAATTTTTATTTTCTAATTCAACATTAATTCCACCCTTGTTACTAATAATAATAGTTAAATTATTATCATGAAC
AAAAAAAGATAAGAAAAATTCTATTTCAATTTTTGAAATTTTTTCTAATAAAAAATACTTAATTTTTTTATTATTAAATTTTGTTTTTTTCCATTTTTTA
TAAAATACTTTTAATTCATTTTTTGAATGAGGTATTAATATACCTCCTTTTTTTTTTCTATAGTTTAAATTAACTTGTATTTTATAAATATAGTTTTTAA
CGTATGTATAATTACTTGTTAAATATGAATTTAAAATTGGTAAATTATATTTAAACAAAATTTTCTTTGATTCAAATTCACATAAATTCATATTTGCGTA
AATATATTTTTTATGTCTATAAATTTTAAATTTGATAAAATACTGCTTATATCAAATAGATTTTTTTCTAGAAAATTAAGATTTGGAATACTTAAATTTA
AATTTATTTTTTTAACGTGATTAAGTACAGTATGATGCATTTTATAAATATTTTTAAAAATATTAGGATCTAAATAATTTATTATATTAATTATATTCAT
AATAAAATTGTTTTTTTTTATTTTTTTTTTTTTAAATAAATAAGATATTTTAATATCATTATTATTATTATAAATGTTGTAAAACATTGTATTATTTATA
TTTTTTTTTTTTATAAAAAAAATTTTTTTATAAATTTTATTAAAATTTAAATTTTTAGCTTTATTTCTTTTTGACAATTTTTTATTATTTATATTTATAT
TAGGTAAATGATGAAAAATAAAATATTTTTTAAAACAGTTTAAAATTAAAACTTGATTTTTTAACTGATTTAACCATTCTTTACCTCTTATAATAATAAA
AGTATTGTTAAATATATCATCAATAATAGACGAAAAATTATATGTTGGTATACCATTTTTTTTAATTATTATTTCTTTATTATCTAAAATAGTACTTTTT
ATGTTTTTGTAACTATTATCATAAAATTTAACTTTTGTAAAATTAATAAATAAAAAAAAAGAATAAGATAAACTAAAATTATTTAATTTAATATTAAAAA
AATTTTTTTTAAAAATTATTTTTTTTTTAAAAATATTTTTTTTATAATAAAAATTTAAATTAGACGTTTGTTTTAAAATTTTTTTAATATAAAAACCTTG
TTTTTTTAAATTTTTTAATATATATAATTTATTTATGTTTTTATTTTTAATTTTATTTGTATCATCAAATCTTAAATAAACATTTCCTAATAATTTATTT
TTTAATATATAATTTATGAAAATTATAAAAAAATTTCCAATGTGCGGAACTCCACTAGGAGATATTGCGATTCTTAAATTAAACACTTTATATTTTATGT
GGGTTTAGTAGGAGTCGAACCTACAACTAAAGAATTATGAGTCCTCTACTCTAACCATTGAGTTATAAACCCAATATTAAATATTTATTTCTAATTTTTT
TTTTTTTTTATAAAAATCCATGATTGTTACCATTAAAATTTATTTTTTTTTATAACGTATGTACATAACATCAGCACTACCTTATATAAATAATATTTTA
CATATCGGACACATTTTCGAAATGTTTTACGCTGAGTATAATTCATTGATTTATAATAAATTAAATAATTTTAAAGTTTTTTCAGGATTAGATTGTCATG
GTTTGATTAAAAAATCAAACTTTAAAAAAATATTTAAATTAAATATAATTAAAATAAATTATTTTAACTTAAATATCGATTTTAATAAAACTATTACATT
AATAAATAAAAGAATTTGTAATTGGATATATTTATTTTTAAATGATAATAATTATTTATTTGGAAATATAAATAAACAACTTTTTAATAAAGAAAAACGA
TTTTTTATTCCGGATAAATATATTAATTATATATGTTTTTATTGTAAATCAAAAATAAATGATTTTTGTTTTAAATGTAAAAACCAAAAATTTTTACTAA
AAATAAAAATTTTAAAGAAGAATATTATATACAGAAAAACATTTAATATATATTTTAAAAATTATAAATTTTTAAATTGGAATATTTCAAGAAGTAAAAA
CTATGTCGGTTTCTTAATATTATCAAAAATTAATATTTATTTTTATGTTTGGTTTGATGCTTTAATTAGCTATATATCTAATAATTTAAAATTTATTAAA
AAAAAATTTTTAAACAAAAAATTAATACAGATTATAGGAAAAGATATTCTTTATTTTCATAAATTATTTAGAGTAATTTTAAAAATTATAAAATTTAAAA
ACAATAAAATTATAATTCATGGTTTTATATTAATTTTAAATAATAAAATTTCAAAATCAAAAAAAAATAATTTAGAAAAAAAAATAAACGTGTTTTATTT
TAAATTATATATTTTATTAAAAATAAAAAATAAAATAAATGATATTAATTTAAATATTAAAGATATAATTTTCTGCAAAAATTTTTTTTTTAAAAAAATA
ATTAACTTATATTTTAGAATTAGGACAATTTTAAACAAGTTTGATAATAAAACATCAGAATATTTCTTCGTAAAAAAACATCATATTGAATTATATAGTT
TTTATAAATTAAATATTTTAAATAAAATACCTAAAAAAAAAGTACAAGAATGTATAAATTTAAATAAAATATTAGAAAAAAATATTTTTTGGAATAATAA
AAATTTATATTTAACACAAATAAAATGTACTTTTTATATGAAAAAATTAATTTCAATAATTAATTTTTTTTATTTTATAATTAATAAAAATAAAATTAAA
AAAAAAATATTAATTAATAGCAATTTATTTAATATAATTAAATTTTATGAAAATTAAAAATTTTTTATTATATATTTTTAAAAAAAAAAAAATTAAAAAT
TGTAAAAAATGTAATATTTGTACAAAAATTTGTCCATTAAATTTAATTTTAATTATTAAAAACAATATATTTAAAAACTGTAAAATATGTAATTTTTGTA
TTTTAAATTGTCCACAAAAATGTATTAAATAATTAATCTTAAACTTTTTATATTAATAAACCCAGATGAATCAGATTGATTAAAAAGATTGCTAATCTCA
TCAAAAGAAGAATTTTTTGAATTAAACAATGAATTAATAGAATTTATACTTGCAATATTTATTTGACCTTTAAAAATTTTTAATTTAATAATTCCATTTA
TACTAGTTTGTGTATAATCAATAATTTTTTGCAATAAAATTCTTTCAGGAGACCACCAATAACCATTGTATACTAACTTTGAATATTTTAATGCTATTTC
TTCTTTAAACGAATAAATTTCTTTATCTAAAATTAAAGATTCTAATTTTTTTCTTGCATACATAATTATACTTGCTCCAGGAGATTCATAACAACCTCGT
GATTTAATTCCTATTAATCTATTTTCAATGATATCTAATCTTCCAATACCTGCTATTGATCCTAAATTATTTAACTTTAAAAATAATTCTTCAACATTAT
AATTTTTGTTATTGATTTTAATTGGATCACCATTTTTAAACGTTAAGCTTATATAGATCGGGTAATCTAATGAATTATAATTTGAAAGAGTATGTTCCCA
CATTGGTTCATCTGGTTCATAATTAATATTGTCTAAATTTCCTCCTTCATATGAATTATGAAATAAATTTTTATCAATCGAATATTTTTTTGTTTTACTA
TCAAATTTAATATTATTTTTAATACAAAAATTTAATAAAGAATTTCTTGAGTTTAAATTCCAAATCCTCCAAGGTGCTATTATTTTAATTTTAGGATTAA
AATACTTAAATCCTAATTCAAATCTAATTTGATCATTTCCTTTTCCAGTTGCGCCGTGAGATACATAATTTGTATTTAAATAATATGATATTTTCATTAA
TTCTTTAACAATTAACGGTCTTGCAATCGCAGTTCCTAGCAAATAATTATTTTCATAAGTTGAACTTGATCTTAAAAAAGGAAACACAAAATTTTTTATA
AATTCTTTTTTTAAATTTTTAACAAAAATATTTTTAATATTTAAGAGTTTAGCTTTTTTTTTTGCTAGAAGAATTTCTTCTCCTTGTCCTAAATCAGCAG
TAAAAGTAATTACTTCAAAATTTAATTCATTTTGTAACCATTTAACAATAACAGAAGTATCTAAACCACCAGAGTATGCTAATACAATTTTTTCTTTAAT
TTTCATTTTAATTTTTGACAATTATATAATATTAGTTTAACATTAAAACGGCTGAATAACATATAGGTTATGTCTTAGGTTGCAAACCTAATTAAATTGG
TTCGAATCCAATTTCAGCTTGCCCATATGGCGAAATGGTAGACGCAAAGGACTTAAAATCCTTGGTTTTCTAAACGTGTCAGTTCAAATCTGACTTTGGG
TAAAATAAAAATTATTTATATCAAAAAAGTATATCAAATTATTTATTAGGCTATTTTTAAACAATATATATTTTATAAAAATAGGATTTTTGTTTAAAAT
AAAATTTTTAATAAATTGTATATTAAATGTTTTAATATACTTAAAGTATTTTGGTAAAAGATTATTAATCACAGAAATAGATCCTGAAAAGTTTAGTTTA
ATATGGTTTAAAATGGTTGAATCTTCACCTGCAAAAAAAATTATTTTAATTTTTTTAATTAATAAATTATTTTTTTTATTATAAGAATTTTTAATAGAAA
AAATTATTTTTGATAATTTTTTCATAAAAAAAAATGAAATATAAATTCCAGTTCTTTTCGGTATATTATATAATATTATAGGAATACCTATTTTACTTAT
TATTTTATAATATTTGTAAATGTACAAATTATTAGGTAAAATGAAACTTATAGGTGATATTAATATTGCTAAAACTTTATTTATTTTCAAAATAAAACAT
ATATCTACAATATCATTTATATTATTTTTATTAATACCAAAAAAACAATTTATTTTATATTTATTAACATATTTAGATAAATCAAAATATGTTTTTTTGT
TAAAAAGATTGGATTCTCCAGTTGTTCCTAAAAATAAAATATTGTTATTTTTATTAAATAAATTGTATATTATTAATAAATTAAATAGTTTCCAATTTAT
TTCTCTCTTTTCATTATAAGGAGTAATTATTGCAACAAAATTATACATTTTATATATTAGTATAAGTATAACAAAAAAAATGAGTTGTATAGGATTTGAA
CCTACGACCAATTGGTTAAAAGCCAAGTGCTCTACCAGTCTGAGCTAACAACCCTTATAAAATTTAATTTTTTTTTTATATTTTATAAATAATAAAAAAA
AGTTTTTTTTAATTTTTTCTTTTAAAAAAAATTTTTTGCAATAAAAACATTTTTTTTTATTAACTATTAAATCTAATTTCAAAAATATTTTTTTATTTTT
TTTAAAAAAATTAAATATTGTTATCTGTTCTACATGAAATCTATTTAAATTTGAATGTAATTTGTTTGCTAAACTAATAACTTTAATTGTATAAAAATAC
AATAATAAATATTCTACTTTTTTTAAAAAATTAAACCTATATTTTAAAAATGTCATATTTTATTTCTATAATTGGTTCTGGTATTATAGGAGCTATAACA
AGCTTAATGTTATTTATTTGTAATAATAAAAATATTGTTATTGTTTTTGAATCAAATAAAAAAGTTTCAATTGAAAATTCAAAAACATTAAATAATGCTG
GTACTGGACATGCAGGTATGTGCGAAAATAATTATGTAATTCAAAAAAAAGAAAATTTTTTTATTAAAAAAAATATTAGAATATATTGTAAGTTTGAAAT
TACAAAAATATTTTTTTCATGGATTAAATATTTAAAAATTTTTAATTTTAAAAAAAGTCTAATAAAAGTACCTCACGTTTCTTTTTTTTTTTTAAAACTA
AATAAAATAAAATTAAAAAAAATTTTTAATAAATTAAAAATATTTTCGAATTCAATAAAATTTACTAGTAATATATACTATATAAATAAAATATATCCTT
TATTATTAAATAATAAAAAAAGTAAAAAAAAATTTACAATTACATATTATAAAAATGGATTTGATATTAATTATCGTTTAATTGTAAAAAAAATTTTTTT
CTTTTTAATAAAGCAAAAAAATTTTTTTTTATATTTAGAAACTGAAGTTTTAAAAATAAAAAAAAAAAATTTTTTTTATTCATTAAATATTAAAAAAAAA
AAATACTTATTTGATTACGTATTAATATGTGCTGGAGGAATGAGTTATAATTTAACTATTGAAAATAATAAATTAAATTTAAATAAGTATTTAAATTTTC
CTATTAAAGGTAATTGGTTAATAAATGAAAAAAAAAAAAATGTTAAAAATCACAATATAAAAGTTTATAGTGAAACTATAAAAAATAATCCACCAATGAG
TACTCCTCATTTAGATTTAAGAAATATTTTAAACGAAAAAAAAATTTTATTTGGACCATATGCAGGTATAACTTTTAATATATTAGTAACTAAAAGAAAA
TTTATATTTAATGATTTAAACATAAAAAATTTTTTTTTAATTATTTTGTTTACAATTAATAATAAAATATTAACTAAATATCTATTGTTTGAAACAATTA
GTACAAAACGAAAAAAAGTATTAAATACTTTAAAATTTTGCAATGTTAAAAAATTTTATTTAAAAAATGCAGGTAAAAGATTACAAATTTTAAAAAAAAA
AAATAATAAAATTGAAATAATTTTTGGCACAAAATTAATTTTTGATAAACATAAACATTTAGCAACAATTTTAGGTGCATCTCCTGGTGCATCAATATCT
GTTTATATTGCTAAAAAACTAATTAAAAATTGGATTAAATTTCCTAAAAAGTTTTTACCTAATTGTAAAAATTTAATAAAAAAAAACAAAATATTTAGTA
AAATACTATACATTTGAAAGGAATTGAACCTTTTACAATCTGATTCGTAATCAAACACTCTATCCGATGAGTTACAAATGTATAGAGAAAGTAGGATTCG
AACCTACGATAGATTCTATACTCCTTTAGCAGAGGAGTGCTTTAAACCGCTCAGCCATTTCTCTTTACAAAAAATTGTTTATTAAAAATTTTATATAAAT
TTTTAATTGCTAATACAGTAAACTTTTTTTTTATTAAAATTGATATTTTTGTTTCAGAAGTAGAAACCAGAATAATATTTATTCCAAGTTTCGACATTGA
ATAAAAAATTTTTCCGATAATATAATTATGTGATCTAAGTCCAATTCCAATAACTGAAACTTTTGAAATATGTTTTTCATATTCTACTCTTCCTCCTATT
TTTGTTATAAAAAATTTTTTAATTAAAAAAATTACTTTTTTTAAAAAAAATTCTTCAATTAAAAATGTAAAATTAGTAAATTTTAAATGATTTAATGAAT
TTTGAATTACCATATCGATACAAATTCCATTCGATATTATTGGTCCTAATATTTTAGATAAAACACCAGATACATTTGGTATATTAGCAACTGTAATTTT
AACTTCATTAGAAGTATATGATATTCCAGATATTAAAACTCTTTCCATAGAATTAACAAATTTTTTTTTTTTTGAAATAAAAGTACCTTTTTTTTTAAAA
AAAGATGATAATAATCTTATGTTAACATTATATTTTCGTGCAAGTTCTATAGATCTAACAAATAAAATTTTTGAACCTAAGCTAGATAATTCTAGCATAT
TTTCAAATGGAAGAAAATTAATTCTATAATTTAAACAAATTCTAGGATCAGAAACAAAAATAGATTTTACATCAGAATATATTTGACATTCTGTTGCTTT
TAAATATATAGATATTGCAACTGCTGAAGTATCAGAACCTCCTCTACCTAATGTTGTTAAATTTCCATTTAAAGTTATACCTTGAAAACCCGTTAAAACT
GGTATTTGTTTTTTTTTTAAAATTTTTTTTATTAAAACTATATTTTTAATAAGAATTATTCTAGCATTTGAATAATTATTATTTGTAATAAAGCCTATTT
GGCTACTTATTAATGGTATAGAAAGTATCTTTACACTATTTAACAATAAATTAAATAATGAAACGCTCATTTGTTCACCAATACATAACATAAAATCTAT
ACTTTTTAATTTTTTATATTTAAAATACTTTGCTAATCTTAACATTTTAGTTGTTTCTCCGCTCATTGCTGATAATATTATTACTATTCTTATATTTAAA
TTTCTATATTTTTTAATAATATTTTTAAGAGATTTAATTCGTTTTTTTGTACCAACTGACGTTCCACCAAATTTTTGAATTATTATCATGAAAATTTTAT
TTTAAATAAATTATATAAAATAGATAAATAAAATCCAAATGTGTTTTTTAAAAAAACAAAATCTTTTATTTTTTTAATATATATATTTTTTTTAATTTTT
TTTAAAAAAAATATTTCTTTTTTAAGAATAAAATTAATAAATTTTTTTAAACAAAAAACAAAATACACTTTTTTTTTTTTTAAAAAAAAATAAAATATAA
AATACCAAATTTTAAAATTATTTTTATATGAAAAAATAAGTATTTTTTTTATCAACTTATTAATAAAAAAAATGTTTTTTTTATTTTTTAAATAACTAAA
AATTAAAAAAATTGAATTAAAAATATCATACAATTTATAATTAATTTTGTTATTTTTTTTTTTAAAAACATAAATTATTCTATTTAATCCTAATCCTGTA
TCTATAATTTTCTTTTTTAAAAATTTTATCTTTTTTTTTTTTTTAATAAAAGAAACATTAACAATATTCCAAATTTCTAATAATTTTTTACTTATTTTTA
TATATATTTCTAAACTAAAACCGCAATATCCGTTTTTATTTATTTTCCAAATATTTTTATTTGTAAATATAACTTGATTAATGTTAATTTTAAGTATTAA
TATAAGAATTATGTTTAGATAATCATTTATATTAATTGTAAAAAATAATTTTTTAAAGTTTAACTTTATTAAAAATAAAATTTTTTTAATTTTTTTAAAA
TTATTTATATTATTTTTTTTTGTAAAATTTCCTAGCATTAAAAAAGATGTTTGATGAATACCATCATTTGTTAATTTAAAATCATTATATATTCCTTGCA
TTCTAACACAATATTGAAATGATGATATTTCAACATTAGATCTAAAAACAAATTTTTTAAGGCTAGCTAAACCAGAATTTACAAAAATTAAAGTTTTATT
AGTTGAATTTATATTTTGAGTACCTAAAATTTTGTAATAATAAAAACTAAAAAAAAGTAAAATTTTTTTCATATCATTGATCTAATATCAGATTCAATTT
TTTTTAATAATTTTTTATTAAACATTAGTTTTTTAATTAAATTTTTTTTGTTATTAAAATATTCGTTTTCTAAAAAAAATTTTTTTTTAATAATTTTAAT
TAAATTAAATTTTAAACCAAAATTAATAATTTCTAAAATTTTGTATATTCCATAACTATAAATTAAGTAAATTTTAGTATCTCTATACGGTTTTGATAAT
TTGTTTTTCAATACATTAACTTTAATTTCTTGTCCAATTATATTATTCTTAGATTTTAGAAATCCAATTCTTTTTAATTCAATTCTAATTGAAGAGTAAA
ATTTTACTGCGTTTCCTCCTGTTGAAATTTCTTTATTATAAAAATTACTAATTTTAATTCTAACTTGATTTATTAATATTAATAAAACATTGTTTTTTCT
TAACAAAGGTATTATTTTTTTTAAATTTTTCGATAAAAATCTAGAATGAGAACCAATATTATTATCATTATTATAAATTTCAAGTTCAGGAATTATTGCA
GCAATAGAATCAATAATAATTAATTTAACATATGTAGAATTTATTAATTTTTGTGATATTTCAAATACTTTTTCTCCATTTTCTGGTTGAAATATTAATA
ATGTTTTTAAATTAATACCTAAGTTTTCAATATAGTTTACATCTATACAATGTTCAACATCTATATAAGCACATATGTCACCAACTTTTTGTGCTTCTTT
TATTATAGAAAAAGCAAACGTAGTTTTTCCGGAAGATTCTTGTCCATAAATTTCAATAATTCTACCATAAGGTAAACCACCAATTCCTAAAATAAAATCT
ACGTTTAAAGATCCAGTAGATATAAATTCTACATTTTTAAGATAATTAATATCATGCATAATAACAATTTTGTTTTTAAAAAAACTTTTAATATTTTATT
ACCTCCCCCAAATATAACAGGAGAAGTACATATTGGACATTTTTATCAATATTTTATTATTGATTTTATTACAAAGTGGAAATTAATACAAGGATATAAA
ATTATTAATAAGTTTGGATTTGATCATGCTGGTATTTCAGCTATAATTAAGTTTAAAAAAAAAAAAAAAATTTTAATTTTTTTAAAAAAAATAAAAATAA
GTTTTCGTAAAAAAATGTATTTTATAAATTTTATTTTAAATAAAAAAATCGAATTCACTTTAAGTAAAGTTTATAAAAAAGTAACAAAAAAAATTTTTTA
TTATTTATTTAAAAATAACATAATATATATTAAAAAAAAAAATATTAATTTTGATTATAAATTAAAATCTATATTATCAGATATTGAAATATCGAAAAAA
ATTTATAGAAAATTTTTATTTTTAATAAAGTATAAGTTGAACAATATAAATATTATTGTTCCAGTATCAAATATATTTTCAATAATAACTAATACTGGAA
TTATAATTAATAAAAATATAAAAAAAAATTCAGTTGCTTTATCTCCTTTTAAAATTAAAGTTAAAATTATAAAAAAAAAAATAAATAAGTTTAATTTCAT
CAAAATATCACCAATATTTAATAATTGTGATTATTTATTAAGTATTAACAATAAAATTGAAATAATAACACTATTAACAAAAAAAAAAAAAATAAAATTA
TTAAATTACAAACATTTAAATAATAAAATAATTGAAAATAAAAAAATTAAAAATTATAAAATTAAAAATAATTTTTTAAAAAAAAAAATTTTAAAATATT
TATTTTATAACAATTATATTATATGTATTAAAAAAATAAAATCATACACAAATGTTAATAAAAAAAACAATTCTAAAATTTTTTATTTATTGATAGATCA
ATGGTATTTAAAAATAAAACACATATTTTCAATAAAAAAAATAATAAACAAAATCTTAATAATACCAAAAAAATATAATAAATTATTAAATAATTGGATA
TTAAATTTATCAGATTGGTGCATTTCAAGACAAATAAATTGGGGATCTAAATTTCCTATTTTTAAAGACCAAGAAAAATTTATTTATTTTAAAAAAACAA
GAATTGAAAAATACAAAAATTTAAACGAAGTTTTAGATACTTGGTTTAATTCTTCAATTTGGAGCATTTATATTTTTAACAAGAATAAAAAAAATCAAAA
TATTTTAATTTCTGGATTTGACATTATTTTTTTTTGGATTTTAAAAATGATAATAATAAACATATATTGTTACAAAAAAATATTATTTAAAAAAATTTTT
TTACATGAAATAGTTAAAGATTATAAAAATAAAAAAATTTCAAAATCAAATAATAATAGTATTCCATTTAATATATTTAAAAAAAAAATTAATAAATACA
AAAATATTTTTATTAATAATATTTCTAAAATAAATATTTTTGAAAAAAAAAAAATTAATTTTTTTATTATTAAAAAAAAAAATATTTTCAGTAATATTAA
AATATTATATTATTTTTATTACTATAAATTTTGCTTGATAAAAAATTTTGAAAATTATAATGTTTTAAAATATAAAAATATTATTAATAAAAGTATTTTA
TCAAAAAAAATATATTTAAAAATACTTAATTTAAGATTTCCTATTAAAAATAAAATTTTAAATAAAAATTGGAATTTTACAATGAATTTTTTATTAAATT
ATAAAATCAAAGCAAATTTTTATTTTATAACTAAGATAAAAAATAAGAATTATTTTATTAAAAACAATTACATTGCAATATTTAATAAATTAAATATTAT
TATTTATGTTAAACAAAAAAATAAATTCATTTAAAAGTATTGCATATTTTAAAAAAAAATTTTTTGAAATTTCTGACAAAGAATTAAAAAACTTTTGGAA
TGTTCTTTTTTTTTTTCCATTTTCATATACTTATATTTGTCCAACTGAATTACTAGAAATATCTAAAAATATAAATTTGTTTAAAAACTTAAAATGTAAT
ATTTATGCAATATCAACTGATAGTCATTATACGCATAAAAACTGGATTGAAAATGAAATTAATTTTATAAATTTTCCTTTTATTTCTGATTTTAATCACA
AAATTTCTAAAAACAATAAAATTTTAGATGAAAAAGATGGAAATTGTTTTAGATCTACATTTATTATTGATCCAAATTTAATTATTAAATCAATAGAAAT
TGTCGATATATCAATTAGTAGATCAATTTTAGAAATTATTAATAAAATTAAAATGTTAAATTTTACATTTAAAAATAAAAACAAACTATGTCCATATAGT
TGGTTAAAAGATAATAAGTCAATTAGTATTTAATTATAAATTATTTATAATATTAGAAATTAATTTAATAAGTTTATTGTGATCAATATCAAATTGTTTG
ATTCCATCATTTAGTAATTTTAATGAAATATTATTTTTTTTATTACAATAATTTTCAATATTTTTATAATTTAAAATATTATTTTCAAAATTAATATTAT
TTAATTTATTATAATATAATGGAGATATAGTAAGAAAATCACAACAAGTTAGATTTAAAATTTGATTTATATTTCTAAAACTTGCTGCCATTATTTTTGT
TTTATAATTATTTTTAATTTTAAAATTAACTAAATTTTTAACAAAATTAACACCAGCATCATATTTACAATAAATTTTATCACTAACTCTACCAACAAAA
GGTGATACAATATATATTCCAGAATCAAAACATTGTTTTGCTTGTTTAAAATCAAATATAAGTGTTAAATTAGAATCTATATTATTACTTTTTAAATATT
TAGCTGCTTTAATACCAGAATTTGTTGCTGGTATTTTAATTAAAATTTTATTTAAGTTAATACCTGATTTTTCACATAAAAATATTATCTTATATGAATA
AAATATAATTTTATTGTAATTAAAAGAAATTCTAGCAGGAATTTCAACTGAAATTTTTTCTTTAATATAAGGAATTAAATTTGATATTATATTAATTAAC
AATTTATCGTATAAACTAATACTATAGTCAGAATCTAGTATATTTTTTCTAGTGTTTTTAAGTAAAAAAAAATATAGTAGTTTATTATAAGATTTATTTA
ATATTGCTTTCAATATTAAACTAGGATTTGTTGTTGCAGCATCGAATAATTTATTTTTAATAAAAAATAAATCAGCACTATCAATTGATATTTTTACTTT
GTTTTTAAAAAAATTAAATAAAAACATAAAAAGAAATTATTTTTTTTTTTTTTTTTTTTGATAATTTTGAAATTTTTAAAAATATATTGTCTAATTATTT
CTTACTCTCACATTAAATGCTACCATAAGCGCTTGATATTTTATCAAAATGAGTTCGGAATGTTTTCATTGTTTTATATCAGCTATAATAATTAGACATA
ATTATTTAAATATATAGATAAAATTCTCACGATAAAATTAGTATGCATTAGCTAAAATATTACTATTTTTACACATTGCACCTATTAAACTTGTTATATT
CAAGTTACCTTAAGAAAGACAACAAATCTTTTGGAGATTTATTTTTGAGATATGCTTCCTACTTATATTGTTTCAGTAGTTATCATTTTAAACTTAGTTA
ACCGGCTATGCTATTGGCATAACAACCGGAACACCAGTGGTTTAACCATCTTGGTCCTCTCGTACTAAAGACAGATTCTCTCAAAATCTCTACATCTACG
GTAGATAGGGACCGAACTGTCTCACGACGTTCTAAACCCAGCTCGCGTACCACTTTAAATGGCGAACAGCCATACCCTTGGGACCAACTTCAGCCCCAGG
ATGTGATGAGCCGACATCGAGGTGCCAAACATTGCCGTCGATATGAGCTCTCGGGCAATATTAGCCTGTTATCCCCGGAGTACCTTTTATTCGTTGAGCG
ATAACCATTCCATTCAGAATTATCGGATCACTAGAACCTACTTTCGTATCTGTTCGATATGTCTATCTTACAGTTAAGCACTCTTATGCTCTTACACTCA
ACATACGGTTACCAAACGTACTGAGAGTACCTTTGTACTCCTCCGTTACTTTTTAGGAGGAGACCGCCCCAGTCAAACTACCCACCACACAATTTTTATA
AATTAGATTATAATTTATATTAGAATATCAAAAATTAAAGACTGGTATTTCAAGAACGACTACTTACAATCTTGTAATTGTAGATCAACATCTCCCAGCT
ATCCTACACAAAAATTACCAATATTCAATATGAAGTTATAGTAAAGGTTCACGGGGTCTTTCCGTCTAGCCGTAGATACACTGCATCTTCACAGCAATTT
CAATTTCACTGAGTTCTAGATGGAGACAGCGTGACCATCATTACGCCATTCGTGCAGGTCGGAACTTACCCGACAAGGAATTTCGCTACCTTAGGACCGT
TATAGTTACGGCCGCCGTTTATTGGGGCTTTTATTAAATGCTTTGTTTTTAAAACTAACATTTTCAATTAACCTACCAACACCGGGCAGGCGTCACACCC
TATACATCCGTTTACACGTTAGCAGAGTGCTGTGTTTTTAATAAACAGTTGCAGCCACCTTACCTCTTAGACTAGTTCTTTAAAACTAAAGTTTTAAATA
TTCTAGTGTATCTTTTCCAAAAGTTACGATACTATTTTGCCTAGTTCCTTCACCTAGATTATCTCAAACGCCTTAGTATACTATACTTATCTACCTGTGT
CGGATTGCGGTACGTGCATTAATATATTAATACTAGAAACATTTCTTGGTAACATAATTCTATTAATTTTATTTATATTATTAAATTATATAACTAGATT
TACTTGGTTAAAAATTTTTTAATATTAAACCTAAAAATCCTATTTTAGGATTAAAGAACTTTTTACATATTTCTTCGTATAAATTAATGGTAAAGAAATT
TTAATCTTTTTCCCATCAGCTACGTTTTTCAACCTCGCCTTAGGTCACGACTAACTCTATTACGATAATCGTAGAATAGAAAACCTTAGATTTTCGGCGA
GAATGATTCTCACATTCTTTATTGTTACTTATGTCAGCATTCGCACTTGTGATATCTCTAAATTATTTTAGAATAATTCTTCTTTGATTTACACAACGCT
CTCCTACCACAAAAATAATGTTCATAATTTCGGTATCTATTTTAGCCCCGTTAAATTTTTTGTATTAGTTTTCTAAATCAATGAGCTATTACGCTTTCTT
TAAAGGATGGCTGCTTCTAAGCCCACCTTTTGATTGTCAAAGAAAACTAACATCATTTTCCACTCAAATAGAATTTAGGGACCTTAATTTATGATCTGGG
TTGTTTCCCTTTTCACAAAGGATGTTATCACCCTTTGTGTGTCTCCTATAATAAAATAATTATTATTCTTAGTTTGTTATGATTCAGTAAAAACTTAATC
AATACAGTGCTTTACCAATAATTATTATTTATAAGGCGCTACCTAAATAGCTTTCGGAGAGAACCAGCTATCTCCGAGCTTGATTAGCCTTTCACCCCTA
TCCACAAATCATCCGAATCTTTTGCAACAGATACCGGTTCGGTCCTCCAGTAAATTTTACTTTACATTCAACCTGTTTATAGATAGATCGCTCGGTTTCG
GGTCTATTTAATTTAACTTTCGCCCTTTCAGACTTGATTTCTCTACGCCTACCTAGTGTTAAGCTTGCTAAAATAAATAAGTCGCTGACCCATTATACAA
AAGGTATATAGTTACTTTTCAGCTTCTATTGCTTTTACGTATATAATTTCAGGTTCTATTTCACTCCTATTTAAAGGTTCTTTTTCATCTTTCCCTCACG
GTACTAGTTCACTATCGGTTAACTATTAGTATTTAGCCTTAGAGGATGGTCCCCCTTTTTTCTGTAAGGATTTCACGTGTCCTAACATACTCATATTATT
AAAATAAAAATTTAAAAATAAAGACTATTACTTTTTAGAGTAAATTATTCAAAATTTTTTTTTTAAATATTTTAAATTTAATTTTTTGGGCTTCTTCCAT
TTCGTTCGCCACTACTTTGGAAATCTCAAATTGATTTCTTTTCCTCGGGTTACTTAGATGTTTCAGTTCACCCGGTTCACTTTTATTTTTAAATAAATAC
TAGTTTTTAAGTAGTAGGTTACCCCATTAAGAAACCTTAAAAAGATTATCGCATATTAGCGTCTTTCATCGTCTTTAGTTACCTAGACATTCGTTATATA
CGATTTTTATATAATTGTATCTATATATTAATATTATCAAAAATTTTAATTATAAAATATACATATTTTTTTTTTTTTAAAAATTTTCCAGCCACAGGTT
CCCCTACAGCTACCTTGTTACGACTTCACCCCAGTTACAAATCATACCGTAATAATACTTAAATTACTTATGATACAATCTACTTCCATGGTGTGACGGG
CGGTGTGTACAAGACTCGAGAACGTATTCACCGTAACATTCTGATTTACGATTACTAGCGATTCCAACTTCATAAAATCGAGTTGCAGATTTTAATCCGT
ACTGAGAATAATTTTGAAAGATTGGCATAATGTTACCATTTAGCATTTACTTTTTGTATTATTCATTGTAGCACGTGTGTAGCCCTACTTATAAGGGCCA
TGATGACTTGACGTCGTCCTCACCTTCCTCCAATTTGTCATTGGCAGTTTCTTATTAGATCTAATATTTTTATAGTAAAATAAGATAAGGGTTGCGCTCG
TTATAGGACTTAACCCAACGTTTCACAACACGAGCTGACGACAGCCATGCAGCACCTGTCTCAAAGCTTAAAAAAGCTTTACTATTTCTAGTAAATTCTT
TGGATGTCAAAAGTAGGTAAGGTTTTTCGTGTTGTATCGAATTAAACCACATGCTCCACCGCTTGTGCGAGTCCCCGTCAATTCATTTGAGTTTTAACCT
TGCGGTCGTAATACCCAGGCGGTCAACTTAACGCGTTAGCTTTTTTACTAAAAATACATAAACTAATTTCTTTAAAATTTTTTTTTATTTTATGTTATTT
AACAAATAGTTGACAACGTTTACTGCGTGGACTACCAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGTATTAGCGTCAGTATTAAAATAGAAAA
ATGCCTTCGCCATTAGTATTCTTCTAAATATCTACGCATTTCACTGCTACACTTAGAATTCTATTTTCTTCTTTTATACTCTAGTATAATAGTATTAACT
GCGTTATTAAAATTTATTTAATAAATTTTACAATTAACTTAATATACAACCTACACACCCTTTACGCCCAATAATTTTGATTAACGCTAGCACCCCTCGT
ATTACCGCGGCTGCTGGCACGAAGTTAGCCGGTGCTTCTTTAATAATTACCGTCAAAAAAAAATTTTTTTTACTAAATTTTTTTCTTCTTTCATTATCGA
AAGTGTTTTACAACCCTAAGGCCTTCTTCACACACGTAGCATAGCTGGATCAAGCTTTCGCTCATTGTCCAATATTCCTCACTGCTGCCTTCCTTAAAAG
TTTGGGCCGTGTTTCAGTCCCAATGTGGTTGTTCATCCTCTAAGATCAACTACGGATCATAGCCTTGTTAAACGTTTATTTTAACAACTAGCTAATCCGA
TATAAGCTCTTTTTTAAGCATATGGCTTTCCCATACTTTATCCATTAGGATTAATAAACTATTAATATATATTTCTATATAGTATTTTTTACTTAAAAGT
AGATTCTTATATTTTACTCACCCGTTCGCTGCTAATACTTTTTTTTAAAATATTCGCACAACTTGCATGTGTTAGGCTTGCTACCAGCGTTCAATCTGAG
CTATGATCAAACTCGTTATTTTTTTTTTTTGAAAAATAATATTTCAATTACAAATAAGTTATTTATTGAATAAATAAAATAGTTAACATTATATATTAAA
TAATTGTTTTATTTTAATTATATTATTTAATATTATATTATTAAAAAATGAATTATTTTTTTTTTTTTTTTTTAAATAAAAATTTAAATTTTTTAAAAGA
TTGATAAAAACTAGTTTTTTTTTTAAAAATTTAAAATAAATATTTTCTACTATATAAATTTTATTTTTATTAAAAGAATTTTTAAAACTAATTTTTGTAT
CTGTTACTAATTCAAAAAAGTTAAATATTGATTTTTTTCTTTTTAAAGTAATAAAAAAATTATTAAAAAACAATAGTAAATTTTTATAATTTGAAAACAA
AAATAAATCGTTTTTATTAGATTTAGACATTTTTTTTTTGTTATAACTGTATAATACTTTATTATTAATTATAAACTTATTTTTTTTTAAAAAAAATACG
TTAATATATTTACTAATTTTATTAACTATTTTATTATATAATTCAATATGTTGAATTTGATCTATTCCAATAAACGTAAAAAAATTATTAATAGTGATTA
TATCAGAACACATTAAAATAGGATAACAAACTTTAGCAATAGATAAGTTAATTTTGTCTTTAACAATTAAAAAATTTTTTGTTTTATTAATATTATAAAA
ACAAGATATTACCCAAAATAAAAATAATATATTTTTATGATTAGATTGATGTAAAAATATTTTTTTTTTAAAAAAAGATTTAAAAATTAAACTAGAAGTT
AATTTGTTTAAATTAACATTTTTGTAGTGTTTATTTTTAGATAAACAATGTAAATCTGCTAAAAAAATATTTTTTTTTTTAATTTTTAAAACTGGTTTTA
TAAGAGTTATATAATTTCCAAAATGAATTAATCCACTTGTATTAATACCTAAAAGAATCATTTTAACTTTATTTCATTAAATAAAACTTTTTTTTTTATT
ATAGGATCAAATTTAAAAATAGATATTTTTTTACTGTTTTTTTTTGATTTTTTTTTAACATAATAATAACTTGTTTTTGAACTTTTTAATTTTATTAAAA
TATACATCTTTTAATTATTTTAATTGATTTTGTAGAAATATTAATTTTAATAAATGCATTTTTTTTCCAAATAGAAATATTCTTAAAATTAGTTTTAGTT
TTTTTTTTATTTTTAATGTTCGAATTAGAAACATTATTTTTTATTAAAGTTTTTTTTTTACATAATAAACAAATTTTAGACATTTAATTACTAAATGTTA
TTTTAATAGAAAAATTGAAAATGTAAACTAAATCTTTGTACGATATTTTATATTTTTTTAACAATTTATTTTTTTTTTTATTTTTATTAAAAATAAATTT
ATTATTTTTATTAAAATAACTTTTTATTTTTACACCATAAAATGGACAGAAAAAAAAATTAATTTTTAAAAAAAAAATATTTTTAAAATTTAAAATAACA
TTTTCTAAACTATTTATATTGTTAGAAAACAATTCTTTTAAAAAAATATTATATTTGAATAAATAATAATTAACTTTAAGTAAATAAATTTTGTAAATCA
TCAATATTATTTAATATACTTAATATTTTTAATTTTATACATAAGCATTAATAATCTATTTATTCCAATTCCAACTCCGTTAAAGTTTTGTAAACCTTTT
TTTATATAAGATAAAAAATCTTTATTTTTTTTATTAAAATTTATTTTTTGTAAAAAATAATCGTTTAATTCATTAAATCCATTTGATATTTCAATACTTG
ATATATACAATTCAAATCTTTTAGCAAACTTAAAATTTAAAAAATTAGGTTTAGATAACAAACTATTTTTTGTTGTGTAAGATTTAACAAAAATAGGTAA
ATTATAATTCTTAATTATTATGTTATCAAACAATTTAAATATTAAACTACTTAAAAAAAAACTTTTTAAAATTTTAATATTTTTTTTTATTAGAAAATAA
AATAATAAATTTGTATTTTTTAAAAATTTTATATTTTTTTTTTTTGCTACTAAACATAATAATTCAATTATACTTTTTTTTTTAAAAAATTGTTTAATTG
AAAAGATAATTTTATAAATATTAAAAAAAAAAATACTATGAATTAAAAAAATATTTTTAATAAATAATTCTAAAAAAGAAATAGAATTTCTAAAAGAAAA
ATTACTAGAATAATATTCTAACATTAAAAATTCAAAATTATGAATATTTGATACTCCTTCATTTCTAAAACAATTACTCAATTCATATATATTATTATAA
TTACAAGAAAGTATTTTTTTAATGCTTAGTTCAGGAGAAATCTTAAGATAAAAAAATTTTTTTTTATAAAAATTAAACGTTTTAAAGTTATTAGAAATTG
AATTAAATTTTTTTGGACAAAGATTATCTGTATTAACATTTAAAAACATAAAAAAATTTAAAAAATTAAAAATTAAATTTTTAATAGAAAAAATTTTTTT
CTTAAAAAAAAATTGTTTTTTAACATTAAAAATACATTTAACAATAACTTTTATAAAAAAAGGAAAATAATAGTAAACTTTTTTTTTTGTTTTTTTTATT
TTACAAAAAGAAGTGATAATATCACCCAATTTAACTTTAGTTTTATTATAAATTTGTATTTTTCCTGAAAAATCTTCAATTTCAATAAATGGATACTTAA
CTCTAGTAATTTTTGAACAAAAAAATTTATATTTTTTACAATTAATAATCATATTTTAATAAAATCAATATTTAAATTAAAATAATCATTTAAATTTATT
TTTTGTTTTTCAACATGATTAATAATTTCATTATTTTCAAAATAATATGTTTTAATATATTTTTTATTATTTTTAAATTTATTTATGTTTTCTTTTACAA
ATAAATTATTTTTTAACAATAATTTATATTCTAATGTCTGCATTGCAAATTTTTTGTTTAGCATTTGAGATCTTTCTGATTGACTTGTAACAACAATATT
CGTTGGTTTATGAACAATTCTAATTGCAGAATTAGTTGTATTAACATGTTGTCCTCCTGATCCTTTTGATCTATAATTTTCAATTATTAAATCATTTTTA
TTAATGAAATTTACATTATTTATAATTTTAGGAACTATATTTAAATTTAAATAAGAAGTTTGTATTTTTTTAGAAGAAATCAGTGGATTTTTTCTAATTA
TTCTGTGTAATCCTGATTCATTTTTAAATAAAAAAAAAGAAAAAAAATTGTTAATTAACAACAATGATTTTTTATATCCGTAAGGACTTTGTTCTAAATT
TATAATTTCTACATTAAATTTATTTTTTAATAACCATTTATAATAAAAATTTAAAAAAAATTTATTTAAATCATGTGTGTCTATACCACCTTGATTTGAA
ATTATTTCAACGAAACAATTTAAATTATTAAATTTAATATTTAAATTATTTATTAAGTATTCACAATTTAATAAAAAAGATTTATAATCAATAATTTTTT
TCTTAATAAATTTTTTTTTAATATAATAGTTTTCAATTAGTTTGTTGTAAAATTCAATCATATTTTTTTTTTTACATAATTTAATATTAAATTATAGTTT
AATGAAAAAATCTTATAGCTTTCAACTAATTTAAAAGTTTTTTTAAATTTAAAATTTAATTTTATACTAGGTATTATTTTTTTTATAGAAAAAAAAAACT
TAATAGGATAAGCTGTTGAAACTACTACATTATTAAATTTTAAATTTTCTTTTATTAAAGATGTTAAACTTGTAATTGTATGCGGATCATATATTTTTTT
ATACTTCAAATAAAAAAATTCTAAAGATAAAATTATATACCTATTATAAATTTTATCAGAAAAAAAAATACTATTAATTTCATTTAATTTTTTAATTTTA
TATTTATAAACATTAATAATTCTTATAAAATTAGACGGTATTGTAATATCAATAGAAGGTGAAATAGTTTTTTTAAAAAAATTATTTTGAATATTTTCAT
ACTTTAAAAAATTATCAATATAATAATTATCATTATTACAAACGATTATTTTTCCAATAGGAAAACCCATTTTTTTTGCAACGTAAACAGATAAAGCATT
TCCAAAATTTCCAGTTGGAACGTAATAATTAATTAAATTATTTTTAAATAATTTAAGTGAAGAAAAACAATAATACACTGTTTGTAATACAATTCTGAAC
CAATTTATAGAATTAACAGATACTAAAAAATTTTTTTTAAACAATTTTAAATTTTCAAAAATTTTTTTAATTAAAAATTGAGCAGAATCAAAATTTCCGA
TTAATGAAATATTATGAATATTATTTCTTTTAATAGTAGTTATTTGTTTTCTTTGAATGTCTGATATCATATTAAAAGGAAAAAAAGTAAACAATTTAAT
ATTGTTAAAGTTTTTAAAACTATTAATAGCAGCAGAACCAGTATCTCCCGAAGTTGCACAAAAAATTATTATTTTTTTATTTATTATTTTAGACACAAAA
TTTAATAAAAAACCTAATGGAACAAGGGCAATATCTTTAAAAGCTAAAGTTTTGCCAGTATTTAAATTTAATATATAATTATCATTTAGTGATTCCAATT
TTATAATTTCATTATTAAATTTTAAATAACTATTGTTTAAAATATTATAAAGTTTTTCATCATTAATATAATCTATAAAATATTTTGAAATATAAAAAGA
AAAATCAATATATGATTTATTTTTTAAAAAAGAAAAGTATTTTTTACTTATTTTTGGTATTTTAAATGGAAAAAATAAAGTTTTATCTTTTGGAATATTA
TTTAAAAATATGCTTATAAAATTTTTTATGTTTAATTTATTATTAATACTATTGTAATACATTATATAAACTTTGATATATAAAAATTACTATTCTTTAT
ATAGTATAAAAAAAACAATATTTTTTTATAAAAATTTTTTTTTGTTTTTAAAACTACTTTTTTATTTATATAAAACTTAATTATTTTAATTTTGAAAATT
TTTATTATATATAAATAATTGTAATCTATTTTTAAATACATTATTAAATTAAAGCAGGAATATATATTGCTTAATAAAAAAAAATTGTTTAAACATTTAT
AATTATTAAAATATTTTTTTTTAAAAACATTAATTAAATCAATATTAACAGTTTCTCCTGTTTCTTCTGCCCCAGCACCTGGAGCTATTAAAATATGTTT
TTTAAAATTTTTGTAATTTAATAATGTTAAATTTAAAGAATTTTTAGTTTTATATAAAAAAACATTTTTTGTTAAAAAAATAGAAATAACAGAAAATATA
TTATTGTTTATATTAACTAAAATAGTTAGGTATTTTTTTATATAAAAATTTTTTTTAAGAAAATTATTAGTTCCAAAAGTACTTTCTAAATTAAACTTAA
AATAACAAAAAAATTTTTTAAATAATAAAGAAAATAAAATAGAATGTTTAAAAAGTAAGTCTAAACCTAGTAAATCATTTGAATAATTATTTTCAGCAAA
ACCTTTTTTAATTGATAAATTGATTAGTTTTTTAAAATTTGATTTAATTAAATTAGTTAAAATAAAATTTGTAGTACCATTTAATATGTTAATACAATAA
AGAATTTTATTATTAAAATAATAATTGATAATTAATCTGATTAAAGGTAAAGATCCTCCTACTGATGCTTCATAAAAAATTTTTATATTATGCTTTTTAA
ATAAAAAATTTAAAAAAAAAGAATATTTAGAAATAAATTCTTTATTAGCAGTAACATAATTGCATTTATTTTTTATTGAATTTAAAATTATTTCAATAAC
ACAATTTACACCACCAATTAATTCTACAAAAATATTATTTTTAACAAATAATTTTTTATATGTAGAAAATGTTTTGTTTAATTTACAATTAAATTTATTT
TTTCTAGTGAATGTAATTACTTTACTTTTGTTTTTAAATATATTATAAAAACTGCTTCCAACCACACCTAAACCAAAAATAGAAATTATCATTTTATTTT
TCTTAATGAAATTAAATTTGGATTTTTAATAAAAAAAATTTTTTTAACATACATGTTTAAATATTTTTTTTTTATTGTAAAAGTATTTTTTCTAATTTTA
ACGATAATACCAATAAAAAAAAATACTTTATTGATACTTAAATACGTAATTAAAACTTTATCATTTTTTCTTATATTCATATAATTTAACATTTTTTTTT
ATAAAAAATGGAAAATAAATAAAAATTTTTTTTTTAAAAATAATAAATTTAAATAAAATAAACTTAAAACTATAATTAAAACATAAATATAAAAAACTAT
ACTTTTTAATAATATAAAAAATTTTATTGTAATAATAAATTTTTAATTTTTTCATATTTTAATAATATTAAATAATTTTTTACTAAATTTGGAACCTGTT
TTAATATAATAATTTACTATTTTAAAATTTAAAAAAAATTTTTTTCCATAATTTATGTGTGAATTGTAATATCCTACTCTTTTTAAAATTTTACCTTTTA
CAGAATTTTTTTTATAAATAATATTTATATAATAAAAATTTTTAACTTTTTTAGAAAGTCTGATTACTATCATTTAATATTTATTTTTTTAAAAAAATAT
TTCAAAATAACAAAATATTATTGCGGGGTAGAGCAGATGGTAGCTCATCGGGCTCATAATCCGAAGGTCGATGGTTCGAATCCATTCTCCGCAATTTATT
TTTTTTTAAAAAAATTAATTAATTATTTTCTTGTAATTCCAACTAAACATTAATCTTATTTTTTTTCCTACACATTCAATTAAATGATTTTTATATTTTT
CTCTAAAAATATTTAATTTTAAAAAATTATTACTTTTATCTAATAAAAATTTTTTTGCAAAATTTCCATTTTGAATATTTTTTAATATTTTTTTCATTAT
TTTTTTAAAAGTATTGAAATAATAATTTTCTGATGTTAAATAATCTCCAAATTCTGCAGTATTTGAAATTGAATATCTCATATCTTTTAATCCACCTTCA
TAAATTAAATCAACAATTAATTTTAATTCATGTAAACACTCAAAATATGCTAATTCGGGAGCATATCCTGCTTCTACTAAAGTTTCAAAACCGCATATTA
TAAGAGAGCTTACTCCTCCGCATAAAACAGATTGTTCACCAAACAAATCTGTTTCTGTTTCATCTTTAAATGTTGTTAAAATTAAACCTATTCTAGAACT
ACCAATGCATTTTGCATAAGCTAATGCTATGTTTAAAGAATTTTTATCATAGTCTTGATATACTGCTATTAGCGCTGGAACACCATTTTTATCTAAAAAT
GTAGATCTAACAGTATGACCTGGAGCTTTTGGTGCTACTAAAATAATATTAAATTGATTTTTAGGAATTATTTGTTTAAAATGAATATTAAAACCATGTG
AAAAAATTAATGTTTTCCCTTTAATAAATTTATTTTCAATATAATTATAATAAAATTCTTTTTGGTCTTCATCTGGAATTAAAATAACTAAGATATCAGA
TAATTTAACTGCATCTTCAGGTTTTAAAATTTCAAAGCCTGCTTTTTTCGCTTTATTAATTGATTGAGAATCTTTTTTTAAACAAATATAAGTTTTAATA
TTTGAATCTTTTAAATTTAAAGCTTGAGCGTGTCCTTGCGATCCATAACCTATAACTGTTATAATTTTATTTTTAATAAAATCAATATTACAATCTTTAT
CGTAAAACATTTAATTTTAATTTTTATTTATAAAAACAAAAAATTATTCAATATACTTGATAAGTATAATTTTTTTTTAAAACTTGAAAATATAAGATAC
TTTTTTTTTTTAAAAATAAAAGTAAAATTTTTATAAAATTTAAATTTTTTAAAAAATTTTAATTTGTTAAAATATATAAAATTATAATTTTTACAAAAAT
AAATTTTATTAATATTAAAACATTTTTTTGAAGATTTAATAAATTGATTAATTTTATAATATTTATCATTAATAAATATTATTTTTTTATGAAATATATT
GGAAAAAAAAAAAATTTTTAAAAATATAATTTTTAAATTGATTTGTATTAAATAATTAAAAAAATTGCTTGTATAATTTTTATTTTTAAATAAATAAATT
AATTTTATTTTCATATTTTATTTTTGAAAAAAAAGCAATAAATTCATTTACTGATTTATTATTAATTTGTATTGGCAAAACATTTTCATTTATATTTATG
TATATATTAATAAAACAAAATATATCTTTTTTAAATATGCTTATAAAAAAATAATAAAATAGTTTAATATTAAAAAAATTAAAACAATGAAACAAATAAG
AAAACATTAATTTTTTAAAATTAGGAATAGAATTCATATATGAACTTGAATACCTTTTATAATAATTTAGTTCTTGCCATTGTTTTACCATTCCTAAACT
TTGATTATTTAAATTTAATATTTTAATGTTTATGCTATATTGTTTACAAGTAGATAATTCTTGCATCATCATTTGAAAACTACTTTCGCTTGTAATTAAT
AAAATACTATTATTTTTGTTTGCAAATTTTATTCCAATACTTGAAGGTAAACCAAACCCCATAGTTCCTAATCCACTTGAAGTTAAAAATCTTTTATAAT
TATAATAATAGTATTTTGCAGCAAACATTTGATGTTGTCCAACATCTGTAGAAATAAAGTATTTTCCTTTGCTAACAAAATATACTAGTTCAATTATTTG
TTGAGGTAAACAAATTTTAAAATTATTTTTATAATTATTTAAAAAAAATTTTTTATAGAATAAAATTATTTTCCACCATTTTTTATTTAATTTTAAACTA
AATTCTTTAATAAAATTAATATTTTTAAAAAAAAAATTACAAAATATTGTTTTAGATATAGAATTATCATTTATATCAAAATGTAAAATTTTAGCATAAG
GCGCAAATAATCTTGGATTATTTGTAATTCTATCGTCTAATCTTGATCCTAAACAAATAATTAAATCTGCAAAATTTAAAATAATATTTGAAACAGAATT
TCCATGCATTCCGATCCAACCTAAATAATTTAATTTTCTATAGTTTATTTTGCCAATTCCCATTATTGTACTTACAATTGGAATTTTTGTATTTAAAATT
AAATTTTCTATATTTTTATTTTTATAATTTTTAAAAAAACCACCTCCAATTATTATTATTGGTCTTTTATAATTATTTATATTAAATGTTTTTTCTTTTA
TATTTTTTTTTTTTTTTAAAAAAAAAGGATATTTATATCTAGATGTATTTTTAAAATACGTTAAATCTTTTGGAAAATTAATAATTATTGGTCCATTTAT
AAAATTTAAAGAATTAGAAAAAGATTCTTTTATTAAATATGGAATATTATAATAGCAATTTAAAGAATAAAATTGTTTAACAATAGGTAACGATATACTT
AAATTATTTAATTCTTGAAAAGAATTTTGTGCTATTAAAATTTTATTAACCTGACCACATAATATTATTATTGAAGAATAATCAATACTAGCAGTTGCTA
TACCAGTTACACAATTAGTATATCCTGGGCCAGAGGTAACAATAACTACTCCAATTTTATTTGATGATCTAGAATATCCATCTGCCATATGAATCGCTGA
TTGTTCGTGTCTTACTAATATATGTTTAATTCCAGATTTAAAAATTTCATCGTAAATATTTAATACAGCTCCACCAGGATAACCAAAAATAAATTCTACA
TTTTCGTTTATTAAAGAATTAATTATTATTTTTGAAGATGAAAACTTCATTTTAAAAAATTAAACAAAAAAATTTTTTTAATATCATTATTTAAACTACA
AAATAAATTTATTTTTTTTTTAATAAAAGGATAGTAAAACATCATTTTTTTAAAATGCAAAGTATTGTTTAATAAGTAAGTTTTATTATAATAAAAATCG
TTAAATATTTTCTTTTCAATGTGCATCATATGCAATCTAATTTGATGTGTTCTACCTGTTTTAATACATATGTTTAAAAAAGAAAAATTTTTTTTTTCTT
TTATTATTTTGTAATAAGTAATGCTTTTTTTATTATTTTTAATTATACTACTTTTTTTTTTTTTAAAAATATCAACATTAACGAAACCTTTTTTTATTTT
TTTATTAACAAAAGATATGTAATTTTTTTTTATTAAATTTTTAACTATTAATTTTTTAAAATAGAAATAAAAATAAAAGTTTTTAGCAACTAAAATTAAA
CCAGAAGTATATTTGTCAAGTCTATTTAATATTCCAAAATTGGGAATTTTATCATTTAATTTTATAAAATTTACATTTAATAAATTTTTAACTTCAATAC
CTATTGGTTTATTTATTACAATAATATTATTATCAAAATAATTTATTTTAAAATTAAACATAGAATCTTAATTCTTCTAAATTGTGTAATTTATTTAAAA
CACAATTTTTACAAATATTTTTATTTAAAACATAATTTTTTTTTTTAACAAAATTTTTTGAAAAAAAACAATAAAAAACTTTCAAATTGTAATAAAAGTA
TTTTTTAATACAAAACCAATTCCAAAATTTGTTGTAAATATTTAAATAAATATCTAAATTTAATAATTCATTTTTAACTTTATTAAAATCTTTTAAAAAA
TTTGAATAAAAAATTGTTTTTTTATTTTTTTTTTTTACAATTTTTAAAAAATGTGAAAATTTCCCTACTTTTTTTTTAAATAGTATTTTTTTATCTTTAA
CAAAATTAATAATTTCGTTAAAAGTAAAAGATAAAACTGGAAAAATTGACTTTTTAATTAACAATAAAATATTATACAATGCAAATAAAACACTATTTCT
TATTAAAGAATTTAACTTGCTTAGATATAATTTATTTTTCGAATAATTAAAATAACTATTACTTAAAAAATTAATTATATAAATTATAATTTTCATAGAT
TTATAAAATTTAAATCTTTTAAAATTATTATCAATTATAAAAATAAATTTTTCAATTTTTTTAATTAACCAAAAATCAAATAATAAAACTTGTTTTTTTT
TAAAATTAAAATTATAAAAATTATTTATTATAAATTTAAAAAAAATTCTTATTTTTTTATAAGAAATAGATAAATTTTTAATTTTAAAATTATCAAAATA
AATATTTTTAAAAAAATTATGTTTTGAAAAATAATATTTTATTATTTCATTACTGTATGTATTAAAAATTTTTTTTATTTTAATATAATTTTCATCTGAT
TTTGAAAATTTTTCACCTTTAATATTAACACAAAAATTATGCATTATTAAAATTTTGATATTTATTTTTAAATATATTAAACATTGAACAATAATACAAG
ATTGAAACCATCCTCTTATTTGATCTTTTCCTTCTATAATAATATTATTATTATTAAATAATAATGAAGAAATTGCAGAATCAAACCATACATCTAAAAT
ATTTTTTTTTTTTTTTTTTAAAAAAAAAAAAACGTGTGAACTAAATGTTTTTGTTAAATTTTTATAAAAAGAAAAATTTTCTTTTTTAATATTAGTTTTA
TAAAATATATAAACACCCCAAATTCTTTGTCTAGAAATACACCAATTTGATCTTAAATTTATCATATTTGATATAATATTTTTTGTTTTTTTTGGAAAAA
AAAAAATTTTTTGCAAATTATGTAGTAAAATTTGTTTAATACGATATTTTTTATACATTAAGTCTAAATCTATAAAAATTTGTTTACTTAAAAAATAAAT
AATTAAATTTTTATGTCTCCAACAAAACATATAATTGTGTATTATTTTAAATTTTTTAATTAAAATATTATGCTTGATCATTAATTTAAACATTATAATG
TTTAAATAAAAAATATTTAAATTTAAAAATAACTTTATTTTTCTTAATAAGCATTTACTTGATATTAAATTATATATTTTACGATTTTTATCGTATAATT
CGTAGTCTTCTATTCCATTTGAAGGAGCAGAATGAATAAAACCTGTTCCTATTTTTATATCTAAATATTTATTTTCTATAAAAAATAGTTCTTTATTATA
AAATTTTAAAATATTTTTTAAAAAAACAAAATATTTTTTTTTAAAAAAAAAAGTACTATTAATATATGATATTATTTTTCCTTTTATATTTATTATTTTT
AATACTTTTTTGTAAATACTATTTAAAAATACTAATAAACAATTTTTAGATTTTAAAATTATATAAATATTGTTTTTAACAAAAAAACAAGCTTGATTAT
TAATAAAAGACCATAAAGAAGTTGTCCATATTAATAAATACAAATTATTTATTTTTATTTTAAAATAAAAAGAAAAAGATTTTTTTTTTTTATATATTAT
TTCTGAAAAAGATAAAGTTGAATTACAATTAAAACATAAAAAATTAGGATAATTTTTTATTTTAATTAATTCTTTAATTAAGAAATATTGAAAAATTTTT
AATTGAAACGATTCATATCGTGGTTCCATTGTATTGTAATTATAATAACTATTAAACAAATTTAGTTTTAAAAAATCATTTTTTTGATCTTTTATAGTAT
TACATACAAAAGATCTAAAATTATTAAATTTTTTTGTAGACAATACTTTTTGTTCAATAGGTAAACCGTGGGTGTCCCACCCTAGATTATTTAATAACAA
AAAACTATTATTAAAATAATTATTTATTATAATATTCTTAATTACTTTATTTAATATATGACCTAAATGAATTTCTCCATTTGCAAATGGAGGGCCATCT
AATAAAAAATAATTTATTAATTTCTTTCTTTTTAAAATATTATAATTTATATTGTTTACAATAATAAGTTGGTTTTTAAGAAGATTTCCTTGCATAAAAT
ATTTTTTTTGTTGAAAATTAATTTTCATATTTTATTAACTTTGATTTGTATTTAGATTTTTTTTTTTTACTAAATTTTGTTTTATCAATAATTGAATAAA
TATTAAAATTACTTTTTTTTTTTTTATATTTTCTTATAATTTTTTTTATTATGAATTTTTCCATAAAAATAAATAATGAATAATTTTTTTTAAACCGAAA
TTAAAAATTTTAGAAATAAAAAATGTTGGAAAAACAATTAATTTAATTTTATTAAATAATTTTAAAAAAAAAATGTTTTTATTAATATTTAAAATTATCC
AATAATAATTATTTTTAAAGTTTTTTAAGTACAAAAAATAAAAATAAAATTTAATTATTTTAATTAATATATGAATTATATTTAAATTTATAATAATTAA
ATAAATTTTAGATATAAAATTAATTTTGAAAACATCTAAATTATTGCTATTTATATTTAAAACAAAACTTTTAATAAAGAATTTATAAACGAACAAACAA
GTAACAAATTTTGTTTTTTCTCCTTTAGTAGCTATTTTTTTATTAAAATAATTTTTAAAATAATGATTTCCTAGTCCTCCTTTTCCACCTTTTAAAATTT
TAACAAAAAAATTATTACAAACTAAATAAATTTTATAATTTTCAATTTCAATATAACTACCTATAGGAAATTTAATAATTAGACTTTGTCCGGTTTTTCC
ATGTCGAATTTTATTTTTTCCAGTTTTGCCATTTTCAGCAGAATAATTTATCTTACTAATAATTTTAGTATTATTATTACAATTAAGGTAAACATCTCCT
CCATTTCCTCCATTTCCTCCATTTGGAAATATTTTATTTCTAATTTTTAAAAAACTTATTAGACCATTTCCACCGTCTCCGCTTTTGATATAAATTATTT
GTTTAATAAACATTTTTAACGCTTACATAAGTTTTTTTTTTTTTTTTTTTAAAAAAAATTAAACCATTAATTTTTGATTGTATTGTATAATTTTTAGATA
TAAAAACATTATTTCCAACTAAATATTTAGTTCCATTTTGTTTAATTATAATAGAACCTTTTTTTACAAAATTATTATGAAAAATTTTAATTCCTAATCT
TTTTGAAGAAGAATCTCTTCCATTTTTAGTACTTCCTCCTGCTTTTTTTTGAGCCATAATTATAATATATTTTGTACAAACAATAAACTTTTTTTATTGT
AATTGAATATTTTTTTTTGAAAATTTTTTCTTCTTTTTTTTTTAAAAGAAACTTTTTTTAAAATAAAGTGTTTTAATAAAATTGCATTTATTTTAAAACA
ATTTATATTAAAAAACATTTTTTTATAAAATAATTTAATTTTATTAAAAATTATATTTTTACCTGGAATAAAATTTATATAATCAACAATTAGATAATTA
TTTATTTTTGCTAAATAATATTTATTTCCAATTATAATTACAATACTTTTATTTTGTTTTTTTATCATTTTATAACTTTATTTATATAATAATTATATAC
TCCAATTTCACCAATGAGTACATTATCACCTTTTCTATTAAATAATTTTATTAAATAACATTTTTTTTTTATTTTTGTAAATAATTTAATTATTTTTTTT
TTAATTTTTTTAAGTCTAAAATAATTACTTTTAATTTTCCCATATTTCAGCATATTTAAATAAATATTACTTAAATACGTTTTATAAAAATTATAAGTCA
TATTCTATTTTTGAATTTAAGTTCAAACCAAAATTTTTTAACGAATTTAAAATTTCATTATATGATTTTTTGCCTAGATTAGGAATATTAATTAAATTAT
TTTTACTAATTTTAATTAAATCGCCAATTAAAAAAATATTATTATTTTTTAAACAATTTAAACTTCTAATACTTAATTTTAAATTATTTATTGAATTTAA
ATAAACGGAATTAATTTTTAAATTTAAATTATTTTTTTTTTCAACATTTATTTTTTTATACTTTTTAAAACCAATGAAACTAAATATTAAATCAAAATAT
TTTTTAATATAAAAAATACAATTTTTAAAACAATCAACGGGTTTTATTGAACCATTTGTTTCAATATCAAAAAACAATTTTTTCAACTTTTTATTAAATA
TTTTTTTGTTTATATAATAATTTAAACTAATTATTGGTGATTTAAAATTATTTAAAAAAATAATATTTTCTTTAAAAATTTTAAATTGAAAAAATTCATC
AGTATAATTTTTAAACAAGGAATTAACGCATTTCATTATGCAATAGAAAACTATATTATTACTAACATTAGCAATAATTTTATTAGGATTAAATATTGTA
ATATTTTTATCTGAAAAAATATCTTTTGCTGTAATTATACAAGGGCCTTTTTTTTTTATAATTAAATTAGCAACATTGTCATTATTAATTTTAATAATAA
TATTATCTAAATTTTTTAAAATTGTTTGTGTATTTTCTATAATACCAGGTAAATCATAAAATTCAGACTTTATTTTATAAATTTTTAAATATATTATTTT
ATAAGAATTAGTTGTTAAAAAGATAACTCTTTTAATAAAATTACCTAACGAATCACAAAAAGAATTATTAAAAGTTTCAATTCTAATTATTGATCTAAAA
AAAGAAATATTTTTAATGCTTATTTTTTTTAAAGTAAATAAATTATTTACAAAATTCATTTATAAATATCTAAAATTAAGTTTTTGTTTAGATTAAATAA
GTTTTTGTTTGAAAATGATAAAAAAATTCCATATCCATCAAAATATTTATTATATAGCCAGTTAACAAATAATATATTTTTATAATCATATAAATATTTT
TTTAATAAAAATTTTTTGTTTATAGTAATTACATCGCCTGGTGATAATAAATAGGAAGGAATTTTATTTTTAAAATAATTAACAAAAACATTACCATGAA
TTATTAATTGTCTTGATTCTTTTCTAGAAATTGAAAAATTAAATCTGTATATTACATTATCTAATCTTCTTTCTAAAAAATTAATTAAATCGATATTATT
AAAAAAATTTTTTTTAATAATAAAATAAATTTTTTTAAATTCTTTTTCAAAGATACAATAATATCTTTTTATTTTTTGTTTTGTTCTTAATATTAAACCA
AAATCAGAAACTTTTCCTTTTTCGGTTCCATTTTCTCCTGGTGTAAGATTAGACCTTGATTTTTCTAAATATTTTTTTTCAGAAAAAAATTCTAAATTTT
CACCTTCTCTTCTACAAAACTTAATATTTTTTGATTTTTTTTTTGTCATACTCTCCTTTTTTTTTTTTTTCTACATCCGTTATGTGGAATTGGAGTTATG
TCAGTAATACTCCGAATAAAAATTCCACTATTATTAATTATTCTCAAAGTTGTTTCTTTTCCTAATCCTAATCCTTTTATATAAATATCTACTACTTTAA
TATTATTTTTTAATATAGAAATTGTAATTTTTTCACTAATTAATTGTGATGCTAAAGATGTACTTTTTTTTGGTCCTTTGTATCCTAGAATACCTGCTGA
ATACCATACAATGGTATTACCAAATACATCAGTTAAAGTAGAAATAGTATTATTAAATGTTGAATGAATATATATAATTCCACGTTTTTTATATAATTTC
ATGATTCATTTTTTTTCTAGTCTTTGCATTTGTTTTGGTTCTTTGTCCTCTACACGGTAATTTTTTTTTATGTCTTAAACTTTTATAACTATTTAAATTT
AATTTTTTTTTGAAATTCTCTTTTATTATAGTTTTTAAATTATTTTCAACATTAATTTTGTTAATAAAATTTTCAATTTTAGCTTTTTCTTCATTAGTTA
AATCTTTAAATTTTTTATTTTTAAAATTTTCAATTTTACTGCATATTTTTAAAGACATAGAA
>CONTIG2 (100bp)
ATATAAAATTATTAAAAACTTAAGTGGAGACATTAGAATTATGTTATCTCAGTATAGAGAATAAGAAGCATTTTCTAAATTTTCATCCGATCTTGATAGT

>CONTIG3 (93bp)
TAAAAAAAATTATTTTAAGTAAATCTTGTTTTTATAATCTTCATGATAATAATAATAATTATTCATTATAAAAAACTTTTAAACAAATAATAA
>CONTIG4 (92bp)
CCTGCTTTTTTCGCTTTATTAATTGATTGAGAATCTTTTTTTAGACAAATATAAGTTTTAATATTTGAATCTTTTAAATTTCAAGCTTGAGC
>CONTIG5 (90bp)
AAAAGTATTTTATGATTAAATTTTTTGAATATAAAATAATAAAATTTAGATTTATAAATTTCTGGAATAATAATTGGAAGATGGGATTTT
>CONTIG6 (90bp)
AAAATTATTCAATATACTTGATAAGTATAATTCTTTTTTAAAACTTCAAAATATAAGATACTTTTTTTTTTTAAAAATAAAAGTAAAATT
>CONTIG7 (87bp)
TATAAAAAAGAATATTTATTAGATTTTACTTTTGGTTGTGGTGAATATATTAAAAATTTAATATTAAAACAAAAATTTAAAATTATT
>CONTIG8 (87bp)
ATTGTATCACAATTAATGTTTTTAGATTCTGAAAATTCTAAAGTTATAATATTATATATTAATTCTCCTGGTGGAGTTGTATCTTCT
>CONTIG9 (87bp)
TCTTTTTAATTCAATTCTAATTGAAGAGTAAAATTTTACTGCGCTTCCTCCTGTTGAAATTTCTTTATTATAAAAATTACTAATTTT
>CONTIG10 (87bp)
ACATCTTGTATAACAACAAGAAAAATAAAAAAAAAAGGTTCAGGATTAATTTCAAATTATATTTGTGATTTAATTGAATTTAATAAT
>CONTIG11 (87bp)
AAAACTTAGAAACACTATTTCTTTTAATTTTAAAGAAAATAGTCTTTTTGTTAAAAAATTTAAATTATTTATTTTTCAAAAAAGTAT
>CONTIG12 (87bp)
CAATAACATAAAATTAAACGTTAATATCTTTTTGAAATACAATGTATTAAAAATTATTAAAATTAAAAATATAATTACAAATGTAAA
>CONTIG13 (87bp)
GCGGATTCAACTACATCATATACTACTTGTTCAATTAAAGCAGCTGAATCTGCTATAAGAGTATTAATTCCTCCTGTTTCTGCTACT
>CONTIG14 (87bp)
AAATACTTTTCAATAGCTAGAGGTATAGATAAGTTAAATATCACTAAAATGACAAAATTTTTTAATACTAACTATCATTATATTGTA
>CONTIG15 (87bp)
TATACTTCACAAATAAAATTAGTATGAGGTTTAATTGTGCCTGCTTTAATTAAAACTTGTCCTCTTTCAACTTCTTCTCGTTTTATA
>CONTIG16 (87bp)
TCCATGACATACAATATTTGAAATTGATAAACAAGAACAAAATGTATAATTTTTAAAATTTATTGTGGAAGATTTTAAGTTTAATAT
>CONTIG17 (86bp)
TTATTATTAAATAGTAATTTTTTAATTTTAGTTAAATAATAACTAAAAGTTTCATAAAAAAAGTTATTTATGTTTTTACAAGTATA
>CONTIG18 (86bp)
AAATAATTAAAAAAGAAAAACCAAAAATTTTAATTTTAGGATTTTCTTCATATCAAAAATATATAGATTGGGATTTTTTTTATTAT
>CONTIG19 (86bp)
TTGTGATTATAAATAATTTTTTTTTTTAATAAAAAATTGTTTCAATTAATTTTAAAATCATTATAAAAACTAAAATGTAAATTAAA
>CONTIG20 (86bp)
GATAGACAAACTGGAAAAACTACAATTTGTATTGATACTATTGTTAATCAAAAAAATAAAAATATTATTTGTGTTTATGTTTGTAT
>CONTIG21 (86bp)
TTTAAAGTTTTAATTTTAAAATTTAAAATTTTTTTGTTAAATCTTGAATATTTAATTACATTTAACAAATAAAACTTATTATTAAC
>CONTIG22 (86bp)
TTGTACTAGAATTTTAATAAAGTATAATAAAATTATAAAAGTAGTTTATATTGGTTCAATGTTTAGATATGAAAAACCTCAATTAA
>CONTIG23 (86bp)
ATTTAATAATATAACAGAATTATCATTAAATTTAATAATAGTTTCATCTAATCGTTTTATTCCAGTTTTGCTTCTAATTATCATTG
>CONTIG24 (86bp)
ATAGATTTTTAAATATTTATAAAAATCTATTATTATTAGGAAGAAGAATTATTAAAAAAAAAAATTATTTAATATTTAAATGTAAA
>CONTIG25 (86bp)
TAGCTATTGTTTCTTATTAGATTATTATTATGTTAAAAATATATTACATAAATCGTTTAACAAAGTTTTTGTAACAAGAAATTTTT
>CONTIG26 (85bp)
TAAAAAAACTGGTATTAGTAAAAATTTTATATTAGGAATAATTCCAAAAAAAAATATAAAAGATATTTTAAATTATACTAATTTA
>CONTIG27 (85bp)
TTTCAGTAATATTAAAATATTATATTATTTTTATTACTATACATTTTGCTTGATAAAAAATTTTGAAAATTATAATGTTTTAAAA
>CONTIG28 (85bp)
TTATATAAACATTTTTTGATTTTAAAAAAAATAAAAATTTTGTTTTTAAAACAAAAATGCTATAAAAACTATTAAACCTTAAGTT
>CONTIG29 (85bp)
CAAATGATACTCTTTGTCCATCTTGTAATGATTTAAAACCATCGACTCTAATTTCTGAAAAATGAACAAATAAATCATCTCCTCC
>CONTIG30 (85bp)
TTTTATTAAAAAAATTAAAAATTTAAATTTTTTTCTAAACAATAAAAATATGTTTCATAAGATATCAAAAATAATTATTTCAATG
>CONTIG31 (85bp)
TTTAGTTTTTAAATTTTACAAAATTAATAATGAATAAAAAATTATTAGTAATAGGAGCAGGTCCTATTTTAGTTGGTCAAGCTTG
>CONTIG32 (84bp)
TAATTCTAATAATTAATGCTTCAAATATTGAATTAATATTTATCTAATTTTTAATGGAACTAATTACGCTCCAACTTCCTGTTC
>CONTIG33 (84bp)
TCCTTGTAATAACATTAAAGACTTATCTGTTTTAATAAATTATCTAAACAAAAATAAACCTAGTTTTGTTTCTGTAACATTTGG
>CONTIG34 (84bp)
AACCTTGTTTTTTTAAATTTTTTAATATATATAATTTATTCATGTTTTTATTTTTAATTTTATTTGTATCATCAAATCTTAAAT
>CONTIG35 (83bp)
TTTTGGAATTAAAATAATTATTGCTGAAAGTTTTAGTGACATTTTTTATGATAATTCTTTTAAAAATAATTTATTTTTAATAA
>CONTIG36 (83bp)
GTTTTCCAAATCTTTTAAAATAAATCAAACTGGAGAAATGGAATTTTTAGTTTTAATTTTTAAATTAATAAATAAGTTTTTTT
>CONTIG37 (83bp)
TTATTATTATTTTTTTAAGAATAATTATTTTTTTTAATGCTAACTTTGATTTTAAACCATTTAATTTAATAATAATTTTAAAA
>CONTIG38 (83bp)
TAATTTTTTAATTTTACAAATTAATTGTAATTTATCAAATAAAAATAATTCATATTGTTGGTACTTTGGTAGTAATATTTATA
>CONTIG39 (83bp)
TCAAATGTAATTAGATCAATAATTATATTATTATCAAAGATAACTTATCTTGTTGTTTTAATATCACCTATTAAATATTGGTT
>CONTIG40 (82bp)
AATATTTTAAATAATTTTCTTAATTTTAAAATAATTGATTTAACTTTAATAATATTATTATTATTTATACATTTAATTGTAT
>CONTIG41 (82bp)
TTTACTCTATCTCCTATCATAATTTTTACATAATTTTTTCCTATTTTACCAGAAATTGACGCAGTAATAATGCTTCCGTTAT
>CONTIG42 (82bp)
TTTTTTTTTATTACAAAATTTAAAATAACAGGTATTGTTAACATAATTAAATTTTTTAAAATTAATATTTTTTATCCTTTTA
>CONTIG43 (82bp)
AAGTTTTATAAAAATATTTTTATTGTTACAAATAAGTGCAAAATATTGCTTAAGTACTTTTTTAATTGAAAATAAGGAATTA
>CONTIG44 (81bp)
TATAATTTTATTTGTTATATTTAGTTTAATTATCATTTTTTTTTTTTAGAAGATTTTTTAGTAAAACTTTTATTTAATAAA
>CONTIG45 (81bp)
TTAAAAAAGAAAATTTCTTTATTTTATACGGTATAATTATAAAAATTAAAAAACTAGGTAAAATAATATTTATAGAATTTA
>CONTIG46 (81bp)
GTAATACTCCGAATAAAAATTCCACTATTATTAATTATTCTCTAAGTTGTTTCTTTTCCTAATCCTAATCCTTTTATATAA
>CONTIG47 (81bp)
TAATGTTAAATTTAAAGAATTTTTAGTTTTATATAAAAAAACAGTTTTTGTTAAAAAAATAGAAATAACAGAAAATATATT
>CONTIG48 (81bp)
AATATGTATATTAGATAATTAGTTAAAAAACTAAAGTATAAAACAATTTTTTTTTCATTAATAAAAATTACTTTTTTATAT
>CONTIG49 (81bp)
TTGAATCATTAGTTGAAGATTTAATTTTAAAATCTTTAAAACCTTGTGAAATTGCGTTAAATGATGCTAAAATAAGTAAAA
>CONTIG50 (80bp)
AAATATTTTTATTCCAGTTTTAGTTTTACAATCAATCTTATAAATTACTAAATCATGAGATACATCAATTAATTTTCTTG
>CONTIG51 (80bp)
ATATATCTACAATATCATTTATATTATTTTTATTAAAACCAAAAAAACAATTTATTTTATATTTATTAACATATTTAGAT
>CONTIG52 (80bp)
AAATAATATTTTCTTTAAAAATTTTAAATTGAAAAACTTCATCAGTATAATTTTTAAACAAGGAATTAACGCATTTCATT
>CONTIG53 (80bp)
TTTAAAAAAAAAAATAAATCAAATTAAAATATTTTGTAATAGTTTATGTTTTAGAAAAGAATCATATAATTTACAAAATA
>CONTIG54 (80bp)
AAATTACCAATATTCAATATGAAGTTATAGTAAAGGCTCACGGGGTCTTTCCGTCTAGCCGTAGATACACTGCATCTTCA
>CONTIG55 (80bp)
AAAAATTGAGATAACTCATTATTACAAAAATAATCTTTTAAACTTATAGTAATTATATCACTATTAACAATAAAATCTAA
>CONTIG56 (79bp)
AAATAATTATTTAATTTAAAACTTAAATTTAAAAAAAACAAGTATTTAATATCTTTTAATTCTATTATTAAAAATAAAA
>CONTIG57 (79bp)
TTTATTTAAAAAAAATAAATTTTAATTTTAAAAAACTTGTTAATAAGTTTAAAAAAAAAAAAATTGTATATAATATTTT
>CONTIG58 (79bp)
CAATTTCTAATATTAAATTTATGAAAAAACAACTAAAATACTTCTCTTTAGATTATTCTAATCTAGAATTTTCAACATG
>CONTIG59 (79bp)
TTTTATACCTAATTTAAACAAATTATTCATATAAAAATATCTCTTATTTTTACTGTTTTTGTTTTTTTATTTTTTTTTT
>CONTIG60 (79bp)
TTTAAATAAGATTTTTTAATAGATTCCAATACTTCACGAGATTGACCTCTACCTAATCCTATATTTTTTTTACGATTAC
>CONTIG61 (79bp)
AAATAATATGTTTTAATATATTTTTTATTATTTTTAAATCTATTTATGTTTTCTTTTACAAATAAATTATTTTTTAACA
>CONTIG62 (79bp)
ATAATAAATAAAAATGTTAAAACAACTTCATTTATTGGAAATACATTCCAATTAAAATTGTTCATTTTTTTTTATAATT
>CONTIG63 (79bp)
TTTTTTTTTTTTTTAAATAAATAATTACTAAATGTAACTCCTCTACTTCTCCAAATTGGAGATTTTTTATCTCCAGCTC
>CONTIG64 (79bp)
TCAATTTTTTTAATTAACCAAAAATCAAATAATAAAACTTCTTTTTTTTTAAAATTAAAATTATAAAAATTATTTATTA
>CONTIG65 (79bp)
AATTTTAGACCTGGACATGCAGACTATACTTATTTTTTAAAGTTTAAATTTAGAGATTATAGAGGTGGAGGTAGATCTA
>CONTIG66 (78bp)
AAATTATTATTTTAAAACATTATTGTAAAAATAAATTATTATTCCATTGTTTTAAACATTTTTTTGAAAAACATCATA
>CONTIG67 (78bp)
CTTGTGCGAGTCCCCGTCAATTCATTTGAGTTTTCACCTTGCGGTCGTAATACCCAGGCGGTCAACTTAACGCGTTAG
>CONTIG68 (78bp)
TTTTTTTTTTTATAAAAAAAATTTTTTTATAAATTTTATTAAATTTTAAATTTTTAGCTTTATTTCTTTTTGACAATT
>CONTIG69 (78bp)
AAAAATTTTTACTAAAAATAAAAATTTTAAAGAATAATATTATATACAGAAAAACATTTAATATATATTTTAAAAATT
>CONTIG70 (78bp)
CTTATTAAACTTTAAAAAAGAAAAAATTTTTAAACAAATAAAAAAGTATCTTATTTCAAAAAAAATAAAATCTATTGA
>CONTIG71 (78bp)
GTTATAAGTTCAAATCTTATATTACCTATAATAAATGGAGCGATAATTTAGTAGGTTAAAATGTTGGCTTGTCACGTC
>CONTIG72 (78bp)
ACTCTTTTTAAAATTTTACCTTTTACAGAATTTTTTTTATAAACAATATTTATATAATAAAAATTTTTAACTTTTTTA
>CONTIG73 (78bp)
TATCATTGTTAAATTTTGCAATTCCTGATCTTATTATAATTTTATCAATATTTGGATTAAAAATATTTTTTTTTTTTG
>CONTIG74 (78bp)
AATTATATTTTTTTATTATTTTAATAATTGAAGAATTATAAGTGGAAAAATAATTTTTTAATATTTTAGAATTACAAC
>CONTIG75 (78bp)
CTTTTTTATTAAATTTGGTAAAAAATTTACTATTTTTGAATAACAAAATTATAATATTTTCATTAGAAATGACTGTTT
>CONTIG76 (78bp)
TTTAAGAGTTTAGCTTTTTTTTTTGCTAGAAGAATTTCTTCTCTTTGTCCTAAATCAGCAGTAAAAGTAATTACTTCA
>CONTIG77 (78bp)
AAGTGTTATTAATAATAAAGTTTTTATTAAAACACGTTTTTTAAAAAAAAATATTAAAAACCTTGGAAGTTTTAAACT
>CONTIG78 (78bp)
GATAGTGTTTTTTTCTTTTTTAAATGTTTTTATAAAAAAAAAACAAACATATTTATTAATTATAATTTAATTTATTTT
>CONTIG79 (78bp)
CGCATCATTATGCGATAAAAAAGGAATCAAAGACCCACCAACAGAAATAATCTGATCTCCACAAATTTCAATTAAATC
>CONTIG80 (78bp)
TAAACTATTAGAAATTATTAACATTTTAATAGGTCTTAATAATGGAATACAAATTTCGTTTAAAAGTAAAATTTTTTG
>CONTIG81 (78bp)
TTTTAGAAATTTGAAAAACGAAAACATCTAAATAATTTTCAATCATATTAAATTGATCAACATTGTTAATATATAAAT
>CONTIG82 (77bp)
ATTAAAAAGTTGTTTTTTAAAAAAAATAAATTAGTATAATTGACTTTTTTATAATTTTTAAAGAACAAGTTTAAATT
>CONTIG83 (77bp)
ACATTTTATATATTAGTATAAGTATAACAAAAAAAATGAGTCGTATAGGATTTGAACCTACGACCAATTGGTTAAAA
>CONTIG84 (77bp)
CTTAGATGTTTCAGTTCACCCGGTTCACTTTTATTTTTAAACAAATACTAGTTTTTAAGTAGTAGGTTACCCCATTA
>CONTIG85 (77bp)
ACCGCATTTAAGAATAATTTTAAAAATTTTTATAATTCTAAATACTTAGCATACGATGTATTTTTAAATGGAATAGA
>CONTIG86 (77bp)
TTATTTTTTTTAAAATTATAAAATACTTTAAAACAAAATTTAAATTATTAATATTAACGTTTAGTTTGTATGAACAT
>CONTIG87 (76bp)
TTCAATTTTTAATTTTAATAAAAAAAATAAAATTTTTTTTTTAAAAAAAAAATTAAATTTTTTTTTTAAAATTAAT
>CONTIG88 (76bp)
AAACAATAGTAAATTTTTATAATTTGAAAACATAAATAAATCGTTTTTATTAGATTTAGACATTTTTTTTTTGTTA
>CONTIG89 (76bp)
AAAAAAATTCTTCTGTAATTTGTTATTATGAAAAAGTATTTTCTTATATTTTATTATTATAAAAAAATAATCTATT
>CONTIG90 (76bp)
TCAACATTGTTAATATATAAATTTATTGTTAACATATTTGATATGATAATATCAATATAACTTAAATTAATAAATA
>CONTIG91 (76bp)
AAATAGTTAAATCATGTATTTATACTGGTTATAAAATAAATATCATTTACAATTATATTGGTGGTATTTCAATCGA
>CONTIG92 (76bp)
TCAAATGGTCTACATTCACTCCCATATAAATCTGAAGAAACTTGAGTTAAAGCAGCAGTCAAAGTAGTTTTTCCGT
>CONTIG93 (76bp)
TAATAAAACAAATGATAATATTAATGTATAAAAATAACGGTATTGGAATTTCTAGTAATCAAATAAATTGTTTTAA
>CONTIG94 (76bp)
ATAATATTTTTTTTGGAAAAAAAAAAAAATTTTATAAAGTATTATACATTTTAAAATTAAATAGAATTAAAAATTA
>CONTIG95 (76bp)
ATGATTAAGAAAATTCTATTTGTTTCTTTATACTATAAAATTTCATTAATTATAGGAACAACAGGTTTTAATTTTA
>CONTIG96 (76bp)
CTAAGCTAGATAATTCTAGCATATTTTCAAATGGAAGAAAATTTATTCTATAATTTAAACAAATTCTAGGATCAGA
>CONTIG97 (75bp)
TATGTTTAATTTATTATTAATACTATTGTAATACATTATATATACTTTGATATATAAAAATTACTATTCTTTATA
>CONTIG98 (75bp)
ATTTTTTTCATTATTTTTTTAAAAGTATTGATATAATAATTTTCTGATGTTAAATAATCTCCAAATTCTGCAGTA
>CONTIG99 (75bp)
TGTAAAATTAATTTTTTATTAATAATACATTTTATTAAAAAAACAATTTTGTTTAATATTTTATTTTCTAAAAAC
>CONTIG100 (75bp)
TTATGAAAAAATTTTTTTTTAGAAAAACAAAATCTATTAAATTCTAAATAATTTATCATTTAATTTCTACTATAC
>CONTIG101 (75bp)
TGTCATACTAGGAATATTACCTAGAAATATTTTTTTATAAACACTTAAATTTTTATTTATTACGTATATACTTAT
>CONTIG102 (75bp)
TTTAATATCATTATTTAAACTACAAAATAAACTTATTTTTTTTTTAATAAAAGGATAGTAAAACATCATTTTTTT
>CONTIG103 (75bp)
TATTTATAAATTTAGTTTTTAACTTTTTTTTTTTAATTGATTGTATTTTTAACAATAAATTGAATGTAAATTTTT
>CONTIG104 (75bp)
ATATTAGTAGAGAATTCAGAATAAATTTTTTAAAAATTATAAATCAAGAATATGAATATTTTATAGAAATTTTAG
>CONTIG105 (75bp)
TGGAATGTTCTTTTTTTTTTTCCATTTTCATATACTTACATTTGTCCAACTGAATTACTAGAAATATCTAAAAAT
>CONTIG106 (75bp)
ATTTTTTTGTTACAGGAGTTTGATAAAAATAAAATTTTATTTACAAATAATAAATTGTTTAAATTAACCATATGA
>CONTIG107 (75bp)
TTTTTGGATTTTAAAAATGATAATAATAAACATATATTCTTACAAAAAAATATTATTTAAAAAAATTTTTTTACA
>CONTIG108 (75bp)
TAAAATAAACTTTTTAAATAGCGATAATTTATATATATTTTTGCTTTTAATTACACATTTAAATAATAAACTATT
>CONTIG109 (75bp)
GATTAGTTGGATCAGTTTTTATAAATAGAATATTTACTTCGAACATAATCAAGTATTTAGAAATATATTTATTTT
>CONTIG110 (75bp)
AATATTAATTAATAGCAATTTATTTAATATAATTAAATTTTATTAAAATTAAAAATTTTTTATTATATATTTTTA
>CONTIG111 (75bp)
AAAAATTTATCATTATGACAAAAACATTTTTCTTTTAAACAACATAATTCAAAATCGCTAAATATTTTTAAACAA
>CONTIG112 (75bp)
TTTTTTTGAAAAAAAAACTTTTGGTAAATCAATTATAACTCTTCTAACACTTATTAAACTAAACATCTATGTTTT
>CONTIG113 (75bp)
TTATTATTTTATAATATTTGTAAATGTACAAATTATTAGTTAAAATGAAACTTATAGGTGATATTAATATTGCTA
>CONTIG114 (75bp)
CTTTTTTTTTTATAATTAAATTAGCAACATTGTCATTATTAATCTTAATAATAATATTATCTAAATTTTTTAAAA
>CONTIG115 (75bp)
ACTCTTGTAGAAATGATTTGTTCATTATTATTTAATAATATAATAGAATTATCATTAAATTTAATAATAGTTCCA
>CONTIG116 (74bp)
GAACTTTAACAATTTAAATATTTTATTTTTCAAATAAAATATTTTAAAATCGTTAACTGCAATTATTGTTCCTT
>CONTIG117 (74bp)
CCAAATAATACACCTCCTATTACTACAGAATGAAAATGAGCTATTAAAAACATTGAATTTTGGAAAATAAAATC
>CONTIG118 (74bp)
ACAAAAAAATATTTTATTAAAATAATAACTTTTAAAAACAAACAACTTTATAATTATAAAAAAAGTATATTTAA
>CONTIG119 (74bp)
AACTTTAAACATAAATCATTTAAATTATTTTTTTTTCTTGGATGTAATTTTCTAAATAACAATAATGAATCTAA
>CONTIG120 (74bp)
ATAATTTTCATAGATTTATAAAATTTAAATCTTTTAAAATTATCATCAATTATAAAAATAAATTTTTCAATTTT
>CONTIG121 (74bp)
GATTTTATTTAGTAAAAATATTAAAGGAATATTATTCAGTTTAAACAAAAAAAATGTAAATATAATTATATTAA
>CONTIG122 (74bp)
GTTTTAACTATAATATACTTTATATAATAACACAAATGATAATATTAATGTATAAAAATAACGGTATAGGAATT
>CONTIG123 (74bp)
CTAGTATATTTTTTCTAGTGTTTTTAAGTAAAAAAAAACATAGTAGTTTATTATAAGATTTATTTAATATTGCT
>CONTIG124 (73bp)
CAATTTTTCCACCACTATTATCTGCAACTTCTAAATAAGTTTGTTCTTGAATCATATTATTTTAACTTTTTCA
>CONTIG125 (73bp)
TAATATTAAAACATTTTTTTGAAGATTTAATAAATTGATTATTTTTATAATATTTATCATTAATAAATATTAT
>CONTIG126 (73bp)
TTTTTATATTATAATATACCTATTTTAAGCGTATGTTTAGGACTTCAAATTATAAGTATTTTAAATAAGTTTA
>CONTIG127 (73bp)
AAATTTTTCATATCTTATAGCTAATTCTATAATTAAATTATCTGACTGTAATTCTCCTATTATGTTAAAACCA
>CONTIG128 (73bp)
CTCACGGTACTAGTTCACTATCGGTTAACTATTATTATTTAGCCTTAGAGGATGGTCCCCCTTTTTTCTGTAA
>CONTIG129 (73bp)
TTTTTTTATAAAAAAATAAAAAATAATAAATTTAAAAATTATAAAACAGTGGATACTTTTATCCAATCATCAT
>CONTIG130 (73bp)
GCACCATTTTCTAATAATCTTCTTACTAAATATGGAAGTAATCCTTTATATTTTCCAATTGGAGCATATTCTC
>CONTIG131 (73bp)
TAAATTTACTTTTTTTTTTAATACTAAATTCAAATATTTTAGTCTTACAATATTTTAAAAATTTAAAAAATCC
>CONTIG132 (73bp)
AATTAGCTTACCAAACTCTTCTATTCAAAATGCGCAATCAACTAATAAAATTGTTAATTGTATTAATTCTGGT
>CONTIG133 (73bp)
TAGGAATAATTTCAAAAAAAAATATAAAAAATATTTTAAATTATACTAATTTATGTAAAATAGATATTCCAAT
>CONTIG134 (73bp)
ATTCAATTTGTTTAAAATCACAAATTAAATATATTTTAGGAGAAATAGATGTTGTATTTAATTCTATATGGTA
>CONTIG135 (72bp)
TTAACACCAGCATCATATTTACAATAAACTTTATCACTAACTCTACCAACAAAAGGTGATACAATATATATT
>CONTIG136 (72bp)
AGAGCTAAATAGTATTAAATTAACAAATTAAGTTTTATTATTACATAATAATTACAAACTAGAATTAAATAA
>CONTIG137 (72bp)
TTTTGCTAATAGTATTAGCAACTTCATTAGAACCTGTAAAAATTATTCCACAAATTTCTTTATGAAACGAAA
>CONTIG138 (72bp)
TTACCATCAAAATATGATCCCGATTTATTCATTAATAAATTTGGTATAAAAGAATTTTTAAATTATTTAAAA
>CONTIG139 (72bp)
GATTATATTTCGGGAGAAAGATTTTCTAGAAGCGTTCCATTTTTAATAGAATTATCTATGGAACATATTTTT
>CONTIG140 (72bp)
ATAAAAAATACAATTTGATATTTGAAGACTATATAATTTGTAAAATTAAAAAAATAGAAAAAAATTTTTTTT
>CONTIG141 (72bp)
TTTATTTTAAATAAATTATATAAAATAGATAAATTAAATCCAAATGTGTTTTTTAAAAAAACAAAATCTTTT
>CONTIG142 (72bp)
AGAAGTTATAATGTTGGATTTTCTGGTGCTGAAAAAAAAAAAAACGAATTTTTATTTTTATTAATTATTAAT
>CONTIG143 (72bp)
AGAATTAAAATTAATTAAAAATATTTATTTAAAAAAACTTAATTTTTTAGATATTTTATATTTAAAAAATTA
>CONTIG144 (72bp)
ATTCTAATATTTTAAAAACAAATATTTTAGAAATTTTAAAAACTAAATATAAAAAAAATATTTTTTTTTCTT
>CONTIG145 (72bp)
AAAACAAAAAAACTTTTAATTATTTTTTTTTCAATTAAAATTCGAAGTAATAACATGCTAATTTTCGAATAA
>CONTIG146 (72bp)
TTATTTTTTTTTTATTATTAAATTTAATATTTATTTGAAAATTTTAATTATTAAATTCAAAAAAGTTTATAA
>CONTIG147 (72bp)
ATTTAAAAAAAAAAACTTTAAAAACACTAAAATTTCAAAAAGATAATTAAAATTTAATAAAATTTATAATTT
>CONTIG148 (72bp)
ATTGATGGATATCCTATCAAATTTGAAAATACTGCAATATAATCACAATATTCTGAATAATTTATATTATCT
>CONTIG149 (72bp)
TATTAATTTTAAAAATGGAAATTAAAAAAATAAAAATGGAAAACATATATTATCAAGTAAATAATAATTTTA
>CONTIG150 (72bp)
ATAATTTAGATATTGTTAATATTAAAAAACACTTTTCAATAGCTAGAGGTATAGATAAGTTAAATATCAGTA
>CONTIG151 (72bp)
TAAAATTGAAATAATAACACTATTAACAAAAAAAAACAAAATAAAATTATTAAATTACAAACATTTAAATAA
>CONTIG152 (72bp)
TAAAAAAAAATAGTTTTTTTTTTTATTTACATAATTATTTTTTTTTTAAAGACTTATTTGCTTCAATAATAA
>CONTIG153 (72bp)
ATCCCGTTTTCGGTATTTTTAAATTTTTATATAATGGATATTATTAAAAACATATATATAGATAATAAAGTT
>CONTIG154 (72bp)
AATATATTTTTAATAAATAAAAAAAATGATATTAATTTAACAATAAAATATTTTATTAAAATAATAAATTTT
>CONTIG155 (71bp)
TTACAAAAAAAAATTTATTTGTTCGAAATTTTATCTAAAATAATTATTTTTAACAAAAGTAAATTTATAAC
>CONTIG156 (71bp)
TTTTAAAGTTTTATAGTAAAATTTATTTTCTTCTAAACAATAATTTTTAACAATATTATTTTCTTTTTTAA
>CONTIG157 (71bp)
ACCATCAACTCTAATTTCTGAAAAATGGACAAATAAATCATCTCCTCCGTCGTCAGGAGATATAAAACCAA
>CONTIG158 (71bp)
GTAGGTAGAATAAAAAAAAGTTTAATACATAATAGAATTTTTAGTTTATTTAGTATAGGAGATCAACTTCT
>CONTIG159 (71bp)
GCAAATTTAATTTTAAATAAATGTATAGTTTGTAATCTCAGAAAAATATATGAAAATAATAAAATTTAAAA
>CONTIG160 (71bp)
GTGTGCGTATAGAACCACCTGTATCGCTTCCTATAGAACCCATAACACATCCACTACTAACTGAAATTGCA
>CONTIG161 (71bp)
ATATTAAGAAAATTTTTTAAAAAAAATGTTTTGTTAAAAAACGTGTATAAATTAATAATAAAAGAATATAT
>CONTIG162 (71bp)
ATAATTTTATAATAAAATTTATCATAAAAAAAAAACAAATTATCATTAGTAAAAAAAGAATATTTTTTTTT
>CONTIG163 (71bp)
TTTTTGCTACTAAACATAATAATTCAACTATACTTTTTTTTTTAAAAAATTGTTTAATTGAAAAGATAATT
>CONTIG164 (71bp)
GCCAGAAGCATATGATTTTGGTACAATGTGTTTAATTAATAATAATGATCTAATAATAATAGATACTATAA
>CONTIG165 (70bp)
AAACTAATACTATAGTCAGAATCTAGCATATTTTTTCTAGTGTTTTTAAGTAAAAAAAAATATAGTAGTT
>CONTIG166 (70bp)
AAATACTTTTAAAATACATTTTTAGAATAAAAATCTTTAGATTTTTTTTTAAAAAAATCAATTGTTTTTT
>CONTIG167 (70bp)
TGTTAAAAAACGAGTATAAATTAATAATAAAAGAATATACATATATATTAGTTGTTAAACAAATAAGTTT
>CONTIG168 (70bp)
CTTTATTAATTGATTTTAAATTGTTACAATTTATTTTATAAAAATTTAAAAAAAAAATCGTAAATACATC
>CONTIG169 (70bp)
TAAATTCTATAAAAGTAATTTCATTAATAGGTGGATCAACAAATTGTTTAATTCATTTATTAGCAATTGC
>CONTIG170 (70bp)
TTTTATTTAATTCGTTTAAAATAAAAAAAGAATCTGGTTATTACGTACAAAATTTTTTTGATGAATATGT
>CONTIG171 (70bp)
TTGGTCTTTTTTAAGAGAAAGACGATGGGGATCTCCTTTTTGTTATAAAAAAATAAAAAATAATAATTTT
>CONTIG172 (70bp)
TTTTGAATTATTATCATGAAAATTTTATTCTAAATAAATTATATAAAATAGATAAATAAAATCCAAATGT
>CONTIG173 (70bp)
AAAATGAACAAATAAATCATCTCCTCCGTCGTCAGGAGATATATAACCAAATCCTTTTGTATCATTAAAC
>CONTIG174 (70bp)
AAACGAAGAATTATTAAAATTAAGTGATAAAATTATAAAAAAAAAATAACAAAGCATTACTCTTTATTTT
>CONTIG175 (70bp)
CAAATAAAATTTTTGGAATATTAAATTTTTCAGATTGATTCCAGACAGTTTCAGTTTGTGGTTGAATACC
>CONTIG176 (70bp)
AATTAATTCAATATGTTTACAATTAACATAAATATTTTGCGGATAATAAATTGAATTTATTTCGTTAACA
>CONTIG177 (70bp)
GTTTATATAAAATATTTCAATATTTATAATAAAAAAAACAAAAAAATATTTAATGGTATTAAATCTAGAT
>CONTIG178 (70bp)
TAAATTACAAAAAAATAGATTTATGTTCTAATGAATTTGAAAACCCTTCTTTGTATTATTATTCTAGCAA
>CONTIG179 (70bp)
AAATTTAATTTTACAAACTCCATCAATTTCTGAAAGTAAAGCTCTCAATTTTGGTATTCTTGCTTCGAAC
>CONTIG180 (70bp)
AATAACAATTTTTTTTTTAATTGTATGCTAATCATAAATATTAATCATAAAATATTTTTATATAAAAATA
>CONTIG181 (70bp)
TAATTTAATAACAAATGATTTTGATTACTAAATTTTTAATACTAAATATTAATAATATTGATACTGATTT
>CONTIG182 (70bp)
CATTAACGTTTTTAATGTGTTATTTTATAAAACAACATAAAAAATTTATTTTAATTTTAAATTTTAATAA
>CONTIG183 (70bp)
GTATTGAAATTTTTTTTAATTCTAAAACAGCTTGATATATTGTCTTATCAATACCTCTTTTTAAATCCAT
>CONTIG184 (70bp)
TTAAATTTCCTAAAAAGTTTTTACCTAATTGTAAAAATTTAATCAAAAAAAACAAAATATTTAGTAAAAT
>CONTIG185 (70bp)
GGATTGTTTTGATCTAAAAATTGAGACAACTCATTATTACAAAAATAATCTTTTAAACCTATAGTAATTA
>CONTIG186 (69bp)
ATTATTGATAATTAAACATGATTTTTTTCGACTAATTCCAATAAACCAAGTTGAATATTTTGAAAATAA
>CONTIG187 (69bp)
ATAATAAAATTATAAAAGTATTTTATATGGGTTCAATGTTTAGATATGAAAAACCTCAATTAAATAAAA
>CONTIG188 (69bp)
TGGAACAGGAATAATAAGTTCATGTCATTATATTAATAACATTAATAAAAAAATAAACTTATTTAATGT
>CONTIG189 (69bp)
TTTTGAACAAAAAAATTTATATTTTTTACAATTAATAATCATACTTTAATAAAATCAATATTTAAATTA
>CONTIG190 (69bp)
TTTCGATATTTTAGAATCAATCGAAAATATTAATTTAATAATCTTAAATTGTAACGTATTATTTAATTT
>CONTIG191 (69bp)
ATAGAAAATAAAACTATTATTGAAATACAAAATTTATTTAATTCAAATTTTGCTAATGTTCAATCTCAT
>CONTIG192 (69bp)
AAAATTTTTTCTTCTTTTTTTTTTAAAAGAAACTTTTTTTAATATAAAGTGTTTTAATAAAATTGCATT
>CONTIG193 (69bp)
TGTTTTGTAATTTCTGACATAAGTCCTATTTCTGGATTAATTGCCTCGGGCTTATACCCTTCACCTTTG
>CONTIG194 (69bp)
ATTCATTGATTTATAATAAATTAAATAATTTTAAAGTTTTTTCTGGATTAGATTGTCATGGTTTGATTA
>CONTIG195 (69bp)
TTAAAAATAACTTTATTTTTCTTAATACGCATTTACTTGATATTAAATTATATATTTTACGATTTTTAT
>CONTIG196 (69bp)
TTGATTGTCAAAGAAAACTAACATCCTTTTCCACTCAAATAGAATTTAGGGACCTTAATTTATGATCTG
>CONTIG197 (69bp)
AAAAAAAATACTATGAATTAAAAAACTATTTTTAATAAATAATTCTAAAAAAGAAATAGAATTTCTAAA
>CONTIG198 (69bp)
AAAAAATTGTTTAAATTAATTTTAAAATCATTATAAAAACTAACATGTAAATTAAATTGATTTAACGAA
>CONTIG199 (68bp)
ATTATTTAAAGTTTGTTTTCCAAACTTTAAATTTGATTCTAACATAAATCCAATTATATTTGTTTTTA
>CONTIG200 (68bp)
CCTGGACCTATACAATCAGGATTTATTAATGATTTTATTAATACAAAAAATAATATTGTAAAAACATA
>CONTIG201 (68bp)
ACAAGATATTACCCAAAATAAAAACAATATATTTTTATGATTAGATTGATGTAAAAATATTTTTTTTT
>CONTIG202 (68bp)
AGAACCTACTCTTATTACGGCAACATCTCCTGCTAATTTTGCCATTCGTTCTTGTAATTTTTCTTTAT
>CONTIG203 (68bp)
ACATTTTGTTAAAAATAATTTTTTCTTTAAAAAAAAAATGAAAAAAATATATTTTTTATGCATTTGTG
>CONTIG204 (68bp)
TGTTAAATAAATAATAGTAACATGAATCAGATATATAAAATCCATAATTAAATAAATTTATAATATTC
>CONTIG205 (68bp)
TAACTAAAATTGATGAATAATTTGAATCTCTACCACCAGAATTTAATAAAAAATTTAACATTAAAATC
>CONTIG206 (68bp)
ATTGGGATAATGTTGAAAATTGTATTCTTTCTAATGAACAAGTCAATAATAATAAAGGATGGAGGTCA
>CONTIG207 (68bp)
ATTGGAATATCAACTCAATGTATTCTTCAAAAAAAATTAAAAAATATGAAAATAGTAATTAATAACTT
>CONTIG208 (68bp)
GGGTAAATTACCTCCTCTTAAAACACAATGACAATTTAAATTTCCTAAAGAACTAGTAAATTTTCTAT
>CONTIG209 (68bp)
TTGAACTATTCTAGCTTAGTAACTACTACTACGCATAAAACATCAAGAGGAATTAAAGGAGGAATTAT
>CONTIG210 (68bp)
TAGTTTATTATAAATTTAAATTTTTTTATGCAAAAATTAAATGATTATTTTATTTAAATTTGATAAAT
>CONTIG211 (68bp)
GGCGAAATGGTAGACGCAAAGGACCTAAAATCCTTGGTTTTCTAAACGTGTCAGTTCAAATCTGACTT
>CONTIG212 (68bp)
AATTATCGAATAAAACCAATAACAATTTATTTTCATGTAGTAACTAAAAAAGTTGAAGAATATATTAA
>CONTIG213 (68bp)
GTAGAATAAAAAAAAGTTTAATAGCTAATAGAATTTTTAGTTTATTTAGTATAGGAGATCAACTTCTT
>CONTIG214 (68bp)
TTTATTAAATAAGAAAAAAAATAACTTTTTATTTTTTTTTTTTTTTTTAAAAAAAAAATGAAAATTCA
>CONTIG215 (68bp)
ATTATTAAAATTAAACAAAAAAAAACACAAAAATGTTATTGCGCATGTAGTATCTGGATCAGAAAATA
>CONTIG216 (68bp)
GATAGTTTATTTGCTTTAATTAAATTTTTTTTATTATTAAATAAAAAAAAAATTTTAGTACCAGATGT
>CONTIG217 (68bp)
CTTTTTTATGGTAATACGTTATCATGTATAATTGATTTTAATAATTATAATTATTTGCATTTTAAAGT
>CONTIG218 (68bp)
GTTAAAATATATAAAATTATAATTTTCACAAAAATAAATTTTATTAATATTAAAACATTTTTTTGAAG
>CONTIG219 (68bp)
CAAATTAAATTAATGTATAACATTTGTTGGCCTCCGTAAAACGATGAATAAAAATGAGTACTTAAGGT
>CONTIG220 (68bp)
ATAATAGTATTAGGAGGTGGATGTATACAATATTTAAAAAAAACAAAAATTAAAAAAAGTATTTTAAT
>CONTIG221 (67bp)
TATATTTTTATGATTAGATTGATTTAAAAATATTTTTTTTTTAAAAAAAGATTTAAAAATTAAACTA
>CONTIG222 (67bp)
AAATATCCATTCTCTTTCATACCCGAAATTACAACATTTGAATTAAGTGAAATTATAGAATTTTTTC
>CONTIG223 (67bp)
AAATATTAACTAAATATCTATTGTTTGAAACAATTAGTCCAAAACGAAAAAAAGTATTAAATACTTT
>CONTIG224 (67bp)
TAAACATATATTGTTACAAAAAATTATTATTTAAAAAAATTTTTTTACATGAAATAGTTAAAGATTA
>CONTIG225 (67bp)
TTTCTAAACTATTTATATTGTTAGAAAACAATTCTTTTAAAAACATATTATATTTGAATAAATAATA
>CONTIG226 (67bp)
TAATTGCCTCGGGCTTATACCCTTCACCTTTGAACTATTCTAGATTAGTAACTACTACTACGCATAA
>CONTIG227 (67bp)
AATTGGAGCGATAATTTAGTAGGTTAAAATGTTGGCTTGTCACTTCAAATATCGCGGGTTCGATTCC
>CONTIG228 (67bp)
TGGACTAGTAATTATTTGTATAATAAACCAATATGGGTAAGTGAAGATTTAAGAGATGGAAATCAAT
>CONTIG229 (67bp)
TGTTAATTAACAACAATGATTTTTTATATCCGTAAGGACTTTGCTCTAAATTTATAATTTCTACATT
>CONTIG230 (67bp)
AAATTGAAAATAAAAAATGTTTTTTAATACCTACATCAGAAGTCATATTAAATTCTCTTTCGTTTTT
>CONTIG231 (67bp)
ATTTTTATTAATGAAATTTACATTATTTATAATTTTAGGAACTTTATTTAAATTTAAATAAGAAGTT
>CONTIG232 (67bp)
ATCTGATTTAATCGTTAACATTTCTTTTATCAAAAATGCGGCACCATATGCTTCTAAAGCCCAAACT
>CONTIG233 (67bp)
GTATAAAATATTGTTCTAAGTATTTAGCTTTAATACTATCATACAATATGAGTATTGGTATCAATAT
>CONTIG234 (67bp)
TTTTTTTAAACAATTTTAAATTTTCAAAAATTTTTTTAATTAACAATTGAGCAGAATCAAAATTTCC
>CONTIG235 (67bp)
AAGTAGTATTTTTTCTACTTCAGTCCACAATGTTTCAACTAATTTTAAACATAGATTTTTTTATAAT
>CONTIG236 (67bp)
TACAATTTAAAATTTTGTTTTTATATATTTTTTTTTGTTTATCTAATAATTTATTGCAGAATAATGA
>CONTIG237 (67bp)
AAGATTTGGAATACTTAAATTTAAATTTATTTTTGTAACGTGATTAAGTACAGTATGATGCATTTTA
>CONTIG238 (67bp)
CATTATTAAATTCTATTTTTTTTAAATATTTATTTTTAATATACCTTTTATTTGATGGTTTTGAAAA
>CONTIG239 (67bp)
TACTGTTTTTTTTTTAGTAAAATCTTCTTATAATTTATATTTTTGTTTATTTTGCAAGATATTTATT
>CONTIG240 (67bp)
TTATATTATAAAAATAATTTTATCAAAAAAATTAAAAATTTAAATTTTTTTCTAAACATTAAAAATA
>CONTIG241 (67bp)
AAATAATTAAAAAAAAAATAATAATTAAAATTTTAAATATTTTAAGCGTTAAATTGAATGAAAATGT
>CONTIG242 (67bp)
ATTCCAAACATTAAAGTTGGTTATTTATACATAAGACCTTTATCGCTTGGAATAGGAGGAATTCTAG
>CONTIG243 (67bp)
TATAATTCTTATATAAAAAAAGTATTTTATGTTTAAATTTTTTGAATATAAAAAAATAAAATTTAGA
>CONTIG244 (67bp)
TAATTTTTTAAAATTAGGAATAGAATTCATATATGAACCTGAATACCTTTTATAATAATTTAGTTCT
>CONTIG245 (66bp)
TATTTATTAAAATATTATTGTTAATCATTCCTATTCCTATTATTTCTATCCATAAATTATTATAAA
>CONTIG246 (66bp)
ATAAATTTAAAAAAAATATTTAAAATTTTTTCTTTACCAATATATATATATGATTTTAAAAAAATA
>CONTIG247 (66bp)
CATTATTACAAAAATAATCTTTTAAACCTATAGTAATTATATCCCTATTAACAATAAAATCTAAAT
>CONTIG248 (66bp)
TATTATAATATGAATCTATTAATTCGCCTAATAATTTGTTTTTATACCCGTTTTCATTTCTTAGTT
>CONTIG249 (66bp)
GAATTAATAAGTTTTTTTATATTAAACAGTATATTATGATTATTAAAGAAATAAATAGTAAAATAA
>CONTIG250 (66bp)
ATTCCAATATTATTATCAATTTAATTAATGATTTAAACAAAAAAATAATTAAAAAAAAAATAATAT
>CONTIG251 (66bp)
AAAAATTTTTTTTAAAATTTAATAAATAATATATAATACTAGACTTAATTATACAACTATCTAAAT
>CONTIG252 (66bp)
TTAATGTTGGATTTAGAAAATTACCTAATATTGGTTTTCCTTCTATAGAATAAACACTCTTTTTGT
>CONTIG253 (66bp)
AGTATTCCATTTAATATATTTAAAAAAAAAATTAATAAATACACAAATATTTTTATTAATAATATT
>CONTIG254 (66bp)
ATAGGACCAGAAATAATTAAACAAGTAATTAAAATAGTTAAATGATGTATTTATACTGGTTATAAA
>CONTIG255 (66bp)
TTAAAATTATCTTTTTTTAAAACAAAACACAAAATAAAAAATTTATTATTTAAGTGTTTACCATCA
>CONTIG256 (66bp)
TTAGTTTTAAAAATAATTTATTGGATAACAAATTGTCTGTATTATTGCCACAAAAAATATTTGTAT
>CONTIG257 (66bp)
ATTCGTTCTTGTAATTTTTCTTTATCGTAATCTGATAACGATTATGATATTTGTTTTTTAATAGTA
>CONTIG258 (66bp)
GTTATATAATAGATTAATTTTTCCAATACTAAATGAAACAGGTTTTTTAATAGGTATAGGTTTAAA
>CONTIG259 (66bp)
GAAAGAGCAAAAAAAATAAGAACATATAATTCTACCAATAATAAAATTACAGATCATATAAATAAC
>CONTIG260 (66bp)
ATAAATAAAGACAAAATAATAGTAAATGAAGCATAACGTTTAAAAATTAAGATTGCGTCTTTTTTA
>CONTIG261 (66bp)
TGCAAAAAAAATTATTTTAATTTTTTTAATTAATAAATTATTTCTTTTATTATAAGAATTTTTAAT
>CONTIG262 (66bp)
AATAAAACATATATCTACAATATCATTTATATTATTTTTATTATTACCAAAAAAACAATTTATTTT
>CONTIG263 (66bp)
TTTTTAAGTAAAAAATTAAATATTTTATATTTTAATACAAATATATTTTTATTAAATTATAGAAAT
>CONTIG264 (66bp)
AAGACATTATTCCTATATCATTATCAATAACAAAATTATAATTTTTTAAATTACTATTTCTAATTT
>CONTIG265 (66bp)
TTTTAAACAAAAAAAATGTTTTTTGTAAACATAAAATTATTTTGTACAAAAATATTAGCAAATATT
>CONTIG266 (66bp)
TAAAAAGTAAGTCTAAACCTAGTAAATCATTTGAATAATTATTCTCAGCAAAACCTTTTTTAATTG
>CONTIG267 (66bp)
TAAAAAAAAAATTGTTTTTTTACATAAATTTTTTAACAAAGAAAGTGAAAATTATATTATAAAAAA
>CONTIG268 (66bp)
TTTAAATAATATGATATTTTCATTAATCCTTTAACAATTAACGGTCTTGCAATCGCAGTTCCTAGC
>CONTIG269 (66bp)
GAATTATTTTTTTTTTTTTTTTTTAAATAAAAATTTAAATTTTCTAAAAGATTGATAAAAACTAGT
>CONTIG270 (66bp)
AAATTAATCAAAAAATATAATTCAAGTTTTAAATTAATTATTGAATTAGGAAGGTTTTTTTTTGGA
>CONTIG271 (65bp)
AAAATTATATTAGTATACCTTTTTATAACAAATTCAATTTTTCCTTTTGCTCTGTAATTTATTTT
>CONTIG272 (65bp)
ATTGTTTTTAACAAAAAAACAGGCTTGATTATTAATAAAAGACCATAAAGAAGTTGTCCATATTA
>CONTIG273 (65bp)
GATTTATTGTGTCCATTTTTAAGAGGAGGAAAAATTGGTTTATCTGGTGGAGCAGGTGTTGGTAA
>CONTIG274 (65bp)
AGTTCGTTACTAGAATATCCAAATTTTTTTATGTTTAACACATCCGTAAAGTATTTTTTATATTT
>CONTIG275 (65bp)
CATAATAATTAAAATTACTTTTTGAGAACACTATTTTAAATAATTTTTCTTCGGGAGAAAACTTT
>CONTIG276 (65bp)
TATTGATGATATTGATATTAACAATAAAAAAAAAATAAAATTTAAAAATTTAGAAAAACCAAATT
>CONTIG277 (65bp)
TTGGAATTACACCAAAGATTTTAAAAATTTATACTAATATAATTGGTATTTCTAACAAAAACATT
>CONTIG278 (65bp)
CTTTCAAAAACAGAAAAATTATAAAAACAAGAATTTAAGTTAAAACAAAAAAAAACAATTAAATA
>CONTIG279 (65bp)
TTAAATAATTTTATAAATTATAAGTTTAAAACAGTATATTTAATCATTATAATACATTTTTTAGA
>CONTIG280 (65bp)
CAATGCAAAGATTAAAAGAAGCAGCAGAAAAAGCTAAAATTGATTTATCTAGTGTAGAACAAACT
>CONTIG281 (65bp)
ATAGATTAATTTTATTAAAAAATTTATTTTTCATATTATGATATATATTTCATTTTCTAAAAAAA
>CONTIG282 (65bp)
ATGTTAAAAACAAAATATTAAGAAAAATAAAATGAATTTAAAATGTCCATTTCATAATGATACAA
>CONTIG283 (65bp)
AATTTAAAACAAATTGAAATTCTTCCAACTTCAAATATAAATATAACTTTTTGTTTTTCAATAAT
>CONTIG284 (65bp)
ATTTTTTAACAATAATTTATATTCTAATGTCTGCATTGCAAATCTTTTGTTTAGCATTTGAGATC
>CONTIG285 (65bp)
TACTTTAATTGATCAGGTTTACTTTTCATTATTTTTTTTAAATAACCATATTTTTACACCTAAAA
>CONTIG286 (65bp)
AATTCTGGATGTTCAATGTGTCTAGCAATGAACGAAGATAAATTAAAACCAGGAGAAAGATGTGT
>CONTIG287 (65bp)
ATGCAACTAAAAGTATTTGCTTTGGTTTTAAGTTTAAAAATAAAAACTCATTTCAAGTAAAAGAT
>CONTIG288 (65bp)
TAGTTTTATTACTCCAAATTTCTTATCTTTAAAAATAATTTCAGATTCTGAACCAGAAGATCTTG
>CONTIG289 (65bp)
GAATGAGGTATTAATATACCTTCTTTTTTTTTTCTATAGTTTAAATTAACTTGTATTTTATAAAT
>CONTIG290 (65bp)
AATATGGTTTAAAATGGTTGATTCTTCACCTGCAAAAAAAATTATTTTAATTTTTTTAATTAATA
>CONTIG291 (65bp)
ATAAATGCATTTTTTTTCCAAATCGAAATATTCTTAAAATTAGTTTTAGTTTTTTTTTTATTTTT
>CONTIG292 (65bp)
AAATTTTATTTTGTAAGTATACCAAAAAAGCTATTTTATTCGGAATAGAAAAAATGTCAAAATCA
>CONTIG293 (64bp)
AATATGATATTTTCATTAATTCTTTAACAATTAACGGTCTTTCAATCGCAGTTCCTAGCAAATA
>CONTIG294 (64bp)
TAAATTTATTTTTTAATAACAATTTATAATAAAAATTTAAAAAAAATTTATTTAAATCATGTGT
>CONTIG295 (64bp)
AATATTTATGTTGATCAATTATTAGGAATTGATACTAATAATATTATAAATATTATTAAAAGTA
>CONTIG296 (64bp)
AAAATAGTTAGGTATTTTTTCATATAAAAATTTTTTTTAAGAAAATTATTAGTTCCAAAAGTAC
>CONTIG297 (64bp)
TTTATAAAATTAAAGCATGACAACATTGTAAAACTTATAATTTCACAATTAATTATTAATAAAA
>CONTIG298 (64bp)
TTGGGGTTATAATTTAATTGTTAAAATATTTGGTTTGCAACCAAAAAAAAAGAGTTCAATTCTC
>CONTIG299 (64bp)
AGATTGATTAAAAAGATTGCCAATCTCATCAAAAGAAGAATTTTTTGAATTAAACAATGAATTA
>CONTIG300 (64bp)
CTTTTATTATTTCAAAGTATGATATTTTTTTATTTTTAATTAACAAAAAAATTTTTTTTATACC
>CONTIG301 (64bp)
CCTTTGTATCCTAGAATACCTGCTGAATACCATACAATGGTATCACCAAATACATCAGTTAAAG
>CONTIG302 (64bp)
TTTAAATTTAAATTTCATACCTTTTTTTTTTAAAAAATTTTTATTATAATGTAATTATTATTTT
>CONTIG303 (64bp)
AAAATAAAAATAAAAAAAAATAAAAAATTATAATTATAATTAAAAAAAATAAATTTATTAAATA
>CONTIG304 (64bp)
TTTCAATAAATTTAAATATTTTTAGTATTCATATAGAATCAATACAAGAATTATTTAAAATTTT
>CONTIG305 (64bp)
TTTTATTAATAATTTTTTTATTTAAACTTTTTTTTTTGGATTTATTCAATAATAAATGAGTTTT
>CONTIG306 (64bp)
ATTTATAATATATCCTTTTTTATATTGTAATAATTTTTTATTTATTTTTAATTTAAAATTTTTA
>CONTIG307 (64bp)
AGTTAATTTTTCAGTAAAAACACCTAATTTATCTTTATAAACTTATAATAAAATTATACTATTT
>CONTIG308 (64bp)
AAGTATTTTTTTGTTTTGATGGTGATCATTCTGGATATTTAGGTGTATTAAAATTATCTTTTTT
>CONTIG309 (64bp)
TTTTTAAACAATATATATTTTATAAAAATTGGATTTTTGTTTAAAATAAAATTTTTAATAAATT
>CONTIG310 (64bp)
TATAGTTTGAGCTAAAACAGTAGCTGTTGTTGCACCATCACCTGCAACATCTGATGTTTTTGAT
>CONTIG311 (64bp)
TTGCATCGTCGCCAAATTTTATTTTCTTATAACCCATTTTCTCTAATAGCGATTACATCTTCTT
>CONTIG312 (63bp)
AAGATGTAATTGTTTTAAATATAACACCATTGTTTATTTATAAATAAAAGTAAATAAATGTTT
>CONTIG313 (63bp)
AAACCTTCCCAATTTTTACTTATTTTAGAATAAAAAATATTATTAAAAAGAATTAGAGATAAT
>CONTIG314 (63bp)
ATAACTATTACTTAAAAAATTAATTATATAAATTATAATTTTCTTAGATTTATAAAATTTAAA
>CONTIG315 (63bp)
TCAATTAATAATATTTTATGATTAGAAAATTTTAATAATTTTAAAAAAAAGAATAATAATTTA
>CONTIG316 (63bp)
TATTATGTGATGTTATAAATAAAATATTTTTAATATACTTTTTTTTATAATAAAAGTCATTAT
>CONTIG317 (63bp)
AATATATAATTTTCTCAAACAGCTTAAAGTATATTCTTTATAATAACCAAATTTTGAAATAAT
>CONTIG318 (63bp)
ACTAAATCTTTGTACGATAGTTTATATTTTTTTAACAATTTATTTTTTTTTTTATTTTTATTA
>CONTIG319 (63bp)
AAAAATTTTAAATTTGATTGTAATCGCTAAAAATTATTTAGGATATTTAAAATTAATTAAAAT
>CONTIG320 (63bp)
GCCGACATCGAGGTGCCAAACATTGCCGTCGATATGAGCTCTCTGGCAATATTAGCCTGTTAT
>CONTIG321 (63bp)
TCGCCTACTAATTTGTTTTTATACCCGTTTTCATTTCTTAGTTCAGAATTTTTTACTATCCAT
>CONTIG322 (63bp)
TTTTAATATTTTAGAATTACAACTTAATAATTTGTTTTGAATGCTATAAATATTTTTTATAGA
>CONTIG323 (63bp)
TTACTAAAAATAAAATTATAGTTAAAATAACATAATTTGATAGAAAATTAAAAAAAATATCAT
>CONTIG324 (63bp)
AAAATAATTATAATTGTTGTTAGAATTGGAACTTTATTTTTAAGTAAAAAATTAAATATTTTA
>CONTIG325 (63bp)
ATATATTTATTTAAATTAAAAAAATTTTTTTTTGTTCTTATAGGACCTATATCATAAAATAAA
>CONTIG326 (63bp)
TCTGTTTTAATAAATTATTTAAACAAAAATAAACCCAGTTTTGTTTCTGTAACATTTGGAAAA
>CONTIG327 (63bp)
TGTTTATGTTATAATTAAAAACTTATATATTAATAGCTTTGTAATTAATTATAATAATATAAT
>CONTIG328 (63bp)
AAGTTTTTATTTATTTTTAATGATCAAAAAATCCTTCTTTATTGAATAAAAAGTAAATATTTT
>CONTIG329 (63bp)
AAATTATAATTTTTACAAAAATAAATTTTATTAATATTAAAACTTTTTTTTGAAGATTTAATA
>CONTIG330 (63bp)
TTCTTCAAAACTAGAAGGAATAGCGAGAAATATAAGTAAACATTCAGGAGGTGTAGTGATATG
>CONTIG331 (63bp)
ACTTATATTACTATGCGAATATTTAATATACCATATATGTACTGCTGGATAAAAAAGTTTTAT
>CONTIG332 (63bp)
GAATTTTTAAAATCAAAAAAAAAATTATATTCATTTTTTATATTTATAATCAAAATTTTTTTT
>CONTIG333 (63bp)
ATTTTTTTAAAAAAAAATTATTAATTAAAAAGTTGTTTTTTAAAAAAAATAAATTAGTATAAT
>CONTIG334 (63bp)
TGGTTTAATATCTAACTTATTTAAATTAACAAAATTGGGTAAACTACCTCCTCTTAAAACAAA
>CONTIG335 (63bp)
TATATTTTTTTTACGATTACCAACAATAGAAGTTGCTGTATCGGAAAAAATTCTTCCTCCTTT
>CONTIG336 (63bp)
TGGTTTAACCATCTTGGTCCTCTCGTACTAAAGACAGATTCTCCCAAAATCTCTACATCTACG
>CONTIG337 (63bp)
AAATTAAATAATGAAACGCTCATTTGTTCACCAATACATAACACAAAATCTATACTTTTTAAT
>CONTIG338 (63bp)
AGTTTTTATCGGATCATGTACAAATTCTAGATTTGAAGATTTACTAGTTTGTTCAAAATTATT
>CONTIG339 (63bp)
ATGGTTTTAAAATTTTTCAACTTGATGAACCAACAATTAAAGAATGTTTACCTATAAATATTA
>CONTIG340 (63bp)
CAATTAATTTTCTTGTTAAATAACCAGAATTTGCAGTTTTAAGTGAAGTATCAGTTAGACCTT
>CONTIG341 (63bp)
TTAATCCACCTTCATAAATCAAATCAACAATTAATTTTAATTCATGTAAACACTCAAAATATG
>CONTIG342 (62bp)
TAAAAATTCATAATTTTTCATAATCATAGAATTATAAAAAACTGTTGAAGCAATTTTTATAT
>CONTIG343 (62bp)
TTGTAAACCTTTTTTTATATAAGATAAAAAATCTTTATTTTTTCTATTAAAATTTATTTTTT
>CONTIG344 (62bp)
TCTTATTTCCATTTCTACTAATTCTAATAATTCTTTATCTTTAGCACAATCTGCTTTGTTTA
>CONTIG345 (62bp)
AAAATTTTTGAAAGTCATGTAATTAAAAAATATAATCACTTATATATTTTATATATAGATAA
>CONTIG346 (62bp)
AATGAATCTAATATATTAGCACAGTTTTGTATTTTTTTAATTTTAAAATTAGTTAAATAAAT
>CONTIG347 (62bp)
TAGCTCTACCTGTTCCTTTTTGAATCCAAGGTTTTTTACCAGACCCAAAAACCAGTCTTCTA
>CONTIG348 (62bp)
TAACATATTTCTGTTTTAATAGATAATTTTTGAGAAGCAACATAAAATGCTTTTTTGACAAT
>CONTIG349 (62bp)
TCCTTTGCTAACAAAATAGACTAGTTCAATTATTTGTTGAGGTAAACAAATTTTAAAATTAT
>CONTIG350 (62bp)
GTAATTACATCGCCTGGTGATAATAAACAGGAAGGAATTTTATTTTTAAAATAATTAACAAA
>CONTIG351 (62bp)
TAATCAGCTATCATTGCTTCTATAGTAATTAATAATCCTCCAATTGAACCTGCTGATTGTAG
>CONTIG352 (62bp)
AGATGTATTCAACACGAAATTGATCATTTAAATAGTAAACTAACATTAGATTATTCTAAAAT
>CONTIG353 (62bp)
CTCGCTTAATTATATTAATAACAATAAAAATTTAGATTTAATTTTAATAGGAAAAAAAGGAA
>CONTIG354 (62bp)
AAATTAAAAATAATTTTTTAAAAAAAAAAACTTTAAAATATTTATTTTATAACAATTATATT
>CONTIG355 (62bp)
TGTTATCATTTTGATTTACTGGAATTTTCATATGTAAAAAATAATTCATTTTCAATTAATGA
>CONTIG356 (62bp)
AATTTTTTTTTTTCCAATAAATTATATTGTTAACATAAATTTATCATTTAATATTATAAATT
>CONTIG357 (62bp)
TTTTATACAATCTAATAAAGCAATTGAATAATGACTATCTGCACCTGGATTTATAATAATAT
>CONTIG358 (62bp)
CTAAATTTAATAATTCATGTTTAACTTTATTAAAATCTTTTAAAAAATTTGAATAAAAAATT
>CONTIG359 (62bp)
ATAGTTTTTAACGTATGTATAATTACTTGTTAAATATGAATTTCAAATTGGTAAATTATATT
>CONTIG360 (62bp)
AAAATTACATTCTTTTTTTAAATTTTAGAAATTTGAAAAACGAAAACATCTAAATAATTTTC
>CONTIG361 (62bp)
AAATTTTGTCTTGATAAAATTAGAGACGAACATTTATTTAAAAATTTCAAAATTAATATCCA
>CONTIG362 (62bp)
AAAAATATAAATATTAATATTTGTAGTGGATGTCCTCCTTTTTTTACAAAAAAAAAAAACAA
>CONTIG363 (61bp)
ATTAGGTAAAATGAAACTTATAGGTGATATTTATATTGCTAAAACTTTATTTATTTTCAAA
>CONTIG364 (61bp)
TTCTAAAAAATTATTTTTAAAATTTGTTGAATAATTACATAAATTAAAAATTAACATAACT
>CONTIG365 (61bp)
TTTAATAGCAGAATTAATTAAAATTTTATTTGTTTCTTGAAAAGAAGCAGCAGAAAAAAAC
>CONTIG366 (61bp)
AAAAATTAAAAATTGTAAAAAATGTAATATTTGTACAAAAATTCGTCCATTAAATTTAATT
>CONTIG367 (61bp)
TGTTGATATACTAAAAAATATATTGGATAAATAACCAGATTCCATTGTAAAAAAATGGTGT
>CONTIG368 (61bp)
TTTTACTAAATTTTTTTCTTCTTTCATTATCGAAAGTGTTTTAGAACCCTAAGGCCTTCTT
>CONTIG369 (61bp)
CACAACAGACATAGAAATTGGAATATCAACTCCATGTATTATTCAAAAAAAATTAAAAAAT
>CONTIG370 (61bp)
TTGAATACCTTTTATAATAATTTAGTTCTTGCCATTGTTTTACGATTCCTAAACTTTGATT
>CONTIG371 (61bp)
ATTAAACCAATAATAGGAGCAGAATGTTTTTTATTAATAAATGGAAAAATTTTAAATTTGA
>CONTIG372 (61bp)
TTATTGCGGGGTAGAGCTGATGGTAGCTCATCGGGCTCATAATCCGAAGGTCGATGGTTCG
>CONTIG373 (61bp)
TTTTAGAATATAAAAATGCAATTTTTATTTTTTTAAATTTTAATTTTTTTTTTTTTTATAA
>CONTIG374 (61bp)
AATTTAAAGTTACTCAAGATATGCTTGAAATAATAAATGGAAGCAATTTATGATAGGAAGA
>CONTIG375 (61bp)
TTTTTTTATATCCAAATTTGATTCCATCATACCTGCAAGAATTTGTATAAAATTCTTTAGA
>CONTIG376 (61bp)
ATTAATACGTTAAAAATCATTTGGAAATTTTTAATAAAGTTTTTTAACAAAAATATTTCTT
>CONTIG377 (61bp)
TTTACCAGAAATTGACGCAGTAATCATGCTTCCGTTATTTACTTTTACTTTAAATATACTA
>CONTIG378 (61bp)
TTTTGGAATGAAATATTCATTACTATCAACAGATATAATTTCAAATTCAATAGAAGTTGTT
>CONTIG379 (61bp)
ATTTATCTTTAAAAGATCTACAATTTGAAACTAATAACTTACCAACGTTATCTGATTTAAA
>CONTIG380 (61bp)
TTTAATAAAAAATTTAACATTAAAATTTATGAATAATTCTTCCTTTAAACAAGTTATACTG
>CONTIG381 (61bp)
TAATTACTTCTCCAAATTTATTAATTACACATTTACATTTTTTGAATTTTACAAATCCAGA
>CONTIG382 (61bp)
ATTAACTAGAACTTTTTTAATGACAGGATTATCAGATGCTATTCCAATGATAAGCGTAGGT
>CONTIG383 (61bp)
CTGATGTAATTGTTATTCCACGTTCTTGTTCTTGTTTCATCCACTCTGTAATTGTGTTACC
>CONTIG384 (61bp)
ACATATATTAAATGATTTTTCTTTTAAAAAGTATGAAATTTTATCTTTATTTGAAATTAGT
>CONTIG385 (60bp)
ACAAAATTTGATAGAAAATTAAAAAAAATACCATTAAATTTAAAAAAAACATATAAAAAA
>CONTIG386 (60bp)
GATTAAAGAACTTTTTACATCTTTCTTCGTATAAATTAATGGTAAAGAAATTTTAATCTT
>CONTIG387 (60bp)
AAATATATGTATAAAAACAATATCAGAATTAAATGTTATCAAACAAAATTTAATTAATGA
>CONTIG388 (60bp)
TAATAGTTTTAAAAAATAATTGTTATGAAAACATTAAAATATACAACGAACCCAAAAAAA
>CONTIG389 (60bp)
TATAATAATTTTGAAATGTTTTGTAAAATATTTAGTGTTTTTTAAATTAAAAAAAAAAAA
>CONTIG390 (60bp)
TTTTAAATATTTTAAAGTTAATTTTTTGGGCTTCTTCCATTTCGTTCGCCACTACTTTGG
>CONTIG391 (60bp)
ATTTTTTTTTAATTTTAAAATATATCATTTTAAAAAAAAAATTTAATATTAAATAACTGT
>CONTIG392 (60bp)
AAAAATTTATAAATAATAAATATAACTTTATTTTTAAATTAATTACAAAAAATAGTGTTA
>CONTIG393 (60bp)
TAAACGAAAATGATTTATTCTATGAAAAAACTAAAAAATATTTTATAGAAAAAAAATTTT
>CONTIG394 (60bp)
AAAATTTGATTTAATTAAATTAGTTAAAATAAAATTTGTAGTAACATTTAATATGTTAAT
>CONTIG395 (60bp)
ATTTTTGATAATTTTTTCATAAAAAAAAATGAAATATAAATTCAAGTTCTTTTCGGTATA
>CONTIG396 (60bp)
AAAAGAATATTTAGAACTAAATTCTTTATTAGCAGTAACATAATTGCATTTATTTTTTAT
>CONTIG397 (60bp)
ATAAATAAAAAAGATTTATAATGTAAAATAGTAATTTCAAAATTATTATGTAAGATAAAA
>CONTIG398 (60bp)
AAATCAAAATAACTAATAGTTATTCCTGAATATGTTGATATTATAAATCCTATTTTTTTT
>CONTIG399 (60bp)
CAATATTTTTTTATTAATTTTCATTTTTTGATTTTTTAAGTAAAAATTTTGTACTCCTAA
>CONTIG400 (60bp)
GCTAATAAAGCAGACCCAATAATAATTTTAGTATTATTACCATCAAAATCATATTCTGTT
>CONTIG401 (60bp)
TTTAATTTTAAATATATTGTTAACATTTATATTATAAGAATTAATTCTTTTTTTTTTAAT
>CONTIG402 (60bp)
TTTAATAATAAAAATACAATTATGTTTTTATTTAATTTACAATTGTTAATTTTTAAAATA
>CONTIG403 (60bp)
ATGGATATCCTATCAAATTTGAAAATACTGTAATATAATCACATTATTCTGAATAATTTA
>CONTIG404 (60bp)
ATCTAGTGTAGAACAATCTGATATTAATTTACCATATATTACTGCTGATAAAAATGGACC
>CONTIG405 (60bp)
TAGAATTATTAATAGATAATTATAATAATATTATAGTTATATGGTGTATTGAAAATATTG
>CONTIG406 (60bp)
CTAGGAGATATTGCGATTCTTAAACTAAACACTTTATATTTTATGTGGGTTTAGTAGGAG
>CONTIG407 (60bp)
CTTATTTTATAGGTTATAAAGAAACAATATTGAAATTTTATAACTTTTCTTTTAGAACCA
>CONTIG408 (60bp)
AAGTGGAAAACTTCACGTTGGACATGCAAGAAGTTATATTATACCCGATGTTATTTCTAG
>CONTIG409 (60bp)
GAAAAAAAATTTTTTTTTAAAAATCCCAAATTTGCAAAGAAATCGTTTAATTTCTAAAAT
>CONTIG410 (60bp)
AAAAAAAACATATTTATTAATTATAATTTAATTTATTTTAAAGAATATTTTTTGCAATTT
>CONTIG411 (60bp)
ATGATTTAGGTGGCGGTACATTTGATATCTCAATAATTGAAACAGCAAATGTTGATGGAG
>CONTIG412 (59bp)
ATTTATTCAATAAAATAACTATTTTAAGAAAAATACTTTTAAACTACATTTTTAGAATA
>CONTIG413 (59bp)
TTAAATCTGCAAAATTTAAAATATTATTTGAAACAGAATTTCCATGCATTCCGATCCAA
>CONTIG414 (59bp)
TAAAAAATTTATTATTTAAGTGTTTACCATTAAAATATGATCCCGATTTATACATTAAT
>CONTIG415 (59bp)
AATATAATAAAAATTTAAATAATAAATAACAATTTATCAAGAGCAAACAATTGTGGAAA
>CONTIG416 (59bp)
ACCTTTAATTAATATTTTTTTAAAAAATTACAAGTTGAATCTAAAAGATTCATCTAATA
>CONTIG417 (59bp)
AAATAATATTTATTCATTTTAATAAGATTAAATAAAAATAATTTTTTATTTTTATAAAA
>CONTIG418 (59bp)
TTTATTATAATATTCATAATTACTTTATTTAATATATGACCTAAATGAATTTCTCCATT
>CONTIG419 (59bp)
TAAAATACTATTATTCTTGTTTGCAAATTTTATTCCAATACTTGAAGGTAAACCAAACC
>CONTIG420 (59bp)
TTTATTATAATTATAATCTAAAATTATGTTTTGTTCCCTTTTATTGAGTAATAATATTA
>CONTIG421 (59bp)
AAAATTTTATTTTTAAGCTAAAACGTTTCTTCAATATTGTTATTTCTAAATATGTTTAA
>CONTIG422 (59bp)
AAAAACTTTTAATTATTTTTTTTTCAATTAAAATTCCAAGTAACAACATGCTAATTTTC
>CONTIG423 (59bp)
TGGAACACAATTCTTCGTTTCTTCCTATTTTAAAAATAATTAAGTTTAAAAATTGTAAA
>CONTIG424 (59bp)
AAAAATTATTTATTAAGTGTAATAATTTTTATATTTTTAGAAACTTAAATTTTTTTTTT
>CONTIG425 (59bp)
TTGAAATGATTCTGGAATTCCAGAATTAGCATCATTAATTCCTGTAATTATATTTTTGT
>CONTIG426 (59bp)
CTACATTATTAAATTTTACATTTTCTTTTATTAAAGATGTTAAACTTGTAATTGTATGC
>CONTIG427 (59bp)
AGCATTTTTAATGATCATTATAACAAAATTAAAACTTTTGATGATATTATTTATTATGA
>CONTIG428 (59bp)
CCTGCTTCTATAGAAGATAATAATATTAATTCTATAAAATTATGTTTTATTAATAAAAC
>CONTIG429 (59bp)
ATATAATAAAAAAAAATATTTAAAGTTAAAAAAAAATTCGTTTTGATTTCAAAAATCTA
>CONTIG430 (59bp)
GAAATTTTAATTTCAATAAATTTAAATATTTTTAGTATTAATAGAGAATCAATACAAGA
>CONTIG431 (59bp)
CTGTTTGCTCCCCACGCTTTCGTGTATTAGCGTCAGTATTAAACTAGAAAAATGCCTTC
>CONTIG432 (59bp)
AAAGAAAAAAAATTGTTAATTAACAACAATGATTTTTTATCTCCGTAAGGACTTTGTTC
>CONTIG433 (59bp)
ATTAATATTATTAATTCTTATTAAGGTTGGATTTAGAAAATTACCTAATATTGGTTTTT
>CONTIG434 (59bp)
TTAGTTTTAAAAATAATTTATTAGATAACAAATTTTCTGTATTATTGCCACAAAAAATA
>CONTIG435 (59bp)
ATTATTGTTATAATATTTTATTAAACAAAATAAATAAAAAATATTATTTTTATCCTATA
>CONTIG436 (59bp)
TTGTTTTAGAAATTTCTGTATTGGTAAATGAATATTTTTAATATAATTATTATTGGAGC
>CONTIG437 (59bp)
AAACATCATTTTTTTAAAATGCAAAGTATTGTTTAATAAGTTAGTTTTATTATAATAAA
>CONTIG438 (59bp)
TTATACAAACATAATTTGTATTAACATACTTTTTTAAATATTTTAATGAAAAATAATTT
>CONTIG439 (59bp)
TAACTTGCTTAGATATAATTTATTTTTCGAATAATTAAAATAATTATTACTTAAAAAAT
>CONTIG440 (58bp)
GAGCATAGATACGCTCCTAATCTTACTCTATTTGTATTTGATCTTCTTACATTATTTA
>CONTIG441 (58bp)
AAAATAATATTTTAGACTTTAACATTAATAATTTTTATTTAATATATAATTTTTATAA
>CONTIG442 (58bp)
ATTGATTGAGAATCATTTTTTAAACAAATATAAGTTTTAATATTTGAATCTTTTAAAT
>CONTIG443 (58bp)
ATTTGATAATTTTAATCTGCAAACTTTTCTTAAAGCAGAATTTGGTTTCTTGGGAGTT
>CONTIG444 (58bp)
TTATTTATATTATTAAATTATATAACTAGATTTACTTGGTTAAGAATTTTTTAATATT
>CONTIG445 (58bp)
AATTTTGAATCTTTATTTTTAATTACCAATAAAAACTTTATTAAAATTTTTTTAAAAT
>CONTIG446 (58bp)
GTTTATTAAACCAATTTCAAATATCTCAGAATTATAATAAGTATAAAAACAATATCTT
>CONTIG447 (58bp)
TATTCAATAGCTAATTCAATATCATATGCTATTTCTCCTTTTTCTTCTATTGTGTATA
>CONTIG448 (58bp)
ATTTATAAATAAAAAGCAAAATGTTAAAAACAAAATATTAAAAAAAATAAAATGAATT
>CONTIG449 (58bp)
CAAAATCACCTTTTAACACTTTAGTATATTGATTTAATATTGTATTTTGATTTGTTCT
>CONTIG450 (58bp)
AAAAAATTATTCTATATTTAATGAAAAAATTATTGTTAAAATAATTAATAAAAAAAAT
>CONTIG451 (58bp)
AAAAAATTTATATTTAATTAAAAAAACAAAAATTTTTTGTTTTCATTCTATAAAAACA
>CONTIG452 (58bp)
ATAGAATCAATAATAATTAATTTAACATATGTAGAATTTATTAGTTTTTGTGATATTT
>CONTIG453 (58bp)
AGATTTTTCATCCGGATAATCAGCTATCATTGCTTCTGTAGTAATTAATAATCCTCCA
>CONTIG454 (58bp)
AAAACTAATCAAATCACAATATTATTTGGAAACTTATCAATAAATGGTTGTATATCAA
>CONTIG455 (58bp)
AATATAGCTCAGCTGGCAGAGCAATAGTTTGTGGTACTATTGGTCGCGGGTTCAAATC
>CONTIG456 (58bp)
AATCTTTATCAATTTCTCAAAAAATGACTTTATGCAATATGTCTATTGAAGCTGGTTC
>CONTIG457 (58bp)
TACCAATGCATTTTGCATAAGCTAATGCTATGTTTAAAGAATTGTTATCATAGTCTTG
>CONTIG458 (58bp)
AAAATTAAATATTTGATATTTTAATACAAATAAATTTTTATTAAATTATAGAAATCAC
>CONTIG459 (58bp)
TTAATAATTTAAATATTTTTGATCTTAAAAAAATTATAATTTCCGTTTTATTATTAAA
>CONTIG460 (58bp)
TATATCAAATAGATTTTTTTCTAGAAAATTAAGATTTGGAATATTTAAATTTAAATTT
>CONTIG461 (58bp)
AATTGCTAGTTTATATTGAAAAAAAAATTTTCAAGCATCATAATATGTAAGTGCAGCT
>CONTIG462 (58bp)
AGATTCTGAACCAGTAGATCTTGAAAAAATTCCTCCAATTTTTGGAAAAATTTCAACA
>CONTIG463 (58bp)
ATTTAAAGATATTTTTTTTACTTTGTTAAAAAAATTATTTGTAAAATTACAATTTTTT
>CONTIG464 (58bp)
TAGCGCAAAACATACTAAATTTATATAAAATAAACCTTTTAAATAGCGATAATTTATA
>CONTIG465 (58bp)
AAAATATGTAATAAAATAAAAATAATCACTTCAAATTTAGAAAACTTAGGTATAATGT
>CONTIG466 (58bp)
CGACAAAATTCATATTTTAATAATTTTTATAAGAAAAATTTTATTCAATTTATCAGAA
>CONTIG467 (58bp)
TACTTTTTCTATTTTAAATAAAAAAATAAAATATAATCTTTTTTTGAACTTTATTAAT
>CONTIG468 (58bp)
TCCTTGAGTTACATCAAATGATACTCTTTGTCCATCTTGTAATTATTTAAAACCATCA
>CONTIG469 (58bp)
TTCCAAGTTTCGACTTTGAATAAAAAATTTTTCCGATAATATAATTATGTGATCTAAG
>CONTIG470 (58bp)
CAAATTTTAAAATTATTTTTATATGAAAAAATAAGTATTTTTTCTATCAACTTATTAA
>CONTIG471 (58bp)
CATCCCCTTTTTAAGAATAGATATTTTTAACTTGAATAAAAAATTTAATTTAAAAAGC
>CONTIG472 (58bp)
AAAAATTAAGATTGTGTCTTTTTTAGATTCTAGCGATAATTGTTCAAAAATCGATTTT
>CONTIG473 (58bp)
AAATTTTATTTTTAACAATTTTAATTCTACGTTTTTTTTTATAAAATAAAATTTTGCT
>CONTIG474 (58bp)
GATATAATTTCAAATTCAATAGAAGTTGTTGGCAAAGCACCAAATTTTGATGGTATTA
>CONTIG475 (58bp)
ATTATTGTTGTGCAAAATTTTCCAAATTAATTAAAAAAATAGGCTTAATTTCTACAAT
>CONTIG476 (58bp)
CATCAGGAACAGCATTGATTTTATATTCAAAACTTAAAAATATCAATTTAAATATATT
>CONTIG477 (58bp)
CCATTTTGAATATTTTTTAATATTTTTTTCATTTTTTTTTTAAAAGTATTGAAATAAT
>CONTIG478 (58bp)
ACAAAGCAGAAAAAAACTCTATTTCTTATTCTTTTAAAAAAATAAATTTATTTATTTC
>CONTIG479 (58bp)
TTATAAGAATTTCTTGAAAAAATGTAAAATTTTTAATACAAAAAATTGTAAAAAAAAA
>CONTIG480 (58bp)
TACTTGCGGAACTGCCGTTTTTATTAAAAATTTAAAAAGTATTCTGTTTAAAAATAAA
>CONTIG481 (58bp)
TTTGATTTAATTAAATTAGTTAAAATAAAATTTGTAGTACCATCTAATATGTTAATAC
>CONTIG482 (58bp)
AACAAAAAAAATGTGTTTTGTAAACATAAAATTATTTTTTACAAAAATATTAGCAAAT
>CONTIG483 (58bp)
TTTTATTAGAAATAAAATATTTTTGACAGGTCGTCCTGGCTTTATTAAAATAATTATA
>CONTIG484 (58bp)
TTAATTTTAGTTAACTTAAAAAACAAAAATATTTTAATAAATACAAAAGGAGGCGGAA
>CONTIG485 (58bp)
TTGAAAATTCATATAAATAAGAAAATATTTTTATAAATTTTTTCAACTTTTTATTTTT
>CONTIG486 (58bp)
GGAATAGAATTCATATATGAACTTGAATCCCTTTTATAATAATTTAGTTCTTGCCATT
>CONTIG487 (58bp)
AATATACAGAAATTAAACTTTAATTTGATTTTAACAACTTATTCTGAGTATATAAACA
>CONTIG488 (58bp)
TTATGAACGCATACATTTAATTGATTTTTAATATTTTTTAAAATTTTTAATACTAAAG
>CONTIG489 (57bp)
GAAAAATTTTAATTGATTCAAATTGTGTAAGAGTTTCAACTTTAAGATGTCACTCTC
>CONTIG490 (57bp)
AATTTCTTTAAAAGTTTTTTTTATTTTATGTTATTTAACAAATAGTTGACAACGTTT
>CONTIG491 (57bp)
CCCACCTTCTTCTCTTGATAATATATATACTTCACAAATAAAATTAGTATGAGGTTT
>CONTIG492 (57bp)
AATAAGCAATAAGTTAAAAATAATTATTACAGATTATAATGTTTATATTAACTTTTA
>CONTIG493 (57bp)
AATAGTTGTGTTGCTAGAAAAAATCCTCAAATAAAAATTTTAAAAGAATTAAAAAAA
>CONTIG494 (57bp)
TACATTCATATAAAAACAAATTTGGTTCAACACTAATTACTGTTCTTCCTGAAAAAT
>CONTIG495 (57bp)
TACCAAACGATAAAAATAAATCTTTAAACCCAAGTGGTATTAGAATTGGAACATCTT
>CONTIG496 (57bp)
AAAATTTTTTAAAATATAAAAACATTATTTTATATTTTTTATACTATAATATACCTA
>CONTIG497 (57bp)
AAAAAAATATATTAAAAAAAGTATATATTATATGTTACCAAAAAATTATAATAGATT
>CONTIG498 (57bp)
TAGAATAATATTCTAACATTAAATATTCAAAATTATGAATATTTGATACTCCTTCAT
>CONTIG499 (57bp)
CTTTTTTTTTTTTCAAAAAAAAATAAAATATAAAATACCAAATTTTAAAATTATTTT
>CONTIG500 (57bp)
CATTCAAGATTATTAGAAAGATCATCTAAAGTGAACAAATTTTCTGTAAATAAAAAA
>CONTIG501 (57bp)
AAAAATTTTAATGAAAAATTTGAAATTGGTACACAATATTTAATTGCTATTTTTTCC
>CONTIG502 (57bp)
TGTATTTATTAATTGTTTTTTAAATTTAATTTGGGGAATTTTTCATTTTATAATTTT
>CONTIG503 (57bp)
AACATCCATGGCTCTTATCTTTAGTTTTAAAATTCCATTTTTTAATAACCCCTGCAA
>CONTIG504 (57bp)
ATAAACCCAATATTAAATATTTATTTCTAATTTTTTTTTCTTTTTATAAAAATCCAT
>CONTIG505 (57bp)
TTATAAACCAGATATTAACCTATAAGATGAAAATTTTAAAATTTATAAAAAAGGTAT
>CONTIG506 (57bp)
TAAAATATTATAAAGTTTTTCATCATTAATATAATCTATAAAAAATTTTGAAATATA
>CONTIG507 (57bp)
CATGTAAAATTGGCCCAAATTGTTCATTTTTTAAATCATTAATATTATTAATTCTTA
>CONTIG508 (57bp)
AAGAAAACGAAGAATTATTAAAATTAAGTGATAATATTATAAATAAAAAATAACAAA
>CONTIG509 (57bp)
AAGTTCAGAAAACCGAGATGCTTGGGTTGAAGTAAAAGATAAAAAATTAGCGCCTCC
>CONTIG510 (57bp)
TAAAGTAAAATAACAATTTAAAAATATTTTTATTCCAGTTTTAGTTTTACAATCAAT
>CONTIG511 (57bp)
AAACAACTAATATCAAAAATAAATATTTTGTTACTATAAATTATCAAATTTATTAAA
>CONTIG512 (57bp)
CAAATGAAGTTTCATTAGATAAATCAATTATAGGTTGGAAAGGATTTGAATTAGAAT
>CONTIG513 (57bp)
AAATAAAATCCAATTGTGTTTTTTAAAAAAACAAAATCTTTTATTTTTTTAATATAT
>CONTIG514 (57bp)
AATGAATTTGAAACACCTTCTTTGTATTATTATTCTAGCAAAAATTTATATTTTAAT
>CONTIG515 (57bp)
AGATCAAAGCTTGAATCATTAGTTGAAGATTTAATTTTAAAATTTTTAAAACCATGT
>CONTIG516 (56bp)
CTGGAGTTTTTTCTTTTTTAACAGATTTTTTTCTCTTAAATTTTAATATTTGATTT
>CONTIG517 (56bp)
AAAATTTTATATCTAAAAATAAAAATTAAGAACACAAATTTTATAAATCATAACGA
>CONTIG518 (56bp)
ATTTTGTCTTTAACAATTAAAAAATTTTTTGTTTTATTAATATGATAAAAACAAGA
>CONTIG519 (56bp)
TTTTTTTTTTTTTAAGATTATTTTTAAAAAAAAATATATTTCAGATATTATTTGAA
>CONTIG520 (56bp)
ATATTTTTATATAAAATTTTTTTTCTTTTTTTTAGAATAAAATTTTTAAAATATAT
>CONTIG521 (56bp)
TTTAAGTATAAAAAAATCAAGCTGGAAAGACTGAAATAAATGAAATTTACTAAAAT
>CONTIG522 (56bp)
TTATATATTTATATTACATGTTTTATTAGTTAATATTTTGTATTTTTTAAATTATA
>CONTIG523 (56bp)
CATGCTAAATCTCCTTTTTTTGATTTAAAATTCGAAATTGGTTTTTTAATTTTAAT
>CONTIG524 (56bp)
TCAGTTTAAGAATACTTCTTTTTATTAAGAAATTATTTAAAAATAAAATTAAAATT
>CONTIG525 (56bp)
TTTCAACAATTAAAGTATCATTTTTTATTAAAATACTTATAAATTTTGGTATAAAT
>CONTIG526 (56bp)
TTAAAAATGGAACTTAAAAAAATAAAAATGGAAAAAATATATTATCAAGTAAATAA
>CONTIG527 (56bp)
TGTTAAGTTAGATTTTGTTAAATCTGCGCTTCCTCCAAATGTTTCATTTATTGCAT
>CONTIG528 (56bp)
TTAGAAAAACGTTTACAAAAAAAGTATTTTTCTTATTTCCCATTTTTCCTGGCATT
>CONTIG529 (56bp)
TAAACTTAAGAAAACCATTTACTTCATTATTTTCAACACAAAGCAGAGAACAAGAC
>CONTIG530 (56bp)
TTTTATATTGGTTCAATGTTTAGATATGAAAAACCTCAATTAACTAAAATTAGACA
>CONTIG531 (56bp)
AAAATTTTATCTAAAAAAATAAAAATGATGGATACAGCAGGTTTAAAAATTAGAAC
>CONTIG532 (56bp)
ATACCTCCTGCACATCCCAAAGTTTTTAATTTATTATTTATTAAAATATTATTGTT
>CONTIG533 (56bp)
ATTATCTAGATATGGATCGTAAATATATCCTTTCCATCCAATACTTGATCTTGGTT
>CONTIG534 (56bp)
TTTAACTTTAAATATTAAAATCGATATATTACCTTTTATAAAAAAATTATTGCTTT
>CONTIG535 (56bp)
TCAATAATTGAAATAGCAAATGTTCATGGAGAAACACAATTTGAAGTTTTATCAAC
>CONTIG536 (56bp)
TAAAATTCTTTTTGGTCTTCATCTGGAATTAAAATAACTAAGACATCAGATAATTT
>CONTIG537 (56bp)
ATCTAGTATATTCTTTCTAGTGTTTTTAAGTAAAAAAAAATATAGTAGTTTATTAT
>CONTIG538 (56bp)
AATAACAAAATTTATAAAATTTCACATGGTAAATCAGTTTTAATAGGAATTATTTT
>CONTIG539 (56bp)
ATGGAAAAATAATTGTTTTAAGTATAGTAGAACATATAGAAAGAACTGGAATTCAT
>CONTIG540 (56bp)
AACAATACATTAAAACTTATAAAGATATGAAAAAGTTATTTTTTGATATTTTAGAA
>CONTIG541 (56bp)
TTAAACAAAAAATAAAATTTTTAAATTTAATAACATATTGTAAAAAAAAAAATATA
>CONTIG542 (56bp)
CGATAATATAATTATGTGATCTAAGTCCAATTCCAACAACTGAAACTTTTGAAATA
>CONTIG543 (56bp)
TATTATTATTTTAACATTTTTTATTAATTTTAAAATTAAATTTTCTTTTGATAATT
>CONTIG544 (56bp)
AGAATTAAAATTTTTGTTTTAATAAACAAAATATCAGAAAAAATAATAGAAAAAAT
>CONTIG545 (56bp)
ATGGGATTAAAAAAATAATTAAAAAATTACCTATTGAAAATATTTTTGATCTAATA
>CONTIG546 (56bp)
CATTGGATTATAAATTCTATTTATGATTTTTCTCAACCTATAAACGATTATTTAAA
>CONTIG547 (56bp)
TGATATTATTAACTAAATTTTACCAATTTTTTTATGATTAAATGTTGTTGAATAAT
>CONTIG548 (56bp)
TTTTAAAGTTTTAACTTTAAAATTTAAAATTTTTTTGTTAAATTTTGAATATTTAA
>CONTIG549 (56bp)
TTGTTTTAAAAACGATGTATTAATTACAACATTATTTATTTTTAAAGGAACACTAG
>CONTIG550 (56bp)
TTAAATTACCATTATACATATTGCAATCTGAAATAAAATCGTTTTTAATAATTTTA
>CONTIG551 (56bp)
ATTCGATACTATCTAATTGTAATAAAAATATAAATAGTAATATTTGTAGTGGATGT
>CONTIG552 (56bp)
AATAAAAATTTAGATTTAATTTTAATAGAAAAAAAAGGAATAGATTTTTTTAATAA
>CONTIG553 (56bp)
AATTATTAAATAGTATCTGTGGAAGAATATTTCAAAGTAAATTTAATATTAAAAAA
>CONTIG554 (56bp)
AACTTTCTTTTTTTTTTTGGTATCATTTTTCCAACAATTACATGTTTTGAAAAAAC
>CONTIG555 (56bp)
AAAGTATCACCAGTAAATGAATTCTTTAATCCGATTAATACAATTATATCTCCTGC
>CONTIG556 (56bp)
AGCATCATATCCAAAATTATTTGAAGAAGATTTAATACTGTTTAATACAATTGATG
>CONTIG557 (56bp)
TTTTTTTTTTAAAAAATTTTTATTATAATGTAATTATTATTTTCTAAATTTTTTTT
>CONTIG558 (56bp)
TTTTTAAGTATAAAAAAATCAAGCTGGAACAACTGAAATAAATGAAATTTACTAAA
>CONTIG559 (56bp)
AAAAACACCTAATTTATCTTTATAAACTGATAGTAAAATTATACTATTTTTGTTTA
>CONTIG560 (56bp)
TTGCACAATGTTTTGCAAAAAGAAAATAAGGAATAATATTTCTTGAAAAATTTAGT
>CONTIG561 (56bp)
ATGTTTTTTTTATTTTTTAAATAACTAAAAATTAACAAAATTGAATTAAAAATATC
>CONTIG562 (56bp)
ATTTAAAAAAAATTCTTATTTTTTTATAAGAAATAGATAAATTCTTAATTTTAAAA
>CONTIG563 (56bp)
GTTTTGTATAAAAAATGGATATAAAAGTATTAATTTTTCAAATCATTTAAATCCTG
>CONTIG564 (56bp)
CTTATTTCAGATAAATAATTTATTATTTTTTTGTCTAAATTTACTTTTTTTAATTC
>CONTIG565 (56bp)
ATTCTAATAATTAATGCTTCAAATATTGAATTAATATTTATGTTATTTTTAATGGA
>CONTIG566 (56bp)
TATTTGATTTGTTTATTATTTCTATAGTTAGAAAACTACTTAATAATAACAAAATA
>CONTIG567 (56bp)
AGACGAAAAATTATATGTTGGTACACCATTTTTTTTAATTATTATTTCTTTATTAT
>CONTIG568 (56bp)
TACTAAATTTAATTGAGATATTGAATCTAAAACATTATTTACAATATTGTTTTCCA
>CONTIG569 (55bp)
ATTTTTATATACTTATTATAACTATTTATTATATAATTATAAACAAATTTGTTAT
>CONTIG570 (55bp)
TTGCAGCCACCTTACCTCTTAGACTAGTTCTTTATAACTAAAGTTTTAAATATTC
>CONTIG571 (55bp)
CCTAAACAAATTATAAATACATTAATTTCAAAATTTATTTTTTTATCAAATACTA
>CONTIG572 (55bp)
TTTGTTTTATAAATTTTAAATTTCATAATTTAACATTATTATTATTTAATAATAA
>CONTIG573 (55bp)
CAATAATAAGTTGGTTTTTAAGAAGATTTCCTTGGATAAAATATTTTTTTTGTTG
>CONTIG574 (55bp)
CCAATACTTATAATTTTATTTCCACCTCCAGATAATGTTCCCATTATCGAAGATG
>CONTIG575 (55bp)
TTTTTCTAATAAAAGATACTTAATTTTTTTATTATTAAATTTTGTTTTTTTCCAT
>CONTIG576 (55bp)
AAGCTACTGGTCATCCAATTGCAAAAATTTCAACTAAATTATCGATTGGATATAG
>CONTIG577 (55bp)
GTATCGCTTCCGATAGAACCAATAACACATCCACTACTAACTGAAATTGCAGATC
>CONTIG578 (55bp)
TATTAATCATTTAGAAAACAATATGTATAAAATTAATAATATTGCTTTATCAACT
>CONTIG579 (55bp)
ATAAACAATTTAGAAATTTTAAATGATATTAAGTAATGTATTTTACACATTTCTC
>CONTIG580 (55bp)
TTACTTAAAAATTAAATAATTAAATTTTTATGTCTCCAACAAAACATATAATTGT
>CONTIG581 (55bp)
TATTCAAAAAATTTTAAACTTAAGAAAACCATTTACTTCATTATTTTCAACACAA
>CONTIG582 (55bp)
CTGTCCCTCTTCCAGAAATTGAAAAAACATCTTCTATAGGCACTAAAAATGGTTT
>CONTIG583 (55bp)
ATACATTTCTTCTAATATGATATTTTACACCTGGTAAATCTTTAACACGGCCACC
>CONTIG584 (55bp)
ATTCTAATCCTTCAACAATAATGACTGATTATGAAATTGCAAATACTGTTTATAT
>CONTIG585 (55bp)
CAAAATAATAGCAAATGAAGCAAAACGTTTAAAAATTAAGATTGCGTCTTTTTTA
>CONTIG586 (55bp)
ACAGGCGCACCTGTATTAACGTTAATAAAAAATAATGATAAACAAAGTTCAGATT
>CONTIG587 (55bp)
TAATAATTTTAAAATTATTAATAATAATATATTTTTTTTCTTACAAAAATACTAA
>CONTIG588 (55bp)
TTTTTAAGTAAAAAATTAAATATTTTATATTTTAAGACAAATAAATTTTTATTAA
>CONTIG589 (55bp)
GCTAAATTTAAAGAATTAAGTGAATCAGAATACATTTATATTCCTATTAAACCTA
>CONTIG590 (55bp)
TATGATATAAATAATATTTATCACGTAATTAAAATAATTAAAATAAATAAAATAA
>CONTIG591 (55bp)
TTTTCTTTAGAATATTCTAATCTAGAATTTTCAACATGCGATATTAAATTTTATA
>CONTIG592 (55bp)
ACTTTAATATTGTTTTTTAATATAGAAATTGTAATTTTTTCACTAATTAATTGTG
>CONTIG593 (55bp)
GTCCGTATTGTCCACGTCCACCTGTTTGTTTTATATAATTTCCTTCTTGGATAAT
>CONTIG594 (55bp)
ATTTATTTTTTTTTAAAAAAATGAATTAATTATTTTCTTGTAATTCCAACTAAAC
>CONTIG595 (55bp)
TAATATTTGTTGAAAAAAAAATAAAACGTAGAGATGAGTTATTTTTAAATTATAA
>CONTIG596 (55bp)
TTATTGCTTGGAATAGGAGGAATTCTAGGAGTTAAAAACTCTATAAAATTTTGTT
>CONTIG597 (55bp)
AAATTAGGTTTCGATAACAAACTATTTTTTGTTGTGTAAGATTTAACAAAAATAG
>CONTIG598 (55bp)
TGTATTTAATAAAATTTTATTTCCTAATATTCGTTCTAATCCAATAAAATGTTTG
>CONTIG599 (55bp)
AAAAAACAGTTCGTTTTTTTAATTTTTAATTTAAAATCTTGCGTATCAATGTTTA
>CONTIG600 (55bp)
TAATTTAAAAACATGTTTTAATTTTAATAACGTTATTTTAATAAATTATTTGTTT
>CONTIG601 (55bp)
TTTTTGAGATATGCTTCCTACTTATATTGTTCCAGTAGTTATCATTTTAAACTTA
>CONTIG602 (55bp)
ATAAAAAAATTACATATCCATCAAAATATTTATTATATAGCCAGTTAACAAATAA
>CONTIG603 (54bp)
ACCAGAAGTATATTTGTCAAGTCTATTTAATATTCTAAAATTGGGAATTTTATC
>CONTIG604 (54bp)
TTAGATTTTGAAATTTCAAATTTAGTATTTGAATCAGATTACTATTTAGTTAAT
>CONTIG605 (54bp)
GTAAAATGAAACTTATAGGTGATATTAATATTGCTAAAACTTTCTTTATTTTCA
>CONTIG606 (54bp)
AGTTTCAATTGAAAATTCATAAACATTAAATAATGCTGGTACTGGACATGCAGG
>CONTIG607 (54bp)
CCTCCATTTCATCCATTTGGAAATATTTTATTTCTAATTTTTAAAAAACTTATT
>CONTIG608 (54bp)
TACTTTTTTTTTTAATACTAAATTCAAATATTTTAGTTGTACAATATTTTAAAA
>CONTIG609 (54bp)
TAATAACGTTATTTTAATAAATTATTTGTTTGAAAAAAAATAACATTAATTAAT
>CONTIG610 (54bp)
TAGATTTAAAAAAAAAATTTTAATCTTTAACAAATTAAAATTTGTTAATATAGC
>CONTIG611 (54bp)
AATAATTAAAGACATGAAATTAATAAAAAATAAAACTAGACAATTTAAAGTTAC
>CONTIG612 (54bp)
TCAATTAATATATCAGATATATTTAATTTATCTTTTATTAGAATTAAAATATTT
>CONTIG613 (54bp)
TGATGAAATTAATTTGTTAAAAAGTAAAGGTATAAAAAAAAAACATTGTTATAA
>CONTIG614 (54bp)
GTTTCTGTTTCATCTTTAAATGTTGTTAAAATTAAACCTACTCTAGAACTACCA
>CONTIG615 (54bp)
AACTTCATAAGATCGAGTTGCAGATTTTAATCCGTACTGAGAATAATTTTGAAA
>CONTIG616 (54bp)
CCAGAAGCATCTGATTTTGGTACAATTTGTTTAATTAATAATAATGATCTAATA
>CONTIG617 (54bp)
TATAAAAAGAATATTGGTTTAATTCTATATTATCCTTGTATTTATAACATTTGT
>CONTIG618 (54bp)
ATTAGTAATTTTTTTAAATTCGTTTATTAAAATTCTTATTTTACTTTTGACGTT
>CONTIG619 (54bp)
ATAATATTCATTATGATTATTTTGACCCTAGGTGTTTAAAAATTTCTTTAAATA
>CONTIG620 (54bp)
TATTTATTAAATTTCTTATATAAAAATTATTAGATATAATTCTAATAATTAATG
>CONTIG621 (54bp)
ATAAAATTTTTTTTAAAGAATTAAAAAATAAAATCATTAATTTATATGTAATAT
>CONTIG622 (54bp)
TAAATTTTTAGATAATTCAACTAAAGTAAAATAAGAATTTAAAAATATTTTTAT
>CONTIG623 (54bp)
TATTCCATAAGATAAAATAGTATTATTTATCATCTCGTAATTTTACATTTATAA
>CONTIG624 (54bp)
AGTGGTGATTTTTTAGCAATGCAAAGATTAAAAGAAGCAGCAGAAAAAGCTAAA
>CONTIG625 (54bp)
CGCAATTAATATTATTAAATATAAAAAAAAAAAATATATTTTTAAAGTTTTTAG
>CONTIG626 (54bp)
AAATTTTATTATTCAAAATTATGAAAATTATTTAAAAAAAAAAGTTAAAATAAA
>CONTIG627 (54bp)
TAATTTTTCAGTAAAAACACCTAATTTATCTTTATACACTGATAATAAAATTAT
>CONTIG628 (54bp)
TCCTATTATTTCTATCCATAAATTATTATAAAATATATCTATTACATAAGAATT
>CONTIG629 (54bp)
AAGAACAAAACTTATAATTTTTAAAATTTATTGTGGAAGATTTTAAGTTTAATA
>CONTIG630 (54bp)
ATTTTTTAGATTCTTTAAACAGAAACAATTTAAATAAACATATCAGTAATTTAC
>CONTIG631 (54bp)
TTAAAATATATGAATTTAAAAAGTAACAATTCTTCATATAAAATAAAAATCGAT
>CONTIG632 (54bp)
CTGGAATTCCAGAATTAGCATCATTAATTCCTTTAATTATATTGTTGTATAACT
>CONTIG633 (54bp)
AATATTTTGAAAATTTTGCATTCAATATTTTAAAAAAAAAAACAACAATAGAAT
>CONTIG634 (54bp)
ACTGGAATTAAAATTAATAAAAATATAAAAAAAAATTCAGTTGCTTTATCTCCT
>CONTIG635 (54bp)
AAAAAAAAATTGGATTTTTAATGTATTTTAAAGTATCAAAAGTAATTTCTTCTT
>CONTIG636 (53bp)
TTTTTAATTAAAAATCCATTTCTAAACTTTATTTCTGTACCATAATAAAATTT
>CONTIG637 (53bp)
CATTTTATTCCAATTATTAAAAATATTAAAACTTATTATTTATTTTTTTTTTA
>CONTIG638 (53bp)
TTACACCTGTTAAATCTTTAACACGGCCACCACGTACTAATACATTAGAATGT
>CONTIG639 (53bp)
CACCAAACAAATCTGTTTCTGTTCCATCTTTAAATGTTGTTAAAATTAAACCT
>CONTIG640 (53bp)
AATTTTAACAATAAAATTTTTTTTGTATAATTGTCCTAAGTTAGAAAATAACT
>CONTIG641 (53bp)
CCCCAAATTCTTTGTCTAGAAATACAACAATTTGATCTTAAATTTATCATATT
>CONTIG642 (53bp)
TATTTATTAATTTATAATAAAAATTCCACATCAAATTTTTTCCTAAAGTAATA
>CONTIG643 (53bp)
TTATTTTTACGAATTAATAAGTATAATGATTTTTTATAATTAAACTTAACGAT
>CONTIG644 (53bp)
TAAAAATAAATTATTATATAAAAAAATAATAATAAAAAAGTTTGTTTAAAAAA
>CONTIG645 (53bp)
TTACTAATACTAATTTGATTATATTTAATTTTAACAATAAAATTTTTTTTGTA
>CONTIG646 (53bp)
TTACACATCTTTTTATCCCTTAGTTTTTTATTCTTTAAAATTAATAGGACACA
>CONTIG647 (53bp)
TAATATTAACACAAAAATTATGCATTATTCAAATTTTGATATTTATTTTTAAA
>CONTIG648 (53bp)
AAAAATAATATTTGTTATAAATAAAAAGTTTTTAAAAGAGAAATATTCGTTTT
>CONTIG649 (53bp)
ATTTTTTTTCTTTTATTATTATTAGTTTATTATTCTTAAACAAATAATTTAAA
>CONTIG650 (53bp)
TGTGTAAAAAAACACCTTTTTTAACTGATCTACTCATTTTTTTTTTGTTTTAT
>CONTIG651 (53bp)
ATCAGATGTACCTCCTAAAAAATTTATTTTTGTTTTAATATTTAGAATTGTAT
>CONTIG652 (53bp)
ATAAAAAGACTATTTAGGATTTAATGCTGATAATTTTATTGAAATTGATGGTA
>CONTIG653 (53bp)
TTTCATCGTTTCCGTTTTTTTCAAAACTATTTAATAATAAAAAAAATAAATTT
>CONTIG654 (53bp)
ATTTAGAAAATCCTTCAATTTTTATATATGGTGGTACAATTTTACCAGGCAAA
>CONTIG655 (53bp)
CCGATATAAGCTCTTTTTTAAGCATATGGCTTTCTCATACTTTATCCATTAGG
>CONTIG656 (53bp)
AATTTGGACAAAATTTCGTTTTATTTCAAAAAAAAAAAAAAATAGATTTTTGT
>CONTIG657 (53bp)
TTTTATTTAATAAGAAAAAAATTCTTTTTTTTCCTCTTATTATGTTTATTTTT
>CONTIG658 (53bp)
TCATGTCTCCAGTAATTTGTTCAATATTATAAAAAGAATATTTAGGATTTAAT
>CONTIG659 (53bp)
AATAGTTTAATATCTTTAAAAATAGATATAGTAATAAACTCATCATCATACAT
>CONTIG660 (53bp)
AACTCCTTTTTATTTAAATTACATATTTAATTATTTCTTAAAACTTAATGAGC
>CONTIG661 (53bp)
TCATATTTATATTGATTATGCAACAATTGCTTTAATAAAAAATTTTAACAAAT
>CONTIG662 (53bp)
CAGGTTCACAAATAGGATGTACAATAATTACAAAAAATAAACAAAATATAATT
>CONTIG663 (53bp)
ATTTGAATATATAATAATTAACTTTAAGTAAATAAATTTTGTAAATCATCAAT
>CONTIG664 (53bp)
ATAAATGTTAACATTATAAAAATTAAAAAAAATATTTTTTTAGGATTATATGA
>CONTIG665 (53bp)
ATATATAGACATTGCAACTGCTGAAGTATCAGAACCTCCTCTACCTAATGTTG
>CONTIG666 (53bp)
ATATTTTTACTATTTAAGAGTTTAGCTTTTTTTTTTGCTAGAAGAATTTCTTC
>CONTIG667 (53bp)
ACCATTAAAATTTATTTTTTTTTTTAACGTATGTACATAACATCAGCACTACC
>CONTIG668 (53bp)
ATTTTTTTATACAAAAAAAATGTTTTTTGTAAACATAAAATTATTTTTTACAA
>CONTIG669 (53bp)
AAGAAATATTTTTAATGCTTATTTTTTTTAAAGTAAATAAATTTTTTACAAAA
>CONTIG670 (53bp)
TTTTATAATCAATAAATAAATTTTTTAAAAAACAGATACTGTTTATATCAATA
>CONTIG671 (53bp)
AGAAATAAGAATATTAAGAAAATTTTTTAACAAAAATGTTTTGTTAAAAAACG
>CONTIG672 (53bp)
CTATATGCACCTGGATTTATAATAATATAATTAAAGTTTATATCTTTTTGAAT
>CONTIG673 (53bp)
ATAAAAATTAAAACTACCTTTATTAATAAATATCATATAAAATCATAATTTAT
>CONTIG674 (53bp)
TTTTTTGAAAAAATTTTTAAAAAAAAAAAAATAATTTTATTTACTTCGTTTAA
>CONTIG675 (53bp)
GCTGCTATTAAATACAATATAAATCCAGAAAAATGGACAATTTTTAATATTAA
>CONTIG676 (52bp)
AATATTTTGAAATATAAAAAGAAAAATCAATTTATGATTTATTTTTTAAAAA
>CONTIG677 (52bp)
TTAATATGTTAATACAATAAAGAATTTTATTATTCAAATAATAATTGATAAT
>CONTIG678 (52bp)
CACTAGAACCTACTTTCGTATCTGTTCGAAATGTCTATCTTACAGTTAAGCA
>CONTIG679 (52bp)
AAATATTATAATTTTATTAATTTTTATATTAAGGTTAGTACTCATAATATGG
>CONTIG680 (52bp)
TGTAATTTTTATGCTATATCCTTCTAGTAAAAAACAGTTTGTTATTTTAATT
>CONTIG681 (52bp)
TCTTTAGATGATAAAATTGTATTTTCATAAAATAGTACTGACAAATCAATTT
>CONTIG682 (52bp)
TTTAGATAAGTAAGCTAAAGTATTTATTAATCCAATATTGTGTCCTTCTGGA
>CONTIG683 (52bp)
ATTTTTTATTGAATTTTTTTATAACTTCATAATTACAATAAGACATTATTCC
>CONTIG684 (52bp)
GTTGCAATAATACTTTATATAATGATTTTGAAAATAAAATATATATTACTTT
>CONTIG685 (52bp)
AAAAAAGGAAACACAAAATTTTTTATAAATTCCTTTTTTAAATTTTTAACAA
>CONTIG686 (52bp)
TAAAATTGTTTTTTTTTATTTTTTTTTTTTTAAATAAATATGATATTTTAAT
>CONTIG687 (52bp)
TTTTAAAGTTTTTTAAGTACAAAAAATAAAAATAAAATTTAATGATTTTAAT
>CONTIG688 (52bp)
AATAAAAATAATATATTTTTATGATTAGATTGATGTAAAAATACTTTTTTTT
>CONTIG689 (52bp)
TTAATATAATCTATAAAATATTTTGAAATATAAAAAGAAAAATGAATATATG
>CONTIG690 (52bp)
CTATAAAAAATTAGTTAAAAATTTAAATAAAAATATTTTACCTGGTACTAAT
>CONTIG691 (52bp)
AAATTATACTTCCTGTGATTCGATATTAATAGGTATGAGTAAAATTTATACA
>CONTIG692 (52bp)
AAGTTTTTACCTAATTGTAAAAATTTAATAAAAAAAAACAAAAGATTTAGTA
>CONTIG693 (52bp)
TTTTTTAAAAATTGATTTTATCATAAAATTTTAATATCAAAACCTAAAGATT
>CONTIG694 (52bp)
TTCTTTCATTTCTATTTCTGTTGTAGAACCTACTCTTATTACGGCAACACCT
>CONTIG695 (52bp)
ATCTTGTTCTTATAATTTTCATGATAATAATAATAATTATTCATTAAAAAAA
>CONTIG696 (52bp)
ATCTTTAAAACCATGTTAAATTGCGTTAAATGATGCTAAAATAAGTAAAAAT
>CONTIG697 (52bp)
ATAGATTTTTTGTGTCCATTTTTAAGAGGAGGAAAAATTGGTTTATTTGGTG
>CONTIG698 (52bp)
ACATTTAATTGATTTTTAATATTTTTAAAAATTTTTAATACTACAGAATTAT
>CONTIG699 (52bp)
TTTCTTTCTAAAACTATTTTGTATAATTCATTAATTTCAGAAGAAGCAAATT
>CONTIG700 (52bp)
TTTCAAAAACTATGAATTTTATTCAAAAAAAAAAATTTCAAAAGAAAATAAT
>CONTIG701 (52bp)
TTTTTAAAATATCAAAAACAGATATATTTTCTGTAATGTAAATAATGATTTT
>CONTIG702 (52bp)
TTTTTTTGTTTTTTTTTTTTTTTTTTAAATAAATAATTACTTAATGTAACTC
>CONTIG703 (52bp)
ATTTATAAGCGAAATTTACAAATTTTATAATTATAGTAAAAAATTAATTTGT
>CONTIG704 (52bp)
TAATTTTAAATAATTAAAATGTATAAAATATTGTTCTAAGTATATAGCTTTA
>CONTIG705 (52bp)
AAAAATATTTACAATAGTTCATTACAATAATAGTATTATTATAGTTGATGCT
>CONTIG706 (52bp)
TTAGACGCGTTTTATTTTTTTTATCTATATATTTTTTAAATAATGATTTTAT
>CONTIG707 (52bp)
AACACAATTTGAAGTTTGATCAACAAATGGAGATACTTTTTTAGGTGGAGAA
>CONTIG708 (52bp)
ACTTTATTTTTAAGTAAAAAATTAAATATTTTATATTTTAATATAAATAAAT
>CONTIG709 (52bp)
TTTAGATTTCCACTTGATCCAAATGGAAATTTACATTTTGGTCATACTTTTA
>CONTIG710 (52bp)
ACTTCTTTTATTCCAACAAATGTTATTTCTATAACTGATGGTCCAATTTTTT
>CONTIG711 (52bp)
AGTTTTGTTATATGGTTTTTGCACAATTGATATGAAAAAATATATTAGAAAT
>CONTIG712 (52bp)
AAAAAACAAAATCTTTTATTTTTTTAATATACATATTTTTTTTAATTTTTTT
>CONTIG713 (52bp)
TAAAATAGTTAAATCATGTATTTAGACTGGTTATAAAATAAATATAATTTAC
>CONTIG714 (52bp)
TTAATTTTCTAGATATTTTATATTTAAAAAATTATTCTATATTTAATGAAAA
>CONTIG715 (51bp)
TGATGAATCTATTTTAAGTTTAAAATTAAATAATTATTGGAAAAAAATATT
>CONTIG716 (51bp)
GGTAAAATTTTTATTTATAAAAAAATAGAATTTTCAAAAATAGATATTGAT
>CONTIG717 (51bp)
TCTGGTATTGAGAGTATAATTAATTTTAGAATAAATCAGGAATTAAATTTT
>CONTIG718 (51bp)
AGATATGCCTGAGTGGATTAAAGAATATGTTTGGAATACATATAAAACTAA
>CONTIG719 (51bp)
CTATTTCTAATTTTTATTTTTTTGTTCTGAATATAATTTCTCGTTGAAAAA
>CONTIG720 (51bp)
ATATTGTCATTGTTTTTGAATCAAATAAAAAAGTTTCAATTGAAAATTCAA
>CONTIG721 (51bp)
TATTACTTCAAATTTATTATCTTTAAAAATAATTTCAGATTCTGAACCAGA
>CONTIG722 (51bp)
CCATTGAGTTATAAACCCAATATTAAATATTTACTTCTAATTTTTTTTTTT
>CONTIG723 (51bp)
ACACGGTAGTTTTTTTTTATGTCTTAAACTTTTATAACTATTTAAATTTAA
>CONTIG724 (51bp)
ATTATATCACCAACTTTAACAAAATCACCTTTTAATACTTTAGTATATTGA
>CONTIG725 (51bp)
TGTGATTTTAAATAATTTTTTTTTTTAATAAAAAATTGTTTAAATTAATTT
>CONTIG726 (51bp)
AAGATTTGCAAAAGCAGGCGAATTTTTAGAAAGACGATATTTAAATGGAAA
>CONTIG727 (51bp)
TTTAGTAGTTCCTCCTGCTTTTTTTTGAGCCATAATTATAATATATTTTGT
>CONTIG728 (51bp)
TATTGTAAAAACATATCATCCATATTCAGATTGTGATTTTTCACAATCAAA
>CONTIG729 (51bp)
ATAATTTTTACATTAAATTTATTTTTTAATAACCATTTATAATAAAAATTT
>CONTIG730 (51bp)
CAAAAATTTCAACTAAATTATCGATTGGATATAGTTTATTAAACATATTTA
>CONTIG731 (51bp)
AATATTATGGTTATATTAAAAAAAAATATTTTAAATAATTTACTTAATTTT
>CONTIG732 (51bp)
ATTTCTTTTAAATAAAGATCAAATTTGTTTTTTTCAGCAATTTTTTCGTTG
>CONTIG733 (51bp)
AAAACAAATAATTAATTATAATTCTACCTAAGGTTGTATTTATCTTTTTTT
>CONTIG734 (51bp)
TAGCAACATTGTCATTATTAACTTTAATAATAATATTATCTAAATTTTTTA
>CONTIG735 (51bp)
AACTATTTTCGCTATGTGAAGGAAAGAGGATTCGAACCTCTGAAGCTTTCG
>CONTIG736 (51bp)
AATTGTTCTAAATAAAACACCATTGTTTATTTATAAATAAAAGTAAATAAA
>CONTIG737 (51bp)
TAGATTTTTTTCAAATTTATCATTTATAAATGAGATTAAAGATATAATTAA
>CONTIG738 (51bp)
GTACCAAAAAACTTTAATAATTTAAAGATTGTAAATAATACATCATTCAAA
>CONTIG739 (51bp)
TATCTCCCGAAATGAAATCTACAGGTGAAGAAATGAATTCTGGTTTAAGTA
>CONTIG740 (51bp)
TTTTTAATAATAATAATTTAAATATATTTGAATTACGCTATTTTAAAAAAA
>CONTIG741 (51bp)
TTTTAAACTTAAATTTCATACTTTTTTTTTTTAAAAAATTTTTATTATAAT
>CONTIG742 (51bp)
ACTAATAAATTTGGAGTTAATATTTCAAATATTAAACATTTTTTATTAATA
>CONTIG743 (51bp)
GTATAATTAAAAAAATTAAAAAACTAGGTAAAATAATATTTACAGAATTTA
>CONTIG744 (51bp)
GATATGGTTCATTAAAAACTTTAAACAATAATGCTAAAAATTTACTTTTAA
>CONTIG745 (51bp)
TAATTTAATTATCGTAGATTGTTCTTCATTATCTCCAAATGATATTGAAAT
>CONTIG746 (51bp)
TTCTGGTAAATTAATAAAAATATTGTTAGAAAATAATATTACACTAAAAAA
>CONTIG747 (51bp)
ATAATTACCATTTAGTGATTCCAATTTTATAATTTCATTATTAAATTTTAA
>CONTIG748 (51bp)
TATTAAAGTATTTTTTTTTAAAAACATTAATTAAATCAATATTAACAGTTT
>CONTIG749 (51bp)
AAAAATACCTGATTTTTTAAACTTACAAATTAATAAAAAATTATTTGATTT
>CONTIG750 (51bp)
TAAATTATCTAATTATTTATATTATTCTATTTACAAAGTTTTGTCATATGG
>CONTIG751 (51bp)
TAACAAATATAAATTTTATTCTTACTATTAATAATTTTAAAATTATTAATA
>CONTIG752 (50bp)
TGAACCTACTGCTGCAGCGCTTGCTTATGCTTTAGATAAAAAAAAAAATG
>CONTIG753 (50bp)
AATGTTCTTTTTTTTTTTCCATTTTCATATACTTATATTTGCCCAACTGA
>CONTIG754 (50bp)
AAATAAGGATTTCCTAGTCCTCCTTTTCCACCTTTTAAAATTTTAACAAA
>CONTIG755 (50bp)
GGTAATCTTTTTCTATTAACAAATATAAATTTTATTTTTACTATTAATAA
>CONTIG756 (50bp)
GAAATATTGAAAATGGTGTTTTTCTAAAATTACACTGGTTGTATAAATTT
>CONTIG757 (50bp)
CTATAGAATAAACACTCTTTTTCTTATAAGTATTAATATATTTATTTAAA
>CONTIG758 (50bp)
TTATAATCGAATTTTTTATAAAAGAAAATTTGGTATTTTTAATTTATACT
>CONTIG759 (50bp)
TAGTTTGTTTTTTTTTAATGATTTTTTATTTATAAAATATACATTTAAAA
>CONTIG760 (50bp)
ATAAAATTTAAACGTTTTAAAGTTATTAGAAATTGAATTAAATTTTTTTG
>CONTIG761 (50bp)
TATTTAATTAAAAAAACAAAAATTTTTTGTTTTTATTCTATAACAACAAA
>CONTIG762 (50bp)
ATGGACGAATAATTGAATACTTTTTTTTTATATTTTTTAAAAAAAAAAAT
>CONTIG763 (50bp)
AAATATAATATTATTCTAAAAAGATAGTTTATTTGCTTTAATTAAAATTT
>CONTIG764 (50bp)
TTGATTTTTTTCTTTTTAAAGTAATAAAAAAATTATTAAAAAAGAATAGT
>CONTIG765 (50bp)
CATAATAAATTCTAATTACTACAACATAATTTTTTGTAAAAATTTTTTTT
>CONTIG766 (50bp)
AATATTTTTCCAGGTTTTATAACTGAAACCCAGTCATAAATTGGACCTTT
>CONTIG767 (50bp)
CAAATTCATTGTCCAATTATTTTTGAAGAAATTTTAATAGCTTCATACGT
>CONTIG768 (50bp)
TTTTAATTTTTGACGATTATATAATATTAGTTTAACATTAAAACGGCTGA
>CONTIG769 (50bp)
TGGTTTAAAAATTTTTCAAATTGATGAACCAACAATTAAAGAATGTTTAC
>CONTIG770 (50bp)
AAAAAACAAAAATTTTAACTATTAAAATATTTCAAAATAAAGTTTTTTTA
>CONTIG771 (50bp)
CAATACTTGATTTTGGTTTTTCATAATAAATTCTAATTATTACAACATAA
>CONTIG772 (50bp)
CTTATAGAATAATAAAAAATTTTCTGTGTTAAATTTTTTTTTCCGTTATT
>CONTIG773 (50bp)
AAATTTTAAATATACACCTTTTATCTCGTTTTGCAATCTTAAAATTAATT
>CONTIG774 (50bp)
AAAAGATTTAAAAAATTCTATTTTTTTTGTATTAAACTATTTATTTAAAA
>CONTIG775 (50bp)
CGGTATAGGAATTTCTAGTCATCAAATAAATTGTTTTAAAAATATTATAA
>CONTIG776 (50bp)
TCTCCTACCACAAAAATAATGTTCATAATTTCGGTATCTATTTGAGCCCC
>CONTIG777 (50bp)
TGTATAGACCTGGACCTATACAATCAGGATTTATTAATGATTTCATTAAT
>CONTIG778 (50bp)
TCCTCCGCATAAAACAGATTGTTGACCAAACAAATCTGTTTCTGTTTCAT
>CONTIG779 (50bp)
TGAAATATAAGAAAAAACACTTTTTACACACATATAAGATATTATTTCTG
>CONTIG780 (50bp)
CCTATTTTTAAAGACCAAGAAAAATTTATTTATTTTAGAAAAACAAGAAT
>CONTIG781 (50bp)
TAAAATTAGAAAATAAAATTGTAAAAATCCAATTAAATATTTTCACACCT
>CONTIG782 (50bp)
ATATAATTATAATAAAACTCTTTTTGGTCTTCATCTGGAATTAAAATAAC
>CONTIG783 (50bp)
TTTTTTATAACGTATGTACATAACATCAGCACTACCTTATATATATAATA
>CONTIG784 (50bp)
GATAAAACTTAATAGCTTAAATTTTCTTAAAAAAAAAGATTTATGTTTTC
>CONTIG785 (50bp)
TAATTTTATGTGTTCTACCTGTTTTAATACATATGTTTAAAAAAGAAAAA
>CONTIG786 (50bp)
CTTTTTTAAAATACACGTAACCAATATTTTCCGAATATATATAAAAATTA
>CONTIG787 (50bp)
AAATGACAAATATTAGTTTTGAAATATTTCCTTGTAATAACATTAAAGAC
>CONTIG788 (50bp)
CAAATCAAAATATTATTTGGAAACTTATCAATAAATGGTTGTAAATCAAA
>CONTIG789 (50bp)
AAAACATCTAAATAATTTTCAATAATATTAAATTGATCAACATCGTTAAT
>CONTIG790 (50bp)
TTTAGGAGAAATAGATGTTGTATTTAATTCTATATGTTATAATAAATTTA
>CONTIG791 (50bp)
TTATATAGTTTTTATAAATTAAATATTTTAAATAAAATACCTACAAAAAA
>CONTIG792 (50bp)
TTTTAAAAAAATAAAAATTTTTCCATATTTAAAAAAAAACTGATAAAAAA
>CONTIG793 (50bp)
TTTTTTCTTTTATATTAATTAACCAATAAAAAAAATTTTTTTTAATTAAT
>CONTIG794 (50bp)
AATAGGAAAATAAATATTGTGGAACAAAAGTAGTATTTTTTCTACTTCAG
>CONTIG795 (50bp)
TATTATGATATTAATAGTCCAATACAAGGAATAATTCATATTATAGCATC
>CONTIG796 (50bp)
ATTATTTCTATGGGAGGATGTGACAAAAATATTCCTGGTTGTATTATTGG
>CONTIG797 (50bp)
AGATAAAGACATCGTATTTATTAGTATTTGATAAGAATTATCATTAAAAA
>CONTIG798 (50bp)
CAAAATCAGACAATCGCTTATTTTTAAGATTTGATAATCCAATATTTAGA
>CONTIG799 (50bp)
TTAAAATTAAATAGCATTAAAAATTATAATTGTATAAATATAAACGATTT
>CONTIG800 (50bp)
TTTAATAGTATTGTTAATCTTATTAAATTTAAAATAAAGTAAAAATGTTA
>CONTIG801 (50bp)
AATTTTTTATTTATTGATAGATCACTGGTATTTAAAAATAAAACACATAT
>CONTIG802 (50bp)
GAAAATAAAATATATATTACTCTGCAAAATCATAATTTTAATATTAAAAA
>CONTIG803 (50bp)
ATTTATAAACATTAAAAACAAAGATAAATTAAATAAATTTTTTCTTTATT
>CONTIG804 (50bp)
AGTTAGCTATTTTTTTATATAAAAATTTTTTTTAAGAAAATTATTAGTTC
>CONTIG805 (50bp)
ACTGGTTTTTTAAAAATAAATATGTTATGAAACCTTTAATTAATATTTTT
>CONTIG806 (50bp)
CATCTGGAATTAAAATAACTAAGATATCAGATAATTTAACTGCCTCTTCA
>CONTIG807 (49bp)
AAGGTGCTTTTTCATCTTTCCCTCACGGTACTAGTTCACTATCGGTTAA
>CONTIG808 (49bp)
TTGATAACAATATATTTTTTTTTTCAAATACATGAAAAAAAAAAATAAA
>CONTIG809 (49bp)
ACATTGACCAACAGAACCTAATGTTCTGTATGAAAGAGAACATGCATGG
>CONTIG810 (49bp)
AATTTATACAAACGAATTTAAAAATAATTCTAAATTTTATCATTATATT
>CONTIG811 (49bp)
ACTTATATTTATTTAATTTTATTATTGAGAGTATTTCTTATATAATGAA
>CONTIG812 (49bp)
ATTGTCGCTTCAACTGCCTCAGATAGTGCAGCGGAGCAGTATATTGCTC
>CONTIG813 (49bp)
AAAAAATATTTTCAGTAATATTAAAATATTATATTATTTTTATGACTAT
>CONTIG814 (49bp)
TTTTTAACTTTAAATATTAAAATCGATATATTACTTTTTATAATAAAAT
>CONTIG815 (49bp)
CAAAACAAAATACAAAAAAAAAATGTTTTTTTAACTTTAAATATTAAAA
>CONTIG816 (49bp)
TGATATGTAATACTGGATTAAATAATTTTACACCAATATTATTCGATGA
>CONTIG817 (49bp)
TATAAAAATGCTTCTATATTAAATGTAAAAACATACTTAAAACCTTTTT
>CONTIG818 (49bp)
AAAAACTTTTTTTAAAAAAAGATAATTGGAATTCTCCAATTTTAGGAGC
>CONTIG819 (49bp)
CGTTTTTTTTAAAAAAAAAAATAAATCAAATTCAAATATTTTGCAATAG
>CONTIG820 (49bp)
TCAAATGTAATTAGATCACTAATTATATTATTATCAAAGATAAATTATC
>CONTIG821 (49bp)
TATCCTATTATACCTAATTTATCTAAAATTTCTTGGTTTAAAAATGGTT
>CONTIG822 (49bp)
TAAGTAGAAAAATAATTTTTTAATATTTTAGAACTACAACTTAATAATT
>CONTIG823 (49bp)
TATTGTAGTATTGTTTGAAATAAGTAATACAAGAAAAAATTATTTTATT
>CONTIG824 (49bp)
TTTAACAATGTTAAAGTAAGTTTTAACTATAATATACTTTATATAATAA
>CONTIG825 (49bp)
TTTAATAAAACATTCATTAATAATTTTAATGTATCAAAATACATGTTTT
>CONTIG826 (49bp)
CCAGCCACAGGTTCCCCTACAGCTGCCTTGTTACGACTTCACCCCAGTT
>CONTIG827 (49bp)
AAGTTTTTGTTTAGATTAAATAAGTTTTTGTTTGAAAATGATATAAAAA
>CONTIG828 (49bp)
TATGGTTAACAAGTTCTAATATAGGAAAAGAACATTTTAATTTAGATTT
>CONTIG829 (49bp)
TTTTAATTGAAACGATTCATATCGTGGTTCTATTGTATTGTAATTATAA
>CONTIG830 (49bp)
TGTATGTTTTTATATATTCTAAAAAAACTATTACAATAAAACAGTTTTA
>CONTIG831 (49bp)
CATATTTAAACCATTTCTTAAATTATCAAAAATCGGATCTTTGATAATA
>CONTIG832 (49bp)
TTCTCAATTATCTTATTTAGAAGGTTGTACCGCTTCAATAAAACAAAAA
>CONTIG833 (49bp)
TTGTTTTAATATCACCTATTAAATATTGGTTTAAATTTTTAATCAAAAA
>CONTIG834 (49bp)
AAATAGATAAATTTTTAATTTTAAAATTATCATAATAAATATTTTTAAA
>CONTIG835 (49bp)
TAAAAGTTAATTTTTAAAACATGAATATTTTAAGAGGATTTAATGATTT
>CONTIG836 (49bp)
GTATTTAATCTTTATCAATTTCTGAAAAAATGACTTTATGCAATATGTC
>CONTIG837 (49bp)
TGTTACCATTTAGCATTTACTTTTTGTATTATTCATAGTAGCACGTGTG
>CONTIG838 (49bp)
AAAAAGATGATTAAAATTAATATTCAATCTGGAAAAGGAGGATTAAAAT
>CONTIG839 (49bp)
TTTACAAAAATGTATTTTTTTAAAAAAAAATTTTTAATTAAAATGTTGT
>CONTIG840 (49bp)
TCTTTATCAATTTCTGAAAAAATGACTTTATGCAATATGTTTATTGAAG
>CONTIG841 (49bp)
ACGATATTATATTTAATTTTACTAATTTTATGTATAAAAAAATTTGAAT
>CONTIG842 (49bp)
TTCTTCTATAATACCATACAAAAATGATATAAAAAAAAGAATATTAAAA
>CONTIG843 (49bp)
TTGCATTTATTTTAAAACAATTTACATTAAAAAACATTTTTTTATAAAA
>CONTIG844 (49bp)
TGATTCTAATCGCTAAAAATTATTTAGGATATTTAAAATTAATTAAAAT
>CONTIG845 (49bp)
ATATTTAAAAAGAAAGTTTTAAAGATTTTTCAAGAACAATAGGTTTAAA
>CONTIG846 (49bp)
AACAATAAAATACTTTTTTTTTTAATTAATATATTAATTAAATTGTTAT
>CONTIG847 (49bp)
ATCAATAATTTTATTATTTATCTTTGTAATTTTTATTTGCAAATTAATA
>CONTIG848 (49bp)
TTCAAATAATATTCGTAAATTAATTGATGAAATTCCTTGCTTATTAATT
>CONTIG849 (49bp)
TTCTTTAAAATTTTTTTTTATTTTATGTTATTTAACAAACAGTTGACAA
>CONTIG850 (49bp)
GTATCTTCTGGTTTAAGTATTTATGATACTATACAATTTTTGATATCAG
>CONTIG851 (49bp)
TTAATACTAAAGAATTATTAAATATATTTAAATTAATATTTTTGTCATA
>CONTIG852 (49bp)
ATATTGAATCTAAAACATTATTTACAATATTGTTTTCCATAATGTTTAA
>CONTIG853 (49bp)
GTATAACTTTTAATATAATAGTAACTAAAAGAAAATTTATATTTAATGA
>CONTIG854 (49bp)
AAGTTTAGATAAACAAAAATTTCTTATTGCACTATACTTATAATTTGAA
>CONTIG855 (48bp)
TATAGAAAATAATTTACTAAATGTATTGCATTTTTTAATTTTAATTTT
>CONTIG856 (48bp)
ATTTTATTTTAATTTTTATTTATAAAAACAAAAAATTATTCAATATAC
>CONTIG857 (48bp)
TATTCTAACTTAAATATCGATTTTAATAAAACTATTACATTAATAAAT
>CONTIG858 (48bp)
TTAATAAATTTTTATAATTATAATTTAAAATGTTTGTATTTAAATAAC
>CONTIG859 (48bp)
TCACTCGGTTCACTTTTATTTTTAAATAAATACTAGTTTTTAAGTAGT
>CONTIG860 (48bp)
TTTAAAAGTTATACTAAAAAAGAAAAAATATCTTTTTTAAAAATAGAA
>CONTIG861 (48bp)
AAAGCTAAAAAAAGTTAATATTTTTTTTAAATCGATTATAAACAAACT
>CONTIG862 (48bp)
TTAGCATAAGCTTGAACAACTATTCCAAAACCTTCCCAATTTTTACAT
>CONTIG863 (48bp)
AGAACCGGGCAATTTTTTAAATTTTAAACATATTAGTATAAAAATTTT
>CONTIG864 (48bp)
TATATTATAATATACCTATTTTAAGCGTATGTTTAGGACATCACATTA
>CONTIG865 (48bp)
CAAAAAACATTTTTTTAACATAATCATATGTTAATGTAAAATTATTTA
>CONTIG866 (48bp)
GTCTTTTGCAATAAGATATAAAAAAGAAAGATTTACTTCTTTTTTAAA
>CONTIG867 (48bp)
TTTTCATTTTTTTTTTAAAATTAATTTTCGGTTTACATTTAATAAAAT
>CONTIG868 (48bp)
ATGTCACTATGATTTAGTGTTATCAAAATTAATTGATTTTCAAAAAAA
>CONTIG869 (48bp)
AACCCTTTTTTTTTTATTTAAAAAAAAAAAATATGCATATAATTAATG
>CONTIG870 (48bp)
CGCTTGTAATTAATAATATACTATTATTTTTGTTTGCAAATTTTATTC
>CONTIG871 (48bp)
AAATATTTGTCTAAAAATTTTTTAATAATACCGTCTTTACTATCTAAA
>CONTIG872 (48bp)
CAGGATGTTATTTTGAAAATCACAAATTTTCTACAAATTTTAAAAACG
>CONTIG873 (48bp)
AATAATTGTAAAAATACTTGTATTATAAAATTAATTAATACTAACTTT
>CONTIG874 (48bp)
TTATATTGTAAGTATAACAAAAAAGCTATTTTATTCGGAATAGAAAAA
>CONTIG875 (48bp)
TGATATGCGAGGAGCTATTACTTCTTCGATTTTATCTTTATTATCTGG
>CONTIG876 (48bp)
TATTTAACTTTAAAAATAATTCTTCAACATTATAATTCTTGTTATTGA
>CONTIG877 (48bp)
AATATTATTTCAGGACCAGGAAGAGGATCTGGTTCTTCTTCGATACTA
>CONTIG878 (48bp)
ACAAAAACAACAAACACCTCTATTTGAAATACAAAATAATCCACTTCT
>CONTIG879 (48bp)
GCTTAAAGTATATTCTTTATAAAAACCAAATTTTGAAATAATACTGAT
>CONTIG880 (48bp)
CACGAAATTGATCATTTAAATAGTAAACTAATATTAGATTATTTTAAA
>CONTIG881 (48bp)
TTATCTTTAAATAATAAAACAAAAAAACTTTTAATTATTTTTTTTTCA
>CONTIG882 (48bp)
ATAAAAAAATAAAATATAATTTTTTTTTGAACTTTATTAATTACTATA
>CONTIG883 (48bp)
ATATCCATTATTGTATTTTATGTTTATTAAACAATTGTTTTTAATAAT
>CONTIG884 (48bp)
AAATAATATTAAACAATTTCAATTATTAAAATTTGGTTCTGTTACAAT
>CONTIG885 (48bp)
TAAAATTGTTAGAAAGATATTAATATATAAAATTTTAGATAATTCAAA
>CONTIG886 (48bp)
CATAATTTCGGTATCTATTTTAGCCCCGTTAAATTTTTTGTATCAGTT
>CONTIG887 (48bp)
GTCAGCATTCGCACTTGTGATATCTCTAAATTATTTTAGAACAATTCT
>CONTIG888 (48bp)
TTTTTTGAATTATTAAATATAATATGTCCAGGTTCAATTTTACCAGAA
>CONTIG889 (48bp)
TGATATAATTAAATTTTACTTAATTAAAAAAAATGAAAATTATCTATT
>CONTIG890 (48bp)
TGTTATGATTCAGTAAAAACTTAATCAATACAGTTCTTTACCAATAAT
>CONTIG891 (48bp)
TTTTAAATAGAAAAAATTATTTTTGATAATTTTTTCATAAAAAAAAAT
>CONTIG892 (48bp)
TAAAAAAAAAATATAATTATTTAATAAAACATAATTTTAATTTTGATT
>CONTIG893 (48bp)
AAAAATATATCCAATATGGAGATTCATATTGAATCCAACTTTTTTTTT
>CONTIG894 (48bp)
ATCAATAAATTTTGCGAAAGATATTTACAAAATGTATGTTAAACTTTT
>CONTIG895 (48bp)
TTTACCTATTGTTAAATTTAAATACTATAATTCTTTAATATACAAATT
>CONTIG896 (48bp)
TCTACCAAGTATATCATCTGATGTAATCGTTAACATTTCTTTTAACAA
>CONTIG897 (48bp)
AATTCAAGAACCTTATAGAATTTATACTTCAAAATCAGACAATCGCTT
>CONTIG898 (48bp)
AAAACATAAACATTTATTAAATGCTTATAAACAAATTTATAATTTTAA
>CONTIG899 (48bp)
TGATTTCTCTACGCCTACCTAGAGTTAAGCTTGCTAAAATAAATAAGT
>CONTIG900 (48bp)
ATTACTATTATTATAAATGTTGTAAAACATTGTATTATTTATATTTTT
>CONTIG901 (48bp)
ATAATTTTTATTTTTACTATTAATAATTTTAAAATTATTAATAATAAT
>CONTIG902 (48bp)
TAAAAAAGAATTTTTAAATACTTTATCTATTAAGTTATTTTTTAATAC
>CONTIG903 (48bp)
CATGTGTTCGAATCCCATTGTCTCTTAAAAATTTTGTTTAATTAGAAA
>CONTIG904 (48bp)
TTTTTATACAAATATTAAAACTATTTCCTTATTATTTTTTAAAAAGAA
>CONTIG905 (48bp)
AGATGGAAACGCAACATCTTTAGAGCCAATTCGTAATGGTATTATTAT
>CONTIG906 (48bp)
AAATAACTATTGTTAAAAATATTATAAAGTTTTTCATCATTAATATAA
>CONTIG907 (48bp)
AAAAACAAATTTATTTATTTGATTTACCTGGTCACTCATTATTTAGTA
>CONTIG908 (48bp)
TGATGTTAGTCAAATATCTTATTTAATTTCTAATAAAATATATGATAA
>CONTIG909 (48bp)
ATTCCAAGAGGCATACCGGGATGTCCAGAATTAGCTTTTGAAACTGAT
>CONTIG910 (48bp)
CTTCTGAAATTATATTTAATATAATTAATAATATGAATGAATTTATTA
>CONTIG911 (48bp)
TGTAAATTAAATTGATTTAACGAATATAAATTATTTTTTTTGTGTTTT
>CONTIG912 (47bp)
AATATTAAAATATACAACGCACCCAAAAAAACATTATTTATTAGAGA
>CONTIG913 (47bp)
GTTAAAATTTAATTTTTAAATTACTTTATTAAACTTGATATTTCTAA
>CONTIG914 (47bp)
GTGGTTCGATTCCACCTTTAGCTATGAGATATATCCAAGTGGTAAGG
>CONTIG915 (47bp)
TTATACCATTTTTTAAATTATCAAAAATCGGATCTTTGATAATATCT
>CONTIG916 (47bp)
TTTACTGTTTTTGTTTTTTTATTTTTTTTTTCTATAATATTTAATTT
>CONTIG917 (47bp)
GTTAAATATTAATTAAAAGTTTCATAAAAAAAGTTATTTATGTTTTT
>CONTIG918 (47bp)
AATATCAAATTTAGTATTTGAATCAGATTACTTTTTAGTTAATAATT
>CONTIG919 (47bp)
TTAAAGTATTTTAAAATAACGTTATTTATGTAAGTTATTAGTATTCC
>CONTIG920 (47bp)
AAAGATGTGTATTAACATCTAATAGAAATTTTGTAGGAAGACAAGGA
>CONTIG921 (47bp)
TGTTTTTTAAAAAAAAATATTAAAAACCTTGCAAGTTTTAAACTAGG
>CONTIG922 (47bp)
TAATTTCATAGAAATTCCAACAACTTTATTATGTCAAATTGTTTCTA
>CONTIG923 (47bp)
TTTTTTTTAATTTTTTTAATTTTAACTAAAACATCATTATTTACTCT
>CONTIG924 (47bp)
AATCAAAAAAATAATATTAATTTTATTGGTAATTCTGGTACAAGTAT
>CONTIG925 (47bp)
ATTTTAAAAGAAAAATTAATTTATTTTTAAAATCACTTATTTTTTTT
>CONTIG926 (47bp)
GTTTGATAATAAAACATCAGAATATTTCTTCGTAAAAAAACATGATA
>CONTIG927 (47bp)
AATAGAATTAAAAGATGAATTTTATTTAAAAAATATAGATGGTCTAT
>CONTIG928 (47bp)
TTATTCTTAAAAAAAAAAAAAGATTTAATAGTAACTTTAATCATTTA
>CONTIG929 (47bp)
AATCAAAAAATGGGATTTTTTTCAATGTATAATGTTTTTTTTTCATG
>CONTIG930 (47bp)
ATTTTTAGATCTATATTCTAAAATGTTAAAATAACGAGTAATTTATC
>CONTIG931 (47bp)
ATTGTAAACAAACTTTAAAATACGGAATTATAGATAATGTTTTTTAT
>CONTIG932 (47bp)
TTATAAGAATTAATAAATAATAAAATTATTTATTTTTCTAAATTATT
>CONTIG933 (47bp)
AATTAAATCAATTAGTTCATATTTGTGATTATAAATAATCTTTTTTT
>CONTIG934 (47bp)
TATTTATAAAATATACATTTAAACAACTATTAAAACTTTTAAAAATA
>CONTIG935 (47bp)
AAAATCTAATATTTTAAAAGCAGGTTCTTTAACTGCATTTCCTCTAA
>CONTIG936 (47bp)
TTAAATTTTTTAATTAAAATATTATTCTTGATCATTAATTTAAACAT
>CONTIG937 (47bp)
TCACTATAAGGAGTAATTATTGCAACAAAATTATACATTTTATATAT
>CONTIG938 (47bp)
CAAGCACCTCTATTTGAAATACAAAATAATACACTTCTAACAAAAAT
>CONTIG939 (47bp)
CTCACATATAAATCTGAAGAAACTTTAGTTAAAGCAGCAGTCAAAGT
>CONTIG940 (47bp)
AATGGTGAATTTGATAATATTAATGAAAAAAATTTTTATATGACAGG
>CONTIG941 (47bp)
ATTTTATTTCTTTTTGAAAAAAAAAAAAAAGAAAAATTTTTTGTATT
>CONTIG942 (47bp)
AAAAAAAAATAACAAAAGTAGAATTCCATGGTATTCAAATTCAAATT
>CONTIG943 (47bp)
AAAAAATTACCTTTAGAATTTAACAACGAAATTGAGAATTAAATATC
>CONTIG944 (47bp)
AATATAATTATTAATAAAAAAATTTTTATCAAAAATGTTTATTTTAG
>CONTIG945 (47bp)
TATTTTAAGTAAATCTTGTTTTTATAATTTTCATGATAATAATCATA
>CONTIG946 (47bp)
AATTAACTTTAGAATAAATCAGGAATTAAATTTTATAAAAAAAAATT
>CONTIG947 (47bp)
ATATTATGTATATACAGATGTTTGGGAATCAATGAATAATAAACATG
>CONTIG948 (47bp)
TAGAGTATGCTTTATTTTGATCTATTATTTTATTTATTTCATCTTTT
>CONTIG949 (47bp)
GCAATTTTTATATTTAAAATATTACTTTTGTTTAAAACTAAATTAAT
>CONTIG950 (47bp)
TAGGCAATTTACACAATTTGGTTTTGAATTGTACAATAATAAAATTT
>CONTIG951 (47bp)
ACTCTTTGTGCTTCTTTTATTATAGAAAAAGCAAACGTAGTTTTTCC
>CONTIG952 (47bp)
ATTTAATTTATTATAATATAATGGAGTTATAGTAAGAAAATCACAAC
>CONTIG953 (47bp)
CATACCGTAATAATACTTAAATTACTTATGTTACAATCTACTTCCAT
>CONTIG954 (47bp)
AGATTTTGAAATAATAGCTAATTTTCCTTTTTTAATAATACTACTAG
>CONTIG955 (47bp)
TTTCAAAAATTAAAAATATTAAAAGGATAAATTTTCTATCATCCAAC
>CONTIG956 (46bp)
TAAAATTAAACCTATTCTAGAACTTCCAATGCATTTTGCATAAGCT
>CONTIG957 (46bp)
CTATAAATATTTTCTAATATTGTTTTTTTACAAACATAATAAGGTA
>CONTIG958 (46bp)
ACTATCAAAATTTCTTATTTTTTTTTTATTATTATTTACATTATGA
>CONTIG959 (46bp)
ATTCCTAATTCTTCACTTATTAAAGTTGAACCAGTTAAAATAGGAA
>CONTIG960 (46bp)
TTTAAATAATGCATCTAATTTAAACAATTTAATAAATATGAAATTT
>CONTIG961 (46bp)
TACTAAAATTGATAAATTAAAACAATTACTTGAAATATATTTATAT
>CONTIG962 (46bp)
ATATGTAGAAAATGTTTTGTTCAATTTACAATTAAATTTATTTTTT
>CONTIG963 (46bp)
TTAGTTTATTTAGTATAGGAGATCAACTTCTTTGGGGAGCTGCCGA
>CONTIG964 (46bp)
AAACACAATTTGTACAAATATTTTTATTTAAAACATAATTTTTTTT
>CONTIG965 (46bp)
TTTTTCTTAAAAAAATACTAAACAGATTTTTCATCCGAATAATAAG
>CONTIG966 (46bp)
TCATTAATGTTGGTATAAATGATTTATTTTTTGAAAAAACTTTTAA
>CONTIG967 (46bp)
ATCTTCTTCTTTTAGAAAATAATAATTAATGGTATCACTTTTATAC
>CONTIG968 (46bp)
TTCAATGTTTATGCTATATTGTTTACAAGTAGATAATTCTTGCATC
>CONTIG969 (46bp)
TTATAATACATTTTTTAGAAAAAAAATTTTTAATTAAAATTGTCAG
>CONTIG970 (46bp)
GCTAAATCGTATTAAATTAACAAATAAAGTTTTATTATTACATAAT
>CONTIG971 (46bp)
ATCTTAAATTAAGTTTTTGTTTAGATTAAATAAGTTTTTGTTTGAA
>CONTIG972 (46bp)
TATATTATCGAATAAAACCAATAACAATTTATTTTCATGTAGTAAA
>CONTIG973 (46bp)
GCAAGCAAGTGATGTGGAAATTCATGCGCGTGAAATGATAATTATT
>CONTIG974 (46bp)
ATTGTTATCGCAATTAAATCGTTTAATCGAACATAGAATTGAAAAA
>CONTIG975 (46bp)
TTCAATAATTCTACCATAAGGTAAACCACCAATTCCTAAAATACAA
>CONTIG976 (46bp)
ATTTAACAGATCCATCTCCAAGTACTACTTTTACTCATTTAGACTC
>CONTIG977 (46bp)
ACATTAATATTAAAAAAAAAAAAATAATATGTTTAAAATTTTACAT
>CONTIG978 (46bp)
TTTATTATTTTTTAAAAAACAAATACTTATAATTTTATTTGTTCTA
>CONTIG979 (46bp)
TTTTTTTTAATTTTAAAATATATTATTTTAAAAAAAAACTTTAATA
>CONTIG980 (46bp)
CATGCTAAATCTCCTTTTTTTGATTTAAAATTAGAAATTGGTTCTT
>CONTIG981 (46bp)
TTTTTTAAATTAAATATGTAAATAAAAGTTTTTTTTTTTTTGCTAA
>CONTIG982 (46bp)
GCAAAACCTTTTTTAATTGATAAATTGATTAGTTTTTTAAAATCTG
>CONTIG983 (46bp)
AATTATTATTGGTGGTGCATTAATGAAAATTTTGTGCAATAATAAA
>CONTIG984 (46bp)
TTCGTTATATATTTAATTATATTATTAATAATATTAAAAAAAAAAC
>CONTIG985 (46bp)
TTTTTTTTTCATTAATAAAAATTACTTTTTAATATAAATTATAAAA
>CONTIG986 (46bp)
CTGCTTTCAATAATAATTTTATTAATTTACAAGCAATTAAAGATGT
>CONTIG987 (46bp)
TTTAATGTATATAGGAAATATAACTTCAAATGTAATTAGATCATTA
>CONTIG988 (46bp)
TTAAATTATTTATTGAATTTAAATAAACGGAATTAATTTTTAACTT
>CONTIG989 (46bp)
AATAGCAATAGTTTAATTATTGTTTATTCTTGTGGAAGTCATAGTA
>CONTIG990 (46bp)
TATTTTATTCTTCATTTTGTAAAATTGTTAATTAAAAAATAAATTG
>CONTIG991 (46bp)
GTTTTATCAATAATTGAATAAACATTAAAATTACTTTTTTTTTTTT
>CONTIG992 (46bp)
ATATTATAATAAAGTAAAACTTTCTCATATTTATATTGATTATGCA
>CONTIG993 (46bp)
GATGTTCATTCTTCTATAATACCATACAAAAATGATACAAAAAAAA
>CONTIG994 (46bp)
ATTATTTTATTTAGAATTACTAATTTTCTAATTGATTCAATAACAT
>CONTIG995 (46bp)
TTGAAAGTTATAATGTTTTAAAATATAAAAATATTATTAATAAAAG
>CONTIG996 (46bp)
AATATTTTATATTGAAATTTTTTTTTTCAGTTTAAGCATACTTCTT
>CONTIG997 (46bp)
ATTAATAGGTTTAGTTGGACTAGTAATTATTTGTATAAAAAACCAA
>CONTIG998 (46bp)
ATATTTAGAATTAAACTTTTTAAATACTGTTATGAAATCAATAAGA
>CONTIG999 (46bp)
TAATTTTTTAATTTTATATTTATAAACATTAATAATTCTTATACAA
>CONTIG1000 (46bp)
GTGCTCCGGATATTGCTAACAAAAATATTGCAAATCCTGTAGGTGC
>CONTIG1001 (46bp)
ATAAAGGATTTGGCGATCCTAAAATTCACAATTTTTTTTTAATTAA
>CONTIG1002 (46bp)
ATTTAGATTTATAAATTTTTGGTATAATAATTGGAAGATGGGATTT
>CONTIG1003 (46bp)
AACTTCCTTCAATTTTTATATATGGTGGTACAATTTTACCAGGCAA
>CONTIG1004 (46bp)
ATAATACCATACAAAAATGATATAAAAAAAAGTATAATAAAACATA
>CONTIG1005 (46bp)
ACATTAGATATTTTTTTATTAGTAATTAAAATTAAACAATTTTCAA
>CONTIG1006 (46bp)
ATCATTAAACCATTTTACTGTTCCGTTTGACATATTTTTTTTTTAA
>CONTIG1007 (46bp)
TTATTTTTACTTCTGGATTAAAATAAGAATGAAATACTCTGCCCGT
>CONTIG1008 (46bp)
TTATTTATTATCTCGTAATTTTACATTTATAAAATAACTTTTTGTA
>CONTIG1009 (46bp)
AGTTTTTTTTATAAATTCTAAATTTAACAATAATCTTTTTACTTCA
>CONTIG1010 (46bp)
AGCAAACATTTGATGTTGTCCAACATCTGCAGAAATAAAGTATTTT
>CONTIG1011 (46bp)
CATTCATTATAACTTTATATATTTATATTACATGATTTATTAGTTA
>CONTIG1012 (46bp)
TTATTAATTTTTTTTAAAAAAAATATTTCTTTTTTAAGAATAAAAT
>CONTIG1013 (46bp)
TTTAAAAAAAGGAGTATGTCTCCCAACTTCTTCTTTTGATAATATA
>CONTIG1014 (46bp)
CCACAAAAAGAACATGTTTCGGAGGACACTCCTTCCATTAAACATC
>CONTIG1015 (46bp)
AAAAAATTGTAAAAAAAAATATAAAAGAAAAAATAGCTACTTTATC
>CONTIG1016 (46bp)
CTATTTTTTCTTTTACTAAAAGTAAATTTTCACTATCTGTCTATAC
>CONTIG1017 (46bp)
TTTAAAAAAAAATTTTTTGAAATTTCTGACAAAGAATTAAAAATCT
>CONTIG1018 (46bp)
GCACCGGTCTCAAAGCTTAAAAAAGCTTTACTATTTCTAGTAAATT
>CONTIG1019 (46bp)
ATTAATAGCAAAAAATTTTAATTGTGTGATAATATCGTTATCGGAA
>CONTIG1020 (46bp)
TAACTAAAAATTAAAAAAATTGAATTAAAAATCTCATACAATTTAT
>CONTIG1021 (46bp)
TTAATAATGGGTTATTGCTTTTATTATGTGATGTTATAAAAAAAAT
>CONTIG1022 (46bp)
AAATTGATAAATTAAAACAATTACTTGAAATATATTTATGTATTTT
>CONTIG1023 (46bp)
AAAACCATCAACTCTAATTTCTGAAAAATGAACATATAAATCATCT
>CONTIG1024 (46bp)
TTAATTTTTTTTATTTTATAATTAATAAAAATAAAATTAAAAATAA
>CONTIG1025 (46bp)
AAAAAAGATATTAATAGTTTTAATAAAACTGATTAAAAACAAACAT
```

**CONTIG1 (159662bp) above is actually a perfect single contig assembly of the Carsonella ruddii genome.**

### Run [QUAST](https://quast.sourceforge.net/) to check our assembly quality against the reference

[QUAST](https://quast.sourceforge.net/) is a quality assessment tool for genome assemblies and what we use to check the quality of our work in this project, it supports aligning to known reference genomes when available.

The QUAST archive (currently a .tar.gz) can be downloaded [here](https://sourceforge.net/projects/quast/files/latest/download).

Expand that and navigate to the folder created before performing the steps below.

```
# we save everything from the line beginning at ">CONTIG1" in the output (e.g. as above) and generate an "output.fasta" for you (see e.g. expected standard output saved in noisyCarsonellaRuddiiReadsAssemblyOutput.txt and compare to the output.fasta generated when you run the command), verify output.fasta was created or copy everything after ">CONTIG1" in the standard output and create it by hand (you can also change the output filename with the `--output_fasta` argument)
# then you can run the following from a freshly unpacked quast folder

$ ./quast.py ~/DeBruijnGenomeAssembler/output.fasta -r ~/DeBruijnGenomeAssembler/util/carsonella_ruddii_complete_genome_AP009180.1.fasta -o ~/carsonella_ruddii_quast_output
/home/liveuser/quast-5.3.0/quast_libs/qutils.py:224: SyntaxWarning: invalid escape sequence '\w'
  value = convert_to_unicode(re.sub('[^\w\s-]', '-', value).strip())
/home/liveuser/quast-5.3.0/quast_libs/qutils.py:225: SyntaxWarning: invalid escape sequence '\s'
  value = convert_to_unicode(re.sub('[-\s]+', '-', value))
/home/liveuser/quast-5.3.0/quast_libs/qutils.py:715: SyntaxWarning: invalid escape sequence '\d'
  version_pattern = re.compile('(?P<major_version>\d+)\.(?P<minor_version>\d+)')
/home/liveuser/quast-5.3.0/quast_libs/qutils.py:755: SyntaxWarning: invalid escape sequence '\d'
  version_pattern = re.compile('version "(\d+\.\d+)')
/home/liveuser/quast-5.3.0/quast_libs/reporting.py:617: SyntaxWarning: invalid escape sequence '\g'
  tex_file.write('\\caption{All statistics are based on contigs of size $\geq$ %d bp, unless otherwise noted ' % qconfig.min_contig + \
/home/liveuser/quast-5.3.0/quast_libs/reporting.py:618: SyntaxWarning: invalid escape sequence '\#'
  '(e.g., "\# contigs ($\geq$ 0 bp)" and "Total length ($\geq$ 0 bp)" include all contigs).}\n')
/home/liveuser/quast-5.3.0/quast_libs/reporting.py:663: SyntaxWarning: invalid escape sequence '\h'
  row = "\hspace{5mm}" + row.lstrip()
/home/liveuser/quast-5.3.0/quast_libs/reporting.py:665: SyntaxWarning: invalid escape sequence '\h'
  row = "\hspace{2mm}" + row.lstrip()
/home/liveuser/quast-5.3.0/./quast.py /home/liveuser/DeBruijnGenomeAssembler/output.fasta -r /home/liveuser/DeBruijnGenomeAssembler/util/carsonella_ruddii_complete_genome_AP009180.1.fasta -o /home/liveuser/carsonella_ruddii_quast_output

Version: 5.3.0, c3eb988a

System information:
  OS: Linux-6.14.0-63.fc42.x86_64-x86_64-with-glibc2.41 (linux_64)
  Python version: 3.13.2
  CPUs number: 12

Started: 2025-09-17 23:21:21

Logging to /home/liveuser/carsonella_ruddii_quast_output/quast.log
NOTICE: Maximum number of threads is set to 3 (use --threads option to set it manually)

CWD: /home/liveuser/quast-5.3.0
Main parameters: 
  MODE: default, threads: 3, min contig length: 500, min alignment length: 65, min alignment IDY: 95.0, \
  ambiguity: one, min local misassembly length: 200, min extensive misassembly length: 1000

Reference:
  /home/liveuser/DeBruijnGenomeAssembler/util/carsonella_ruddii_complete_genome_AP009180.1.fasta ==> carsonella_ruddii_complete_genome_AP009180.1

Contigs:
  Pre-processing...
/home/liveuser/quast-5.3.0/quast_libs/site_packages/joblib3/func_inspect.py:51: SyntaxWarning: invalid escape sequence '\<'
  '\<doctest (.*\.rst)\[(.*)\]\>',
/home/liveuser/quast-5.3.0/quast_libs/site_packages/joblib3/_memory_helpers.py:10: SyntaxWarning: invalid escape sequence '\s'
  cookie_re = re.compile("coding[:=]\s*([-\w.]+)")
  /home/liveuser/DeBruijnGenomeAssembler/output.fasta ==> output
/home/liveuser/quast-5.3.0/quast_libs/html_saver/html_saver.py:324: SyntaxWarning: invalid escape sequence '\/'
  ref_name = re.findall('<a.*>(.*)<\/a>', line)[0]
/home/liveuser/quast-5.3.0/quast_libs/ra_utils/misc.py:212: SyntaxWarning: invalid escape sequence '\S'
  chr_name_pattern = 'SN:(\S+)'
/home/liveuser/quast-5.3.0/quast_libs/ra_utils/misc.py:213: SyntaxWarning: invalid escape sequence '\d'
  chr_len_pattern = 'LN:(\d+)'

2025-09-17 23:21:34
Running Basic statistics processor...
  Reference genome:
    carsonella_ruddii_complete_genome_AP009180.1.fasta, length = 159662, num fragments = 1, GC % = 16.56
  Contig files: 
    output
  Calculating N50 and L50...
    output, N50 = 159662, L50 = 1, auN = 159662.0, Total length = 159662, GC % = 16.56, # N's per 100 kbp =  0.00
  Drawing Nx plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/basic_stats/Nx_plot.pdf
  Drawing NGx plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/basic_stats/NGx_plot.pdf
  Drawing cumulative plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/basic_stats/cumulative_plot.pdf
  Drawing GC content plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/basic_stats/GC_content_plot.pdf
  Drawing output GC content plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/basic_stats/output_GC_content_plot.pdf
Done.

2025-09-17 23:21:36
Running Contig analyzer...
Compiling Minimap2 (details are in /home/liveuser/quast-5.3.0/quast_libs/minimap2/make.log and make.err)
  output
  Logging to files /home/liveuser/carsonella_ruddii_quast_output/contigs_reports/contigs_report_output.stdout and contigs_report_output.stderr...
  Aligning contigs to the reference
  Analysis is finished.
  Creating total report...
    saved to /home/liveuser/carsonella_ruddii_quast_output/contigs_reports/misassemblies_report.txt, misassemblies_report.tsv, and misassemblies_report.tex
  Transposed version of total report...
    saved to /home/liveuser/carsonella_ruddii_quast_output/contigs_reports/transposed_report_misassemblies.txt, transposed_report_misassemblies.tsv, and transposed_report_misassemblies.tex
  Creating total report...
    saved to /home/liveuser/carsonella_ruddii_quast_output/contigs_reports/unaligned_report.txt, unaligned_report.tsv, and unaligned_report.tex
  Drawing misassemblies by types plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/contigs_reports/misassemblies_plot.pdf
  Drawing misassemblies FRCurve plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/contigs_reports/misassemblies_frcurve_plot.pdf
Done.

2025-09-17 23:22:01
Running NA-NGA calculation...
  output, Largest alignment = 159662, NA50 = 159662, NGA50 = 159662, LA50 = 1, LGA50 = 1
  Drawing cumulative plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/aligned_stats/cumulative_plot.pdf
  Drawing NAx plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/aligned_stats/NAx_plot.pdf
  Drawing NGAx plot...
    saved to /home/liveuser/carsonella_ruddii_quast_output/aligned_stats/NGAx_plot.pdf
Done.
/home/liveuser/quast-5.3.0/quast_libs/search_references_meta.py:413: SyntaxWarning: invalid escape sequence '\['
  seqname = re.sub('[\[,/\]]', ';', seqname)

2025-09-17 23:22:02
Running Genome analyzer...
  NOTICE: No file with genomic features were provided. Use the --features option if you want to specify it.

  NOTICE: No file with operons were provided. Use the -O option if you want to specify it.
  output
  Analysis is finished.
  Skipping drawing Genome fraction, % histogram... (less than 2 columns histogram makes no sense)
Done.

NOTICE: Genes are not predicted by default. Use --gene-finding or --glimmer option to enable it.

2025-09-17 23:22:02
Creating large visual summaries...
This may take a while: press Ctrl-C to skip this step..
  1 of 2: Creating PDF with all tables and plots...
  2 of 2: Creating Icarus viewers...
Done

2025-09-17 23:22:02
RESULTS:
  Text versions of total report are saved to /home/liveuser/carsonella_ruddii_quast_output/report.txt, report.tsv, and report.tex
  Text versions of transposed total report are saved to /home/liveuser/carsonella_ruddii_quast_output/transposed_report.txt, transposed_report.tsv, and transposed_report.tex
  HTML version (interactive tables and plots) is saved to /home/liveuser/carsonella_ruddii_quast_output/report.html
  PDF version (tables and plots) is saved to /home/liveuser/carsonella_ruddii_quast_output/report.pdf
  Icarus (contig browser) is saved to /home/liveuser/carsonella_ruddii_quast_output/icarus.html
  Log is saved to /home/liveuser/carsonella_ruddii_quast_output/quast.log

Finished: 2025-09-17 23:22:02
Elapsed time: 0:00:41.851998
NOTICEs: 4; WARNINGs: 0; non-fatal ERRORs: 0

Thank you for using QUAST!
```

![](noisyCarsonellaRuddiiReadsAssemblyQuast.png)

This result is consistent for different generations of noised Carsonella ruddii data (try yourself different random seeds when generating the random reads) and is not super sensitive to the choice of k, e.g. `--k 40` works as well as `--k 45`.

### More examples with larger genomes

See [the dataset generation tool README for instructions on generating synthetic sequencer data for larger genomes and how to trigger assembly for them](../util/README.md#nz_cabfnn0100000011-escherichia-coli-j53-isolate-whole-genome-sequence).

### License

Copyright (C) 2025 William Bruns

This documentation comes with ABSOLUTELY NO WARRANTY and is provided "AS IS".
This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See [the full license](../LICENSE) .