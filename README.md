### RPE RACE

I took retinal pigment epithelium (RPE) RNA from an untreated 211fl/fl mouse (P1) and a 211fl/fl mouse (P2) that was injected with AAV-cre.

I performed nested RACE on these samples and sent the PCR products for sequencing.

It was run on a MiSeq.

Fastq files were aligned to the GRCm38 mouse genome with hisat2 (v2.2.1) and counted with featureCounts (subread v2.0.1).

Downstream analysis was performed with R and the DESeq2 package.

Because there was only one sample in each run as a pilot, DESeq2 didn't really know what to do with the data, so I duplicated each sample and had it run a Wald test.

Unfortunately, when I plotted the data, there were 3000+ DEGs, and the volcano plot also testifies to the fact that my study was underpowered.

I did make a heatmap, though the ENSEMBL IDs were pretty useless as well.

I put the ENSEMBL IDs into [this website](https://biotools.fr/mouse/ensembl_symbol_converter) to generate a gene symbol list:

| ENSEMBL ID         | Gene symbol   |
|--------------------|---------------|
| ENSMUSG00000086503 | Xist          |
| ENSMUSG00000102289 | Gm31258       |
| ENSMUSG00000095186 | Gm10718       |
| ENSMUSG00000070385 | Ampd1         |
| ENSMUSG00000104938 | Gm42840       |
| ENSMUSG00000024843 | Chka          |
| ENSMUSG00000110804 | Olfr1084      |
| ENSMUSG00000112800 | Gm40604       |
| ENSMUSG00000047747 | Rnf150        |
| ENSMUSG00000031358 | Msl3          |
| ENSMUSG00000024654 | Asrgl1        |
| ENSMUSG00000021301 | Hecw1         |
| ENSMUSG00000057110 | Cntrl         |
| ENSMUSG00000041594 | Tmtc4         |
| ENSMUSG00000079262 | Slco1a6       |
| ENSMUSG00000014767 | Tbp           |
| ENSMUSG00000019872 | Smpdl3a       |
| ENSMUSG00000034837 | Gnat1         |
| ENSMUSG00000098557 | Kctd12        |
| ENSMUSG00000024587 | Nars          |
| ENSMUSG00000020100 | Slc29a3       |
| ENSMUSG00000097061 | 9330151L19Rik |
| ENSMUSG00000090362 | Vmn2r79       |
| ENSMUSG00000070644 | Etnk2         |
| ENSMUSG00000105345 | BC030343      |
| ENSMUSG00000028174 | Rpe65         |
| ENSMUSG00000117599 |               |
| ENSMUSG00000026889 | Rbm18         |
| ENSMUSG00000090942 | F830016B08Rik |
| ENSMUSG00000027088 | Phospho2      |

Of note, Gnat1 is highly expressed in the neural retina, and RPE65 is highly expressed in the RPE. Gnat1 may have been a hit because of a sample contamination.

Alla in all, this was probably not the correct way to investigate this question, but I will continue investigating the dataset.