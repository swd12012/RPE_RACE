### RPE RACE

I took RNA from untreated 211fl/fl RPE (P1) and 211fl/fl RPE (P2) that was injected with AAV-cre.

I performed nested RACE on these samples and sent the PCR products for sequencing.

It was run on a MiSeq.

Fastq files were aligned to the GRCm38 mouse genome with hisat2 (v2.2.1) and counted with featureCounts (subread v2.0.1).

Downstream analysis was performed with R and the DESeq2 package.
