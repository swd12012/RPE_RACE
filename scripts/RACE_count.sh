#!/bin/bash
#SBATCH --job-name=RACE_count      ## Name of the job.
#SBATCH -A ecoevo283         		 ## account to charge 
#SBATCH -p standard           		 ## partition/queue name
#SBATCH --cpus-per-task=4  		     ## number of cores the job needs, can the programs you run make used of multiple cores?

module load subread/2.0.1
gtf='/data/homezvol2/swdu/ref/mmu/Mus_musculus.GRCm38.102.gtf.gz'
myfile=`cat /dfs6/pub/swdu/hts.igb.uci.edu/swdu21030345/sample_coding.txt | cut -f3`
featureCounts -p -T 4 -t exon -g gene_id -Q 30 -F GTF -a $gtf -o /dfs6/pub/swdu/hts.igb.uci.edu/swdu21030345/out/RACE_counts.txt $myfile