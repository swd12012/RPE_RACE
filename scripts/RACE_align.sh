#!/bin/bash
#SBATCH --job-name=RACE_align      ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --array=1-2         ## number of tasks to launch, given hint below wc -l $file is helpful
#SBATCH --cpus-per-task=2    ## number of cores the job needs, can the programs you run make used of multiple cores?

module load samtools/1.10
module load hisat2/2.2.1

ref='/data/homezvol2/swdu/ref/mmu/Mus_musculus.GRCm38.dna.toplevel.fa'
file='/dfs6/pub/swdu/hts.igb.uci.edu/swdu21030345/sample_coding.txt'
prefix=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1 | cut -f2`
rawdir='/dfs6/pub/swdu/hts.igb.uci.edu/swdu21030345/raw/'
outdir='/dfs6/pub/swdu/hts.igb.uci.edu/swdu21030345/bam/'

hisat2 -p 2 -x $ref -1 $rawdir${prefix}-READ1.fq.gz -2 $rawdir${prefix}-READ2.fq.gz -S $outdir${prefix}.sam
samtools view -bS $outdir${prefix}.sam > $outdir${prefix}.bam
samtools sort $outdir${prefix}.bam -o $outdir${prefix}.sorted.bam
samtools index $outdir${prefix}.sorted.bam