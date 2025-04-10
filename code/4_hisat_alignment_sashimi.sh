#!/usr/bin/bash
#SBATCH --job-name=hisat2         ## Name of the job for the scheduler
#SBATCH --account=siwase99         ## Your PI's uniqname plus 0,99, or other number
#SBATCH --partition=largemem         ## Choose: standard, largemem, gpu, spgpu, debug            
##SBATCH --gpus=1  			   ## if partition=gpu, number of GPUS needed
#SBATCH --nodes=1                      ## number of nodes you are requesting
#SBATCH --ntasks=1                     ## how many task spaces do you want to reserve
#SBATCH --cpus-per-task=16            ## how many cores do you want to use per task
#SBATCH --time=5:00:00                 ## Maximum length of time you are reserving the 
#SBATCH --mem=120G                      ## Memory limit  
#SBATCH --mail-user=mgavilan@umich.edu  ## send email notifications to umich email listed
#SBATCH --mail-type=ALL                ## when to send email (standard values are:
                                       ## NONE,BEGIN,END,FAIL,REQUEUE,ALL.
   ## (See documentation for others)
#SBATCH --output=./%x-%j               ## send output and error info to the file listed
                                       ##(optional: different name format than default) 

module load Bioinformatics
module load hisat2/2.2.0-r3kcuwd
module load samtools

# MUT_Exon6
hisat2 -p 16 --dta -x /nfs/turbo/umms-siwase/GENOME_INDICES/HG38/ht2_hg38-index -1 DM010_Mut_Exon6_331786w_JD8_R1.fastq -2 DM010_Mut_Exon6_331786w_JD8_R2.fastq --known-splicesite-infile hg38_splice_sites.txt --novel-splicesite-outfile sample_novel_splice_sites.txt --pen-noncansplice 100 --max-intronlen 500000 -S MUT_Exon6.sam

samtools view -bS MUT_Exon6.sam > MUT_Exon6.bam
samtools sort MUT_Exon6.bam > MUT_Exon6.sorted.bam
samtools index MUT_Exon6.sorted.bam
