#!/usr/bin/bash

#SBATCH --job-name=split_fastq       ## Name of the job for the scheduler
#SBATCH --account=siwase99         ## Your PI's uniqname plus 0,99, or other number
#SBATCH --partition=largemem         ## Choose: standard, largemem, gpu, spgpu, debug            
##SBATCH --gpus=1  			   ## if partition=gpu, number of GPUS needed
#SBATCH --nodes=1                      ## number of nodes you are requesting
#SBATCH --ntasks=1                     ## how many task spaces do you want to reserve
#SBATCH --cpus-per-task=16            ## how many cores do you want to use per task
#SBATCH --time=30:00               ## Maximum length of time you are reserving the 
#SBATCH --mem=10MB                      ## Memory limit  
#SBATCH --mail-user=mgavilan@umich.edu  ## send email notifications to umich email listed
#SBATCH --mail-type=ALL                ## when to send email (standard values are:
                                       ## NONE,BEGIN,END,FAIL,REQUEUE,ALL.

module load Bioinformatics
module load seqtk/1.3

for file in *.fastq; do
    base=$(basename "$file" .fastq)  # Extract base name without extension
    echo "Processing $file..."
    seqtk seq -1 "$file" > "${base}_R1.fastq"
    seqtk seq -2 "$file" > "${base}_R2.fastq"
done