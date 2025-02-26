####### Total reads
module load samtools
samtools view -c -h sort_neurons_mut1.bam chr11:45938065-45949400 > neurons_mut1_total.txt

### Exon 13
samtools view -c -h sort_293T_mut1.bam chr11:45948886-45948946 > 293T_mut1_Ex13.txt
### Exon 14 - C
samtools view -c -h sort_neurons_mut1.bam chr11:45945840-45946030 > neurons_mut1_Ex14_c.txt
### Exon 14 - N
samtools view -c -h sort_293T_mut1.bam chr11:45946076-45946098 > 293T_mut1_Ex14_n.txt
### Exon 15 
samtools view -c -h sort_293T_mut1.bam chr11:45938157-45938312 > 293T_mut1_Ex15.txt


######## For introns counts
conda activate bioinfo_env
#####3# To make intron file:
awk 'BEGIN { OFS="\t" }
{
    # Ensure proper calculations between consecutive exons
    if (NR > 1 && prev_chr == $1 && prev_strand == $4) {
        intron_start = prev_end
        intron_end = $2 - 1

        # Only print valid introns
        if (intron_start < intron_end) {
            print prev_chr, intron_start, intron_end, prev_strand
        }
    }
    # Store current exon information
    prev_chr = $1; prev_end = $3; prev_strand = $4;
}' phf21a_exons_adjusted.bed > phf21a_introns.bed


INTRON_BED="/nfs/turbo/umms-siwase/CECILIA/MiSeq_patient_mutation/mini_gene_fastq/split_reads/bam/phf21a_introns.bed" 

bedtools intersect -a phf21a_introns.bed -b neurons_wt1_hisat2_sorted.bam -c > neurons_wt1_introns.txt
