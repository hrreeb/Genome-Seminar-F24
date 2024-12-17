mkdir -p /scratch/biol726311/gouldian/hisat_mapping/namesorted_bam_hisat_align

samtools sort -n -o $2 $1
