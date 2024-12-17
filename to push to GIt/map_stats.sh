mkdir -p /scratch/biol726311/gouldian/hisat_mapping/mapping_stats

samtools flagstat $1 > $2

