#!/bin/bash

mkdir -p filtered_tsvs

for vcf in *.vcf.gz; do
    sample_name=$(basename "$vcf" .vcf.gz)
    echo "Processing $sample_name..."

    bcftools view -i 'QUAL>=30 && INFO/DP>=10' "$vcf" |
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AC\t%INFO/AN\n' |
    awk -v sample="$sample_name" 'BEGIN {
        OFS="\t";
        print "Sample", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "DP", "AC", "AN", "AF";
    }
    {
        ac=$8; an=$9;
        if (an > 0) {
            af = ac / an;
        } else {
            af = 0;
        }
        if (af > 0.05)
            print sample, $1, $2, $3, $4, $5, $6, $7, ac, an, af;
    }' > "filtered_tsvs/${sample_name}_filtered.tsv"
done




