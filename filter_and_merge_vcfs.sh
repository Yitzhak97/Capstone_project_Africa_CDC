#!/bin/bash

OUTPUT_FILE="hq_allvariants.tsv"
echo -e "Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tDP\tAC\tAN\tAF" > "$OUTPUT_FILE"

for vcf in *.vcf.gz; do
    echo "Processing $vcf..."
# Extract sample name by removing .vcf.gz
    sample_name=$(basename "$vcf" .vcf.gz)

# Apply filters and extract fields
bcftools view -i 'QUAL>=30 && INFO/DP>=10' "$vcf" |
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AC\t%INFO/AN\n' |
    awk -v sample="$sample_name" 'BEGIN {OFS="\t"} {
        ac=$8; an=$9;
        if (an > 0) {
            af = ac / an;
        } else {
            af = 0;
        }
        if (af > 0.05) {
            print sample, $1, $2, $3, $4, $5, $6, $7, ac, an, af;
        }
    }' >> "$OUTPUT_FILE"

done

