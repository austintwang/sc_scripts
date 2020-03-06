# module load seq/tabix
# module load seq/bcftools
# export BCFTOOLS_PLUGINS="/home/ag284/apps/bcftools/bcftools-1.3.1/plugins"

PPATH="/agusevlab/awang/plink/plink"

for BED in /agusevlab/awang/sc_kellis/gen/impute/*.bed; do
    PRE="${BED%.*}"
    OUT=$PRE

    # remove strand-ambiguous SNPs, non-autosomes, and indels
    cat $PRE.bim | awk '$1 != "X" && $6 != "I" && $6 != "D" && $5 != 0 && $6 != 0 && !( ($5=="A"&&$6=="T") || ($5=="T"&&$6=="A") || ($5=="G"&&$6=="C") || ($5=="C"&&$6=="G") ) { print $2 }' | sort | uniq -u > $OUT.extract

    # create a VCF, zip, and index
    $PPATH --bfile $PRE --extract $OUT.extract --recode vcf --out $OUT
    bgzip -c $OUT.vcf > $OUT.vcf.gz
    tabix $OUT.vcf.gz

    # use the fixref plugin to identify and flip any nonta-matching SNPs
    bcftools +fixref $OUT.vcf.gz -Oz -o $OUT.fixed.vcf.gz -- -f /agusevlab/awang/refs/human_g1k_v37.fasta -m flip
    tabix $OUT.fixed.vcf.gz

    # finally, double check that there are no remaining errors
    bcftools norm --check-ref w -f /agusevlab/awang/refs/human_g1k_v37.fasta $OUT.fixed.vcf.gz -o /dev/null > $OUT.ERRLOG
done