# FNVC incremental learning
Using the standardized vcf file

singularity exec variant_calling_image.sif \
    perl FNVC_Main.pl -t filter -c gatk \
    -snpm gatk_snp_example.model -im gatk_indel_example.model \
    -i example_standardized.vcf \
    -o example_output.vcf

