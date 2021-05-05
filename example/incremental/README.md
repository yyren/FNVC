# FNVC incremental learning

singularity exec variant_calling_image.sif \
    perl FNVC_Main.pl -t incremental -c gatk \
    -snpm gatk_snp_example.model -im gatk_indel_example.model \
    -tp gatk_example_tp.vcf -fp gatk_example_fp.vcf

output: gatk_snp_example.model.incremental gatk_indel_example.model.incremental
