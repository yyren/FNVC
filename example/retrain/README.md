# FNVC retraining

singularity exec variant_calling_image.sif \
	perl FNVC_Main.pl -t retrain -c gatk \
	-tp gatk_example_tp.vcf -fp gatk_example_fp.vcf \
	-o /example/incremental/gatk_retrain

-o: prefix name of the model <br>
