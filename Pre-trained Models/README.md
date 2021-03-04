# Pre-trained models
The pre-trained models built in the leave-one-individual-out cross-validation step were released in this folder.<br>

The model 'gatk_indel_HG001.model' indicates that we trained the model using all gold-standard Indels from the GATK analysis results of the five individuals (HG002, HG003, HG004, HG006, HG007).<br>
All gold-standard Indels from the HG001 sample was used as testing data.<br>

Algorithm
------------
These models were trained by using XGBoost with the parameters in the 'Reproduce Results/Training and Evaluation/ML_parameters.json'
XGBoost demonstrated a significant improvement than other machine learning (ML) methods in both SNV classification and Indel classification.
![](https://github.com/yyren/FNVC/raw/master/Picture/Comparison_of_different_ML_methods.png)<br>


Evaluation
------------
The filtering performacne may be overestimated by using a subset of WGS variants. Therefore, these released models were tested on all gold-standard variants of one individual WGS resutls.<br>
![](https://github.com/yyren/FNVC/raw/master/Picture/Performance_on_tesing_data_with_different_imbalanced_ratios.png)<br>

Retraining
------------
FNVC support retraining or named as incremental training a new model based on the pre-trained model and the additional data

perl FNVC_Main.pl \
	-t retrain -c gatk \
	-tp TP.vcf -fp FP.vcf \
	-snpm gatk_snp.bin -im gatk_indel.bin \
	-o /home/ubuntu/gatk_retrained