Introduction
------------
This folder contains the scripts used for measuring the performance of different software.

Users can double check the results by using these scripts and the corresponding model (in the '../models' folder) and the [data] (https://drive.google.com/open?id=17AOz1SChj22PTlHGPTLNlBrgjk5ExYHV)

Evaluation
------------
For example if we evaluate the FNVC-gatk model by using 'HG001' SNV data, we can used the:

Data:

'HG001_TP_snp_features.record' in the folder 'example_FNVC_gatk'
'HG001_FP_snp_features.record' in the folder 'example_FNVC_gatk'

Model:
['FNVC_gatk_snp_model.bin'](https://github.com/yyren/FNVC/blob/master/models/FNVC-GATK/HG001/FNVC_gatk_snp_model.bin) in 'models/FNVC-GATK/HG001/FNVC_gatk_snp_model.bin'

Result:
'FNVC_gatk_snp_statistic_summary.txt'

FNVC-gatk model:<br>
		python FNVC_gatk_statistic.py
			--test_TP_file HG001_TP_snp_features.record --test_FP_file HG001_FP_snp_features.record \
			--in_model FNVC_gatk_snp_model.bin --out_file FNVC_gatk_snp_statistic_summary.txt \
			--acc_cutoff 0.5126 --mcc_cutoff 0.5 
			
All the data used for training and test can be download from [here](https://drive.google.com/open?id=17AOz1SChj22PTlHGPTLNlBrgjk5ExYHV)

