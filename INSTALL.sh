#!/bin/bash
# builds the dependent modules for FVC
# developed by ryy, yongyong.ren@sjtu.edu.cn

modules=('numpy' 'pandas' 'xgboost' 'os' 'argparse' 'math' 're' 'collections' 'itertools' 'datetime' 'sys' 'pysam')

for module in ${modules[*]}
do
	python -c "import $module"
	if [[ $? -eq 0 ]]; then
		echo "checking $module: exist"
	else
		echo "Checking : $module not exist"
		echo "Starting Install: $module"
		pip install --user $module
		python -c "import $module"
		if [[ $? -ne 0 ]]; then
			echo "$module: Failed Install !!!"
			break
		fi
	fi

done
