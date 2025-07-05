#!/bin/bash

while getopts "n:i:o:r:" opt; do
	case $opt in
		n) sample_name="$OPTARG" ;;
		i) input_directory="$OPTARG" ;;
		o) output_directory="$OPTARG" ;;
		r) reference="$OPTARG" ;;
		?) echo "This tool is used for RNA-seq analysis. It requires the following parameters:" 
			echo "Usage: $0 -i <input_dir> -o <output_dir> -r <reference> -n <sample_name>"
			exit 1 
			;;
		
		
	esac
done
start_date=$(date) 
echo "Start date: $start_date" 

echo "Sample name: $sample_name"
echo "Input directory: $input_directory"
echo "Output directory: $output_directory"
echo "Reference: $reference"
echo "----------------------"
echo "\n"

# ---- Implementing a function for exception handling 






# Stage 1: Preprocessing

## Trimming (skewer) 
echo "Trimming reads with skewer..."
# ---- Including a file format checker
if find ${input_directory} -name "${sample_name}_R1.fastq.gz" | grep -q .; then
	echo "Found file format in _R1, _R2 format...." 
	skewer -r 0.1 -d 0.03 --min 18 ${input_directory}/${sample_name}_R1.fastq.gz ${input_directory}/${sample_name}_R2.fastq.gz -o ${output_directory}/${sample_name}	

else
	echo "Found alternate format..." 
	skewer -r 0.1 -d 0.03 --min 18 ${input_directory}/${sample_name}_1.fastq.gz ${input_directory}/${sample_name}_2.fastq.gz -o ${output_directory}/${sample_name} 
fi
# I also need to trap the error

skewer_exit_code=$?
if [ $skewer_exit_code -ne 0 ]; then
	echo "There was an error with the skewer. Exiting..." 
	exit 1
fi
 


# Stage 2: Read Quantification (kallisto) 

## Pseudoalignment

## Building the indexes (use the reference index) 
echo "Building Kallisto Index..."
kallisto index -i ${output_directory}/GRCh38.idx $reference
kallisto_index_exit_code=$?
if [ $kallisto_index_exit_code -ne 0 ]; then 
	echo "Error with the Kallisto index generation. Exiting..."
	exit 1
fi


## Quantiifcation algorithm (reference and then index)
echo "Quantifying reads with Kallisto..."

# Format: SRR19383351-trimmed-pair1.fastq
# --- making a directory for the sample name and storing the quant information in there

if find "${output_directory}" -type d -name "results_${sample_name}" | grep -q .; then
	echo "Found exiting directory file for the sample" 
else
	echo "Creating directory to store the abundance results" 
	mkdir ${output_directory}/results_${sample_name}/

fi
kallisto quant -i ${output_directory}/GRCh38.idx ${output_directory}/${sample_name}-trimmed-pair1.fastq ${output_directory}/${sample_name}-trimmed-pair2.fastq -o ${output_directory}/results_${sample_name} 

kallisto_quant_exit_code=$?
if [ $kallisto_quant_exit_code -ne 0 ]; then
	echo "Error with the Kallisto Quantification. Exiting..." 
	exit 1
fi

# moving all the sample-associated results to the same file 
mv ${output_directory}/${sample_name}* ${output_directory}/results_${sample_name}

echo "The program is now complete!"

