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
date 
echo "Sample name: $sample_name"
echo "Input directory: $input_directory"
echo "Output directory: $output_directory"
echo "Reference: $reference"
echo "----------------------"
echo "\n"

# Stage 1: Preprocessing

## Trimming (skewer) 
echo "Trimming reads with skewer..."
skewer -r 0.1 -d 0.03 --min 18 ${input_directory}/${sample_name} -o ${output_directory}/${sample_name}	

# Stage 2: Read Quantification (kallisto) 

## Pseudoalignment

## Building the indexes (use the reference index) 
echo "Building Kallisto Index..."
kallisto index -i ${output_directory}/GRCh38.idx $reference

## Quantiifcation algorithm (reference and then index)
echo "Quantifying reads with Kallisto..."
kallisto quant -i ${output_directory}/GRCh38.idx 	${output_directory}/${sample_name}_1-trimmed.fastq ${input_directory}/${sample_name}_2-trimmed.fastq

# the next stages would then be diffbind

# reference = /gpfs/home/wmp21vtu/scratch/reference_genome_hg38_86/gencode.v25.transcripts.fa


