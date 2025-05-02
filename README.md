# RNA-Seq-Pipeline
I have created this tool as part of my MScR Bioinformatics Project using shell scripting alongside Bioinformatic packages. This is a fast and portable tool that can be executed anywhere with a UNIX/UNIX-like terminal. This tool accepts paried-end fastq files and outputs gene expression profiles in .tsv format. 


<b>NOTE THAT</b>: This work is still in early stages and has not been benchmarked against published datasets. 


## Introduction
RNA-Seq is a powerful tool to quantify the gene expression of samples, employing Next Generation Sequencing (NGS) technologies to create detailed maps of the transcriptome. This can unmask the role of various genes towards phenotypic traits such as hair loss, disease progression, and more. Although powerful, RNA-Seq is commonly used alongside other tools such as ATAC-Seq for a more complete developmental analysis in organisms. This accounts for additional factors such as the epigenome and chromatin accessibility. 

This tool uses Bioinformatics tools within a shell script for RNA-Seq analysis, accepting parameters for the sample name, input path, output path and reference genome path. 


## Prequisites and Caveats
1. This tool is designed for paired-end .fastq files only
2. The required tools and packages are located within the ```<b>requirements.txt</b>``` file 
3. the main analysis code is found within the ```<b>rna_pipeline.sh</b>``` script


## How to use

1. Clone this repository

```git clone https://github.com/Thomas-Nguyen12/RNA-Seq-Pipeline ```

2. Install the required packages
```conda install -r requirements.txt```

3. Run the script
```./rna_pipeline.sh -i <input_dir> -o <output_dir> -r <path_to_reference> -n <sample_name>```

<b>NOTE THAT</b>: if you are running this script in a SLURM environment, you can create a submission script that executes the ```rna_pipeline.sh``` script with a ```.sub``` that contains the command in the above example. For example: 
```sbatch <b>run_rna_pipeline.sub</b> ```

where <b>run_rna_pipeline.sub</b> contains: 
```./rna_pipeline.sh -i <input_dir> -o <output_dir> -r <reference> -n <sample_name>```
with your specified SLURM parameters at the beginning 


## Future Directions
Although promising, this work needs further refinement. First of all, I need to benchmark this work against published data and ensure that the reference data and annotation can be used without prior installation by the user. 

## Acknowledgements
I would like to thank my Primary Supervisor David Monk and my Bioinformatics Mentor for their support in this work. I would not have been able to do this without them. 


## References

