# BAMquantify
A bash script that quantifies .bam files using .bed intervals and outputs results into a .csv table. 

## Usage
*Requires samtools version 1.10 or above.*

Execute the script with the following commands (in order):

1) Path to the directory containing .bam files, requires indexing;
2) Path to .bed file containing coordinates of regions of interest, requires 4 columns (chr start end name);
3) Set the run mode to single-end (SE, for example, ChIP-seq) or paired-end (PE, for example, ATAC-seq), default is PE;

For example: ./BAMquantify.sh /usr/home/path_to_bam_files /usr/home/path_to_bed_file.bed PE

Results will be outputed to final_quantification.csv within the given folder containing bam files. 

## Results
The resulting table (final_quantification.csv) will contain the following columns:

1) file_name: contains the name of each analyzed .bam file;
2) read_count: total number of counts in each given .bam file (library size);

The remainder columns will contain the number of counts within each given region (from .bed file) and will be named accordingly.

If the objective is to calculate reads per million from each .bed region, divide counts from each region to read_count from each file.

## Credits
Written by Luis Abatti, based on Linh Huynh's original script
