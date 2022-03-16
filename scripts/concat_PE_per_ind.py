#!/bin/env python3

import os, sys, argparse, subprocess
from collections import defaultdict

## This script will concatenates sequences files of the same individuals split on different lanes
## Format ex.: 
#/home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_12_S3_L001_R1_trim_001.fastq.gz
#/home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_12_S3_L002_R1_trim_001.fastq.gz
#/home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_12_S3_L001_R2_trim_001.fastq.gz
#/home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_12_S3_L002_R2_trim_001.fastq.gz
## files contain: "trim", should not contain "single"
## new name: _L001L002


def main():

	## arguments with help messages
	parser = argparse.ArgumentParser(description="Concatenate multiple runs into one file. Assuming paired reads with the following structure: gene_ind_sample_lane_readdirection_trim_number.fastq.gz", allow_abbrev=False)
	parser.add_argument("--files", nargs = "*", help = "File(s) to concatenate. Default: No")
	parser.add_argument("--output_directory", type = str, default = "", help = "Ouptut directory. Default: No")
	parser.add_argument("--single_read_info", type = str, help = "String to exclude unpaired reads. Ex. single", default = "single")

	## Parsing the arguments
	args = parser.parse_args()
	#print("Arguments: ",args,"\n")
	dictionary_lanes_filenames = defaultdict(list)

	# Organize files into groups by their individual and 
	for file in args.files:
		# Only consider paired reads
		if args.single_read_info not in file:
			file_base=os.path.basename(file)

			# Get data from file name
			gene, individual, sample_numb, lane, read, trim, junk = file_base.split("_")

			# Add this data to a dictionary where:
			# key = (gene, individual, sample_numb, trim, junk)
			# value = [(lane_0, read_0, file_0), (lane_1, read_1, file_1), ... (lane_N, read_N, file_N)]
			key = (gene, individual, sample_numb, read, trim, junk)
			dictionary_lanes_filenames[key].append((lane, file))

	# Iterate through the groups of files and cat
	for key, values in dictionary_lanes_filenames.items():
		print("="*20)
		print("cat file:")
		# Get file data from key
		gene, individual, sample_numb, read, trim, junk = key
		print(*key, sep='\t')

		# Sort values based on ?
		values = sorted(values, key=lambda x:x[0])
		print(*values, sep='\n')

		# Get filenames to cat
		lanes = [v[0] for v in values]
		filenames = [v[-1] for v in values]

		# Create new filename
		lane_name = ''.join(lanes)
		new_filename = '_'.join([gene, individual, sample_numb, lane_name, read, trim, junk])
		new_filename = f"{args.output_directory}/{new_filename}"
		print("New filename: ", new_filename)

		# Run cat command
		cmd = f"cat {' '.join(filenames)} > {new_filename}"
		print(cmd)
		#process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command 
		#process.wait()

if __name__ == "__main__":
	main()