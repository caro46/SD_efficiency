#!/bin/env python3

import os, sys, argparse, subprocess


def splitText(name):
	""" Splits a name ending in digits into two strings """
	# Iterate through each potential position
	for i in range(1, len(name)):
		# Split string into first and second sections
		first = name[:i]
		second = name[i:]
		# Try to convert second section into an integer
		try:
			# If success, return first and second sections
			test = int(second)
			return first, second
		except ValueError:
			# If fail, thry next iteration
			pass
	# Raise ValueError if no valid split is detected
	raise ValueError


def main():

	## arguments with help messages
	parser = argparse.ArgumentParser(description="Running kallisto (indexing ot quantification) on paired-end rnaseq reads", allow_abbrev=False)
	parser.add_argument("--ref_transcriptome", type = str, default = "", help = "Transcriptome to index (assuming .gz). Default: No")
	parser.add_argument("--fastq", nargs = "*", help = "Fastq file(s) to use. Default: No", default = []) #default = empty list
	parser.add_argument("--Read1_info", type = str, help = "String distinguishing Read1 and Read2 files for Paired-end reads. Ex: _R1. Default: _R1", default = "_R1")
	parser.add_argument("--Read2_info", type = str, help = "String distinguishing Read1 and Read2 files for Paired-end reads. Ex: _R2. Default: _R2", default = "_R2")
	parser.add_argument("--single_read_info", type = str, help = "String to exclude unpaired reads. Ex. single", default = "single")
	parser.add_argument("--threads", type = int, help = "Number of threads. Default: 2", default = 2)
	parser.add_argument("--output_directory", type = str, default = "", help = "Output directory for kallisto quant outputs files. Default: No")
	parser.add_argument("--prefix_out", type = str, nargs = "?", help = "Additional info to be included in the output name. Example: Run1_ Default: empty", const = "", default = "")
	parser.add_argument("--step", type = str, help = "Step to be run. Ex: index, quant. 'index' will produce a [ref_transcriptome].idx file, 'quant' assumes index is [ref_transcriptome].idx. Default: quant", default = "quant")
	parser.add_argument("--bootstrap", type = int, help = "Number of bootstrap samples. Default: 100", default = 100)

	## Parsing the arguments
	args = parser.parse_args()
	# print("Arguments: ",args,"\n")

	if args.step == "index":
		ref_no_gz = os.path.splitext(args.ref_transcriptome)[0]
		ref_no_ext = os.path.splitext(ref_no_gz)[0]
		cmd = f"""kallisto index -i {ref_no_ext}.idx {args.ref_transcriptome}"""
		print(f"Command: {cmd}")
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command 
		process.wait()

	elif args.step == "quant":
		ref_no_gz = os.path.splitext(args.ref_transcriptome)[0]
		ref_no_ext = os.path.splitext(ref_no_gz)[0]
		for filename in args.fastq:
			if filename.endswith(".fastq.gz") | filename.endswith(".fq.gz"):
				if args.single_read_info not in filename:
					if args.Read1_info in filename:
						#print(f"{filename}")
					## Extract individual IDs
						file_base=os.path.basename(filename)
						file_no_gz = os.path.splitext(file_base)[0]
						file_no_ext = os.path.splitext(file_no_gz)[0]
						# Handle bad formatting case
						try:
							Gene, Individual, Sample_numb, Lane, Read, trim, junk = file_no_ext.split("_")
						except ValueError:
							gene_individual, Sample_numb, Lane, Read, trim, junk = file_no_ext.split("_")
							Gene, Individual = splitText(gene_individual)
						R1_filename = filename
						R2_filename = R1_filename.replace(args.Read1_info, args.Read2_info)
						cmd = f"""kallisto quant --thread {args.threads} -b {args.bootstrap} -i {ref_no_ext}.idx -o {args.output_directory}/{Gene}_{Individual}_{args.prefix_out}{Lane}_kallisto_bout_out {R1_filename} {R2_filename}"""
						process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True) #execute the command 
						process.wait()
						print(f"Command: {cmd}")
						print(f"Done with Individual: {Gene}_{Individual}")

if __name__ == "__main__":
	main()
