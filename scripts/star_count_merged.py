#!/bin/env python3

import os, sys, argparse, subprocess
from collections import defaultdict

def main():

	## arguments with help messages
	parser = argparse.ArgumentParser(description="Combine and convert star count outputs for R analysis.", allow_abbrev=False)
	parser.add_argument("--files", nargs = "*", help = "File(s) to merge. Default: No")
	parser.add_argument("--output_directory", type = str, default = "", help = "Output directory. Default: No")
	parser.add_argument("--prefix_out", type = str, help = "Output prefix. Default: No")

	## Parsing the arguments
	args = parser.parse_args()
	dictionary_filenames = defaultdict(dict)
	gene_set = set() #create an empty set to store the gene names
	file_names = []
	for file in args.files:
		file_base = os.path.basename(file)
		file_base = file_base.split("ReadsPerGene.out.tab")[0]
		file_names.append(file_base)
		# Get data from file name
		#gene, individual, Lane_Star_ext = file_base.split("_")
		#lane = Lane_Star_ext.split("ReadsPerGene.out.tab")[0]
		dictionary_filenames[file_base] = {} #value is empty dictionary
		## Populate the empty dictionary
		for i, line in enumerate(open(file)):
			if i in range(4): #(inclusive, exclusive) - default inclusive = 0
				pass
			else:
				gene, count_unstranded_RNAseq, count_1st_read_strand, count_2nd_read_strand = line.rstrip('\n').split("\t")
				key_gene = gene
				dictionary_filenames[file_base][key_gene] = count_unstranded_RNAseq
				# make a set of genes (works like a list except uses .add instead of .append)
				# keep track of every gene (no duplicate)
				gene_set.add(gene)

	## Iterate through the genes
	# Assume file_names is order of columns
	#file_names = []
	outfile = open(args.output_directory + "/" + args.prefix_out + ".tab","w")
	outfile.write("Gene_name" + "\t"+ "\t".join(file_names) + "\n")
	#outfile.write("\n")
	for gene in gene_set:
		#print(gene, end='\t')
		outfile.write(gene + "\t")
		for file in file_names:
			if gene in dictionary_filenames[file]:
				#print(dictionary_filenames[file][gene], end='\t')
				outfile.write(dictionary_filenames[file][gene] + "\t")
			else:
				#print("NA", end='\t')
				outfile.write("NA" + "\t")
		#print(end='\n')
		outfile.write("\n")

	#outfile.write('\t'.join(count_unstranded_RNAseq))
	#outfile.close()

if __name__ == "__main__":
	main()