import argparse, re, os, random, gc, pathlib, glob,json
import pandas as pd
import numpy as np

from read import simulate_reads
from region import merge_close_regions, find_cpg_overlaps
from database import download

def arg_parser():
	# Parse input arguments
	parser = argparse.ArgumentParser()

	parser.add_argument("-o", "--output_dir", type=str, default="./", help="Directory to save the results (default: ./)")
	parser.add_argument("-c", "--cores", type=int, default=1, help="Number of cores for multiprocessing. A larger number increases the computation speed. (default: 1)")
	parser.add_argument("-r", "--f_region", type=str, default=None, help="Selected regions for training MethylBERT. If not given, it automatically selects regions for the given files.")
	parser.add_argument("-g", "--genome", type=str, default="hg19", help="Reference genome (either hg19 or hg38). Currently only hg19 is available. (default: hg19)")
	parser.add_argument("-f", "--f_input", required=True, help="Text file containing a list of .pat files OR path to a .pat file")

	return parser.parse_args()

def filename2ctype(f_name: str) -> str:
	
	'''
	
	Extract cell types from the file name
	It returns the cell type matching to the cell types written in the region information 
	e.g.,) Blood-Monocytes -> Blood-Mono+Macro
	
	'''

	f_name = os.path.basename(f_name)
	file_ctype = f_name.split("_")[1].split("-Z")[0]
	try:
		return df_ctype_match[file_ctype]
	except KeyError:
		print(f"{f_name} cannot be matched to the 39 cell types. NA is assigned for the cell type")
		return "NA"

if __name__=="__main__":
	
	args = arg_parser()

	### Set-ups
		# output directory
	if not os.path.exists(args.output_dir):
		os.mkdir(args.output_dir)

	# load cell tyep match dictionary 
	global df_ctype_match
	f_cell_type_match = download("cell_type_match")
	with open(f_cell_type_match, "r") as fp:
		df_ctype_match = json.load(fp)

	# input files
	f_input_extension = os.path.splitext(args.f_input)[-1]
	if f_input_extension == ".pat":
		ctype = filename2ctype(args.f_input)
		df_files = pd.DataFrame({
			"files": [args.f_input],
			"cell_type": [ctype]
			})
	else:
		# list of files are given 
		df_files = pd.read_csv(args.f_input, header=None)
		if df_files.shape[1] != 1:
			raise ValueError("The input file must have only one column containing .pat file names.")
		df_files.columns = ["files"]
		df_files["cell_type"] = df_files["files"].apply(lambda x : filename2ctype(x))

	# Region selection
	if args.f_region is not None:
		df_region = pd.read_csv(args.f_region, sep=",")
	else:
		f_region = download("unmethyl_regions")
		df_region = pd.read_csv(f_region, sep="\t")
	unique_ctypes = df_files["cell_type"].unique()
	df_region = df_region.loc[df_region["Type"].apply(lambda x: x in unique_ctypes), :]
	df_region["dmr_id"] = list(range(df_region.shape[0]))
	df_region.rename(columns={"Type":"ctype"}, inplace=True)
	print(df_region)
	unique_ctypes=df_region["ctype"].unique()
	print(f"Regions for {unique_ctypes} are selected")

	# Read cpg file 
	if args.genome == "hg19":
		print("Read hg19 cpg files")
		f_cpgs = download("hg19_cpgs")
	else:
		pass
	df_cpg = pd.read_csv(f_cpgs, sep="\t")
	df_cpg = df_cpg.set_index("index")

	# To match 0-based and 1-based
	df_cpg["start"] -= 1
	df_cpg["end"] -= 1

	df_cpg["prev_cpg"] = [-1] + list(df_cpg["start"][:-1]) # cytosine pos of previous CpG
	df_cpg["next_cpg"] = list(df_cpg["start"][1:]) + [-1] # cytosine pos of next CpG

	# Reset columns for finding overlaps
	df_subject = df_cpg.loc[:, ["seqnames", "start", "end", "prev_cpg", "next_cpg"]]
	df_subject.columns = ["chr", "start", "end", "prev_cpg", "next_cpg"]
	df_region.index = df_region["dmr_id"]

	# Find CpGs overlapping with given regions
	print("Find CpGs in the selected regions")
	cpg_overlaps = find_cpg_overlaps(df_region.loc[:, ["chr", "start", "end", "ctype"]], 
									 df_subject, 
									 n_cores = args.cores)
	f_out = os.path.join(args.output_dir, "selected_cpgs.csv")
	cpg_overlaps.to_csv(f_out, header= True, sep="\t", index=True)
	print(cpg_overlaps)
	del df_cpg, df_region, df_subject
	gc.collect()

	for idx in df_files.index:

		# Read input .pat file
		input_file = df_files.loc[idx, "files"]
		print(f"{input_file} is being processed...")
		f_out = os.path.join(args.output_dir, 
							 os.path.basename(input_file).replace(".pat","_reads.csv"))
		df_reads = pd.read_csv(input_file, sep="\t", header=None)
		df_reads.columns = ["chr", "index", "methyl", "n_reads"] # for the consistency 

		# Read simulation 
		res = simulate_reads(df_reads, cpg_overlaps, genome=args.genome, n_cores=args.cores)
		res["ctype"] = df_files.loc[idx, "cell_type"]
		res.to_csv(f_out, header= True, sep="\t", index=False)