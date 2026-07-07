from ncls import NCLS

import multiprocessing as mp
from functools import partial
import pandas as pd

def find_cpg_overlaps_chr(df_regions: pd.DataFrame, df_cpgs: pd.DataFrame) -> pd.DataFrame:
	'''
	Find CpGs located in the given regions (in a chromosome)
	'''

	object_ncls = NCLS(df_regions['start'].values-1, df_regions['end'].values, df_regions.index.values)
	cpg_idx, region_idx = object_ncls.all_overlaps_both(df_cpgs["start"].values,
														df_cpgs["end"].values,
														df_cpgs.index.values)

	# When a CpG falls inside multiple regions, keep only its first match.
	# This mirrors UXM's behaviour (pat_utils.py: "only the first matching region is used").
	if cpg_idx.shape[0] != len(set(cpg_idx)):
		chrom = df_regions["chr"].iloc[0]
		print(f"Warning: some CpGs overlap multiple regions on {chrom}; "
			  f"keeping first match per CpG (consistent with UXM).")
		seen = set()
		keep = []
		for i, c in enumerate(cpg_idx):
			if c not in seen:
				seen.add(c)
				keep.append(i)
		cpg_idx    = cpg_idx[keep]
		region_idx = region_idx[keep]

	df_cpgs = df_cpgs.loc[list(cpg_idx),:]
	df_cpgs["dmr_label"] = region_idx
	df_cpgs["dmr_ctype"] = [df_regions.loc[ii, "ctype"] for ii in region_idx]
	df_cpgs["dmr_coordinates"] = [f"{df_regions.loc[ii, 'chr']}:{df_regions.loc[ii, 'start']}-{df_regions.loc[ii, 'end']}"
								   for ii in region_idx]
	return df_cpgs

def find_cpg_overlaps(df_region: pd.DataFrame, df_cpg: pd.DataFrame, n_cores=10) -> pd.DataFrame:
	'''
	Find CpGs overlapping with the given region
	'''

	list_args = list()
	chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
	for chrom in chroms:
		q = df_region[df_region["chr"] == chrom]
		s = df_cpg[df_cpg["chr"] == chrom]
		if (q.shape[0] == 0) or (s.shape[0] == 0):
			continue
		else:
			list_args.append((q, s))

	with mp.Pool(n_cores) as pool:
		overlaping_cpgs = pool.starmap(find_cpg_overlaps_chr, list_args)

	return pd.concat(overlaping_cpgs)

def merge_close_regions(df_region: pd.DataFrame, distance: int = 150) -> pd.DataFrame:
	"""
	Merge regions very close to each other (<150 bps) when they are from the same cell-type DMRs 
	Assumed that all regions are on the same chromosome
	"""
	df_merged_region = list()

	new_region = df_region.loc[df_region.index[0], ["chr", "start", "end", "startCpG", "endCpG", "ctype"]].copy()
	for i in range(1, df_region.shape[0]):
		if (df_region.loc[df_region.index[i], "start"] - new_region["start"] < distance) and \
		   (df_region.loc[df_region.index[i], "ctype"] == new_region["ctype"]):
			new_region["end"] = df_region.loc[df_region.index[i], "end"]
			new_region["endCpG"] = df_region.loc[df_region.index[i], "endCpG"]
		else:
			# if regions are for different cell types, they cannot be merged
			df_merged_region.append(new_region)
			del new_region
			new_region = df_region.loc[df_region.index[i], ["chr", "start", "end", "startCpG", "endCpG", "ctype"]].copy()
	return pd.DataFrame(df_merged_region).reset_index(drop=True)
