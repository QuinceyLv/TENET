# TENET_single_core_version

A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data.  
Modified from https://github.com/neocaleb/TENET

## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014

## Dependency

	python3
	openmpi (>4.0)
	JPype
	numpy
	statsmodels
	scanpy

This script could run with python==3.10.4, openmpi==4.1.3, JPype==1.3.0, numpy==1.21.6, statsmodels==0.13.2, and scanpy==1.9.1.

#### Usage

With h5ad file (PAGA output) as input:

	python TENET.py \  
		--h5ad h5ad_file \  
		--out out_dir

With expression_file, trajectory_file and cell_select_file as input:

		python TENET.py \  
		--filelist expression_data.csv trajectory.txt cell_select.txt \
		--out out_dir

Optional parameter:

	--hist history_length  Set history length, default=1

#### Input format

###### (1) expression_file (raw count is recommended) - a csv file with N cells in the rows and M genes in the columns.

		GENE_1	GENE_2	GENE_3	...	GENE_M

	CELL_1	

	CELL_2

	CELL_3

	.
	.
	.

	CELL_N

###### (2) trajectory_file - a text file of pseudotime data with N time points in the same order as the N cells of the expression file.

	0.098
	0.040
	0.023
	.
	.
	.
	0.565

###### (3) cell_select_file - a text file of cell selection data with N Boolean (1 for select and 0 for non-select) data in the same order as the N cells of the expression file.

	1
	1
	0
	.
	.
	.
	1

#### Output

	TE_result_matrix.txt - TEij, M genes x M genes matrix representing the causal relationship from GENEi to GENEj.

	TE	GENE_1	GENE_2	GENE_3	...	GENE_M
	GENE_1	0	0.05	0.02	...	0.004
	GENE_2	0.01	0	0.04	...	0.12
	GENE_3	0.003	0.003	0	...	0.001
	.
	.
	.
	GENE_M	0.34	0.012	0.032	...	0

## Downstream analysis

#### (1) Reconstructing GRN
###### Usage
	python makeGRN.py [cutoff for FDR]
	python makeGRNsameNumberOfLinks.py [number of links]
	python makeGRNbyTF.py [species] [cutoff for FDR]
	python makeGRNbyTFsameNumberOfLinks.py [species] [number of links]
	** Note that "TE_result_matrix.txt" should be in the same folder.

###### Example
	python makeGRN.py 0.01
	python makeGRNsameNumberOfLinks.py 1000
	python makeGRNbyTF.py human 0.01
	python makeGRNbyTFsameNumberOfLinks.py human 1000

###### Output file
	TE_result_matrix.fdr0.01.sif
	TE_result_matrix.NumberOfLinks1000.sif
	TE_result_matrix.byGRN.fdr0.01.sif
	TE_result_matrix.byGRN.NumberOflinks1000.sif

###### Parameter
	[cutoff for fdr] - A cutoff value for FDR by z-test
	[number of links] - The number of links of the GRN
	[species] - User can choose [human/mouse/rat]

#### (2) Trimming indirect edges
###### Usage
	python trim_indirect.py [name of GRN] [cutoff]
###### Example
	python trim_indirect.py TE_result_matrix.fdr0.01.sif 0
###### Output file
	TE_result_matrix.fdr0.01.trimIndirect0.0.sif
###### Parameter
	[cutoff] - A cutoff value for trimming indirect edges. Recommended range is -0.1 to 0.1

#### (3) Counting out-degree of a given GRN
###### Usage
	python countOutdegree.py [name of GRN]
###### Example
	python countOutdegree.py TE_result_matrix.fdr0.01.sif
###### Output file
	TE_result_matrix.fdr0.01.sif.outdegree.txt
