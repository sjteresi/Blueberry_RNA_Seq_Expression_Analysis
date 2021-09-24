# Script to recreate analyses
# __file__ Makefile
#  __author__ Scott Teresi, Alder Fulton

ROOT_DIR:= $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA:= $(ROOT_DIR)/..Blueberry_Data
DEV_CODE_FILES:=$(ROOT_DIR)/src

# Step 1 of Workflow:
# Index genome with STAR
index_genome:
	# NOTE this needs to be run on the cluster
	@echo Running STAR
	sbatch $(DEV_CODE_FILES)/Genome_Index/GenomeIndex_STAR.sb

# Step 2 of Workflow:
# Trim reads for STAR mapping
trim_reads:
	# NOTE this must be run on the cluster
	@echo Trimming reads
	sbatch $(DEV_CODE_FILES)/Trim/Trim_All.sb

# Step 3 of Workflow:
# Align reads with STAR
align_reads:
	# NOTE this must be run on the cluster
	@echo Aligning reads
	sbatch $(DEV_CODE_FILES)/Mapping/STAR_Map.sb

# Step 4 of Workflow:
# Quantify reads with HTSeq
quantify_reads:
	# NOTE this must be run on the cluster
	@echo Quantifying reads with HTSeq
	sbatch $(DEV_CODE_FILES)/HTSeq/HTSeq.sb

# Step 5 of Workflow:
# Collate reads with Python into one table
collate_reads:
	# NOTE this must be run on the cluster
	@echo Collating count files to make one large count file
	python $(DEV_CODE_FILES)/CountCollate/count_collate.py $(DEV_DATA)/Counts/

# Step 6 of Workflow:
# Run EdgeR to determine differentially expressed genes
EdgeR:
	# NOTE I ran this on personal computer, file paths are hardcoded into the script
	@echo Running EdgeR
	Rscript src/EdgeR/EdgeR_Blueberry.Rmd
