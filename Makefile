# Script to recreate analyses
# __file__ Makefile
#  __author__ Scott Teresi, Alder Fulton

ROOT_DIR:= $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA:= $(ROOT_DIR)/..Blueberry_Data
DEV_NEW_DATA:= $(ROOT_DIR)/data
DEV_RESULTS:= $(ROOT_DIR)/results
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
	# NOTE ran on personal computer not cluster
	# NOTE I hate package management in R. It is so poorly implemented.
	# I cannot easily get a version greater than R 3.6 in Anaconda
	# So Anaconda is out of the question
	# So I am stuck using this crappy renv. 
	# renv doesn't like me trying to declare paths whenever
	# I have to load or unload packages
	# So we have to run EdgeR via some shitty method.
	# You CANNOT run via Makefile. These notes here are just for reference
	# You need to CD into the requirments folder
	# And then run the command manually.
	#
	# ADDENDUM, run this in Rstudio because it is so hard from command line.
	@echo Running EdgeR
	cd $(ROOT_DIR)/requirements
	Rscript src/EdgeR/EdgeR_Blueberry.R $(DEV_RESULTS)/All_Counts_Blueberry.tsv $(DEV_RESULTS)/EdgeR_Differential_Expression
