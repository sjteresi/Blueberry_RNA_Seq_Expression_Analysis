# Script to recreate analyses
# __file__ Makefile
#  __author__ Scott Teresi, Alder Fulton

ROOT_DIR:= $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_RESULTS:= $(ROOT_DIR)/results

.PHONY: setup
setup:
	mkdir -p doc/ results/ data/ src/ requirements/

# Step 1 of Workflow:
# Index genome with STAR
.PHONY: index_genome
index_genome:
	# NOTE this needs to be run on the cluster
	@echo Running STAR
	sbatch $(ROOT_DIR)/src/genome_index_STAR.sb

# Step 2 of Workflow:
# Trim reads for STAR mapping
.PHONY: trim_reads
trim_reads:
	# NOTE this must be run on the cluster
	@echo Trimming reads
	sbatch $(ROOT_DIR)/src/trim_all.sb

# Step 3 of Workflow:
# Align reads with STAR
.PHONY: align_reads
align_reads:
	# NOTE this must be run on the cluster
	@echo Aligning reads
	sbatch $(ROOT_DIR)/src/STAR_map.sb

# Step 4 of Workflow:
# Quantify reads with HTSeq
.PHONY: quantify_reads
quantify_reads:
	# NOTE this must be run on the cluster
	@echo Quantifying reads with HTSeq
	sbatch $(ROOT_DIR)/src/HTSeq.sb

# Step 5 of Workflow:
# Collate reads with Python into one table
.PHONY: collate_reads
collate_reads:
	# NOTE this must be run on the cluster
	@echo Collating count files to make one large count file
	python $(ROOT_DIR)/count_collate.py $(DEV_RESULTS)/htseq/ -o $(DEV_RESULTS)/count_collate

# Step 6 of Workflow:
# Run EdgeR to determine differentially expressed genes
.PHONY: EdgeR
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
	Rscript src/EdgeR/EdgeR_Blueberry.R $(DEV_RESULTS)/count_collate/All_Counts_Blueberry.tsv $(DEV_RESULTS)/EdgeR_Differential_Expression
