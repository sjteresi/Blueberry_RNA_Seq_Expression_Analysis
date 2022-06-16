############################################
# Install the required packages, please reference README and
# the version control
# NOTE each installation command must be done separately
#install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages(c("dplyr","tidyverse"))

############################################
# Load the libraries
# TODO semi-automate path

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))

# TODO semi-automate path
# NOTE this doesn't need to be run everytime, only if I install new stuff
############################################
# MASTER ARGS
input_counts_path = '/home/scott/Documents/Uni/Research/Projects/Blueberry_RNA_Seq_Expression_Analysis/results/All_Counts_Blueberry.tsv'
output_dir = '/home/scott/Documents/Uni/Research/Projects/Blueberry_RNA_Seq_Expression_Analysis/results/EdgeR_Differential_Expression'
############################################

# Loading the Data Into R (Part 1):
# Here we will load all of our data from HTSeq into R.
# Below I import the data, perform the appropriate merges and reorient the data.
load_counts_as_table = function(count_file){
	count_table = read.csv(count_file, sep='\t', header=TRUE)
	
	# Remove the extraneous 5 rows at the bottom,
  # these are supplementary info from HTSeq
	count_table = head(count_table, n = -5)

	# Set the index
	rownames(count_table) = count_table$Gene

	# Get rid of gene name column because it is now the index
	count_table = subset(count_table, select = -c(Gene))

	# Remove the rows that are completely 0s, uninformative rows
	count_table = count_table[rowSums(count_table[])>0,]
return(count_table)
}

# Loading the Data Into R (Part 2):
# Here we will filter the data and add the appropriate "metadata"
# so that we can easily recognize each sample.
# I write these as a function so that we can easily utilize them later inside the individual sample comparison chunks.

clean_data = function(Blueberry){
	Metadata_Blueberry = Blueberry
	# Add two empty rows at the end of the data frame to be filled with the
	# experimental factors that we later plug in.
	Metadata_Blueberry[(nrow(Metadata_Blueberry) + 1):(nrow(Metadata_Blueberry) + 2), ] = NA

	# Loop through the columns and assign experimental factors based on the
	# sample names, filling the last two rows
	#columns = colnames(Metadata_Blueberry)
	dra_front_str = 'Dra_'
	lib_front_str = 'Lib_'
	stop_control_str = 'dpo_c'
	stop_treatment_str = 'dpo_t'
	 
	for (i in 1:ncol(Metadata_Blueberry)) {
	  for(xi in 1:7) {
	    xi = toString(xi)
			dra_full_C = paste(dra_front_str, xi, stop_control_str, sep = '')
			dra_full_T = paste(dra_front_str, xi, stop_treatment_str, sep = '')
			lib_full_C = paste(lib_front_str, xi, stop_control_str, sep = '')
			lib_full_T = paste(lib_front_str, xi, stop_treatment_str, sep = '')
			dra_name = paste('DRA_', xi, sep = '')
			lib_name = paste('LIB_', xi, sep = '')

			if (grepl(dra_full_C, colnames(Metadata_Blueberry)[i])) {	
			      Metadata_Blueberry[(nrow(Metadata_Blueberry) - 1),i] <- dra_name
			      Metadata_Blueberry[(nrow(Metadata_Blueberry)), i] <- "Control"

	      } else if (grepl(dra_full_T, colnames(Metadata_Blueberry)[i])) {
			Metadata_Blueberry[(nrow(Metadata_Blueberry) - 1),i] <- dra_name
			Metadata_Blueberry[(nrow(Metadata_Blueberry)), i] <- "Treatment"  

	      } else if (grepl(lib_full_C, colnames(Metadata_Blueberry)[i])) {
			Metadata_Blueberry[(nrow(Metadata_Blueberry) - 1),i] <- lib_name
			Metadata_Blueberry[(nrow(Metadata_Blueberry)), i] <- "Control"  

	      } else if (grepl(lib_full_T, colnames(Metadata_Blueberry)[i])) {
			Metadata_Blueberry[(nrow(Metadata_Blueberry) - 1),i] <- lib_name
			Metadata_Blueberry[(nrow(Metadata_Blueberry)), i] <- "Treatment"  
	      } 
	  }
	# Update rows with "metadata" on each sample's identity
	  row.names(Metadata_Blueberry)[(nrow(Metadata_Blueberry) - 1) : (nrow(Metadata_Blueberry))] = c("Identifier", "Experiment_Treatment")
	}
return(Metadata_Blueberry)
}

EdgeR_Func = function(Counts, G1, G2, xi, data_input_type, test_type, output_dir, treatment_groupings_same=FALSE, comparison_order='Simple') {

  setwd(output_dir)  
	rownames(Counts) <- Counts$Row.names
	Counts = subset(Counts, select = -c(Row.names))
	if (isTRUE(treatment_groupings_same)) {
	  my_grouping = Counts['Identifier',] 
	} else {
	my_grouping = Counts['Experiment_Treatment',] # Make grouping for treatment groups.
	}
	Counts = Counts[!(rownames(Counts) %in% 'Identifier'),]
	Counts = Counts[!(rownames(Counts) %in% 'Experiment_Treatment'),]

	# Convert to matrix
	Counts = as.matrix.data.frame(Counts)
	x = rownames(Counts)
	# We need the gene names for later
	Gene_Row_Key = data.frame("Num"=rownames(as.data.frame(x)),"Gene"=x)
	rm(x)
	Counts = apply(Counts, 2, as.numeric)

	D = DGEList(Counts, group = unlist(my_grouping))
	D = calcNormFactors(D)
	D_Samples = D$samples
	D = estimateCommonDisp(D)
	D = estimateTagwiseDisp(D)
	if (comparison_order == 'Simple'){
	  Fish_Exact = exactTest(D, pair=c('Control', 'Treatment'))
	}
	else if (comparison_order == 'Complex'){
	  Fish_Exact = exactTest(D)
	}
	topTags = topTags(Fish_Exact)

	# P = 0.05
	simplified_DGE = decideTestsDGE(Fish_Exact, p=0.05, adjust='fdr')


	# Write summary table
	write.table(summary(simplified_DGE), file = paste(G1, xi, '_vs_', G2, xi, '_Summary.txt', sep = ''), quote = F, sep = '\t')

	# Write direction of differential expression table
	simplified_DGE_frame = data.frame(simplified_DGE)
	colnames(simplified_DGE_frame) = 'Direction_Differentially_Regulated'
	row.names(simplified_DGE_frame) = Gene_Row_Key$Gene
	Direction = as_tibble(rownames_to_column(simplified_DGE_frame, var = 'Gene_Name'))
	print('writing')
	write_tsv(Direction, paste(G1, xi, '_vs_', G2, xi, '_Direction.tsv', sep = ''))
}



# Compare C vs T
c_obj_str = function(my_obj) {
  deparse(substitute(my_obj))
}

run_comparisons = function(data_input_type, output_dir) {
	print('Running comparisons now!')
	dra_front_str = 'Dra_'
	lib_front_str = 'Lib_'
	stop_control_str = 'dpo_c'
	stop_treatment_str = 'dpo_t'
	 
	for(xi in 1:7) {
	  xi = toString(xi)
		dra_full_C = paste(dra_front_str, xi, stop_control_str, sep = '')
		dra_full_T = paste(dra_front_str, xi, stop_treatment_str, sep = '')
		lib_full_C = paste(lib_front_str, xi, stop_control_str, sep = '')
		lib_full_T = paste(lib_front_str, xi, stop_treatment_str, sep = '')
		dra_name = paste('DRA_', xi, sep = '')
		lib_name = paste('LIB_', xi, sep = '')
		
	  DRA_C =  select(Metadata_Blueberry, matches(dra_full_C))
	  DRA_T =  select(Metadata_Blueberry, matches(dra_full_T))
	  LIB_C =  select(Metadata_Blueberry, matches(lib_full_C))
	  LIB_T =  select(Metadata_Blueberry, matches(lib_full_T))
	 
	  # DRA_C vs DRA_T 
	  Counts = merge(DRA_C, DRA_T, by = 'row.names')
	  G1 = c_obj_str(DRA_C)
	  G2 = c_obj_str(DRA_T)
	  EdgeR_Func(Counts, G1, G2, xi, data_input_type, 'fdr', output_dir, comparison_order='Simple')
	  
	  # LIB_C vs LIB_T 
	  Counts = merge(LIB_C, LIB_T, by = 'row.names')
	  G1 = c_obj_str(LIB_C)
	  G2 = c_obj_str(LIB_T)
	  EdgeR_Func(Counts, G1, G2, xi, data_input_type, 'fdr', output_dir, comparison_order='Simple')
	  
	  Counts = merge(DRA_T, LIB_T, by = 'row.names')
	  G1 = c_obj_str(DRA_T)
	  G2 = c_obj_str(LIB_T)
	  EdgeR_Func(Counts, G1, G2, xi, data_input_type, 'fdr', output_dir, treatment_groupings_same = TRUE, comparison_order='Complex')

	  Counts = merge(DRA_C, LIB_C, by = 'row.names')
	  G1 = c_obj_str(DRA_C)
	  G2 = c_obj_str(LIB_C)
	  EdgeR_Func(Counts, G1, G2, xi, data_input_type, 'fdr', output_dir, treatment_groupings_same = TRUE, comparison_order='Complex')
	}
}

Blueberry = load_counts_as_table(input_counts_path)
Metadata_Blueberry = clean_data(Blueberry)
run_comparisons('All', output_dir)

