# Goal:
Align RNA-Seq data to a genome, generate gene expression tables (FPKM, TPM), generate lists of differentially expressed genes

# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/sjteresi) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |

# Workflow Plan:
1. Index the genome with STAR
2. Trim the reads
3. Align the reads with STAR
4. Quantify the reads with HTSeq
5. Collate the reads into one table of counts
6. Run EdgeR to determine differentially expressed genes

# Description:
## Indexing the genome with STAR
First, the genome was indexed using `STAR v2.6.1`, the genome FASTA file, and the gene annotation file. The FASTA file and annotation file were derived from the [genome publication](https://academic.oup.com/gigascience/article/8/3/giz012/5304886). The commands used to perform this can be found in the `src/GenomeIndex_STAR.sb` script.

## Trimming the reads:
Illumina adapters were removed from the raw reads using `Trimmomatic v0.38`. More details can be found in the `src/Trim/Trim_All.sb` script.

## Alignment of reads:
Filtered reads were then aligned to the genome using `STAR v2.6.1` and the script associated with this command may be found at `src/Mapping/STAR_Map.sb`. Multimapping reads were discarded.

## Quantification of counts and calculation of differentially expressed genes:
Count files were calculated using `HTSeq v0.12.4`. Individual count files were then collated using the custom Python script at `src/CountCollate/count_collate.py`. This was then used as input to `EdgeR v3.30.3` `(R v4.0.2)` to determine which genes are differentially expressed in each condition comparison. An FDR correction using the Benjamini-Hochberg method was utilized. The script associated with this analysis may be found at `src/EdgeR/EdgeR_Blueberry.Rmd`. 

# Follow Up Project:
This project continues with [Network Analysis](https://github.com/EdgerLab/Blueberry_Network_Rewiring). There I identify orthologs and examine gene network differences between the two blueberry cultivars in this dataset.

# Version Control:
Refer to the `requirements/` folder. Major packages used are STAR v2.6.1, Trimmomatic v0.38, SAMtools v1.9, R v4.0.2, edgeR_3.30.3, and HTSeq v0.12.4.
