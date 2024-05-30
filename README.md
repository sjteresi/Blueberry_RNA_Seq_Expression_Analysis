# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/sjteresi) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | NA               | <edgerpat@msu.edu> |

# Purpose:
Generate gene expression data for a blueberry genome, to be used in a genomics project later.
- Align RNA-Seq data to said genome
- Generate gene expression tables (FPKM, TPM)
- Generate lists of differentially expressed genes (DEGs).

# Solution:
The workflow is recapitulated in the `Makefile`.
1. Index the genome with STAR
2. Trim the reads
3. Align the reads with STAR
4. Quantify the reads with HTSeq
5. Collate the reads into one table of counts
6. Run EdgeR to determine differentially expressed genes

## Indexing the genome with STAR
First, the genome was indexed using `STAR v2.6.1`, the genome FASTA file, and the gene annotation file.
The FASTA file and annotation file were derived from the [genome publication](https://academic.oup.com/gigascience/article/8/3/giz012/5304886).
The commands used to perform this can be found in the `src/genome_index_STAR.sb` script.

## Trimming the reads:
Illumina adapters were removed from the raw reads using `Trimmomatic v0.38`. More details can be found in the `src/trim_all.sb` script.

## Alignment of reads:
Filtered reads were then aligned to the genome using `STAR v2.6.1` and the script associated with this command may be found at `src/STAR_map.sb`. Multimapping reads were discarded.

## Quantification of counts and calculation of differentially expressed genes:
Count files were calculated using `HTSeq v0.12.4`.
Individual count files were then collated using the custom Python script at `src/count_collate.py`.
This was then used as input to `EdgeR v3.30.3` `(R v4.0.2)` to determine which genes are differentially expressed in each condition comparison.
An FDR correction using the Benjamini-Hochberg method was utilized. The script associated with this analysis may be found at `src/EdgeR/EdgeR_Blueberry.R`. 

# Follow Up Project:
This project continues with [Network Analysis](https://github.com/sjteresi/Blueberry_Network_Rewiring).
There I identify orthologs and examine gene network differences between the two blueberry cultivars in this dataset.

# Version Control:
Refer to the `requirements/` folder. Major packages used are STAR v2.6.1, Trimmomatic v0.38, SAMtools v1.9, R v4.0.2, edgeR_3.30.3, and HTSeq v0.12.4.

# Context and Future Considerations:
This was my first experience working on an RNA-Seq project, and I started this during the first year of my PhD.
I learned a lot about operating on the computing cluster, in particular running array jobs.
If I had the chance to go back and do it again, I would modify the array jobs to read from a manifest file, where each row would be a job.
I think this would be more legible and reproducible than `ls`-ing a directory and feeding that to the job with `sed`; the manifest file could easily be tracked with Git.

I also think the EdgeR script is a little too hard-coded, and would be a pain to refactor or use for a similar project.
If I had the change to re-do the DEG analysis, I would also consider generalizing the EdgeR script to operate on one pair of gene expression columns at a time, reading those inputs from files already saved on disk.
Most of the work (and a major source of hard-coding) was the product of just trying to iterate over my gene expression table, get the name for the comparison, and then perform the comparison (which is pretty simple on its own from a code standpoint).
