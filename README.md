# Blueberry RNA Seq
This is a repository containing the code for generating gene expression and differential gene expression data from RNA-Seq reads in blueberry.

# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/huckleberry-hound) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |


This README will be divided into two parts, generating counts, and performing differential expression analysis.

# Part 1 - Generate Counts:
# TODO General description

## Example Directory Structure:
```
/Input
/Scripts
/EdgeR_Output

/Input/Count/
/Input/Count/Collate
/Input/Genome_Index
/Input/Raw_Input

/Scripts/CountCollate
/Scripts/EdgeR
/Scripts/Genome_Index
/Scripts/HTSeq
/Scripts/Mapping
/Scripts/Trim

/EdgeR_Output/Bonferroni
/EdgeR_Output/FDR
```

# Workflow Plan:
1. Index the genome with STAR
2. Trim the reads
3. Align the reads with STAR
4. Quantify the reads with HTSeq
5. Collate the reads
6. Run EdgeR to get differential expression

# 1) Genome Indexing:
The first step in our workflow is indexing the genome. This step is necessary to improve the speed and functionality of the mapping algorithms in part **3** and **4**. Thankfully we can use **STAR** to perform the genome indexing, which is the same program we use to align the reads in part **3**. 

This step requires two files, which I place in `/Input/Raw_Input`: a `gtf` or `gff` of the annotated genes, and a `fasta` file from the genome assembly. It is important to note whether you have a gtf or gff annotation file because that will impact an option in the genome index script. Check the STAR manual for more notes (section 2.2.3 for STAR v 2.7.3a).

This command, `/Scripts/Genome_Index/GenomeIndex_STAR.sb` runs rather quick. It outputs several genome index files into my designated output folder `/Input/Genome_Index`. Examine the code and manual for more notes.

# 2) Trim the Reads:
The second step in our workflow is trimming the reads. This script is one of the more involved scripts, and you may have to modify the `R1` and `R2` variables quite a bit to get the names matching your naming scheme.

# 3) Align the Reads with STAR:
The third step in the workflow involved using STAR to map the reads to the reference genome. It is important to have the right reference genome for your organism, otherwise things may not align in the best possible way. Alternative programs to STAR include TopHat and HISAT. 

This script, like the read trimming, is a more involved script and you may have to edit it to match the naming scheme of your reads. However, STAR itself is relatively straightforward. The only main option that I edited on the STAR call is that I set `--outFilterMultimapNmax 1`. I did this option because we wanted to eliminate multiple alignments for a read, and avoid statistical anomalies for this set. 

I then run STAR over 8 threads and specify the inputs and outputs.
## To Do
2. Add comments to EdgeR script
3. Generate workflow explanation or diagram
