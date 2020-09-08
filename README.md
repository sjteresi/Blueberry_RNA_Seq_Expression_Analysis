# Blueberry RNA Seq
Generate lists of differentially expressed genes from RNA-seq data in Blueberry.

# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/huckleberry-hound) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |

# Workflow Plan:
1. Index the genome with STAR
2. Trim the reads
3. Align the reads with STAR
4. Quantify the reads with HTSeq
5. Collate the reads
6. Run EdgeR to get differential expression

# Follow Up Project:
This project continues with [Network Analysis](https://github.com/EdgerLab/Blueberry_Network_Rewiring). There I identify orthologs and examine gene network differences between the two blueberry cultivars in this dataset.
