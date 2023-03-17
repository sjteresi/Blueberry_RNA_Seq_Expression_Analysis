#!/bin/bash -login

#SBATCH -J HTSeq
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-84
#SBATCH --output=SubmissionScripts/ArrayHTSeq.%A_%a.log
#SBATCH --output=/mnt/research/edgerpat_lab/Scotty/Blueberry_RNA_Seq_Expression_Analysis/results/htseq/htseq_logs/ArrayHTSeq.%A_%a.log
#SBATCH -A general


# Load Modules
module purge
module load Anaconda2/4.2.0
source activate HTSeq

# BODY
begin=`date +%s`
echo $HOSTNAME
echo "My Task ID:" $SLURM_ARRAY_TASK_ID
mapping_location='/mnt/research/edgerpat_lab/Scotty/Blueberry_RNA_Seq_Expression_Analysis/results/star_map/'
cd $mapping_location

Read_Folder_Path=$(ls -d $PWD/* | sed -n ${SLURM_ARRAY_TASK_ID}p)
Read_Base=`basename $Read_Folder_Path` 

cd /mnt/research/edgerpat_lab/Scotty/Blueberry_RNA_Seq_Expression_Analysis/results/htseq/

# Code
htseq-count --format=sam --order=name --stranded=yes --type=exon --idattr=Parent --mode=union ${Read_Folder_Path}/SortedSam.sam /mnt/research/edgerpat_lab/Scotty/Blueberry_RNA_Seq_Expression_Analysis/data/V_corymbosum_v1.0_geneModels.gff > ${Read_Base}_Results.count
gzip ${Read_Folder_Path}/SortedSam.sam


end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
