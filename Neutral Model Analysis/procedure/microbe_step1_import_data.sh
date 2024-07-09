#!/bin/bash
##Select what type of fastq data to import

echo How much should be trimmed for the left?
read triml
echo How much should be trimmed for the right?
read trimr
echo How many CPUs should be used to run this quality control? Leave yourself atleast 1-2
read speed
echo Please "provide a mapping file for visualization"
read mapping

mkdir working_files
echo Running qiime import tool and creating qza file
#qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path data    \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path working_files/demux.qza

mkdir visual_files
mkdir extra
#Use the output to summarize the data as a .qzv file so we can visualize the results using the qiime visualization website.
#qiime demux summarize \
 --i-data working_files/demux.qza \
 --o-visualization visual_files/demux.qzv

#Use the dada2 denosing command to trim, truncate your sequences.
echo Using Dada deposise paired tool

#qiime dada2 denoise-paired \
 --i-demultiplexed-seqs working_files/demux.qza \
 --p-trunc-len-f $triml \
 --p-trunc-len-r $trimr \
 --p-n-threads $speed \
 --o-table working_files/table.qza \
 --o-representative-sequences working_files/rep-seqs.qza \
 --o-denoising-stats extra/denoising-stats.qza

## Check the merged output with new mapping file
qiime feature-table summarize \
--i-table working_files/table.qza \
--o-visualization visual_files/table.qzv \
--m-sample-metadata-file $mapping

# Check merged rep-seq files
qiime feature-table tabulate-seqs \
--i-data working_files/rep-seqs.qza \
--o-visualization visual_files/rep-seqs.qzv

echo complete



