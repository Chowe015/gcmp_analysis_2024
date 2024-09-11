#! /bin/bash
# import table and rep-seq file for down stream analysis
echo "Please select a mapping file..."
read mapping
echo "How many cpu would you like to use?" 
read speed

echo filtering rep-seq file using silva extended taxonomy

qiime feature-classifier classify-consensus-vsearch --i-query working_files/rep-seqs.qza --i-reference-reads silva_extended_sequences.qza --i-reference-taxonomy silva_extended_taxonomy.qza \
 --p-threads $speed \
 --o-classification working_files/silva_extended_classification_taxonomy.qza \
 --o-search-results extra/silva_extended_classification_results.qza \

#filter out non-target reads from otu table 

echo filtering table...

qiime taxa filter-table  \
 --i-table working_files/table.qza \
 --i-taxonomy extra/silva_extended_classification_taxonomy.qza \
 --p-exclude 'mitochondria,chloroplast,eukaryota' \
 --o-filtered-table working_files/nomito_table.qza \

echo filtering rep_seq

qiime taxa filter-seqs \
 --i-sequences working_files/rep-seqs.qza \
 --i-taxonomy extra/silva_extended_classification_taxonomy.qza \
 --p-exclude 'mitochondria,chloroplast,eukaryota' \
--o-filtered-sequences working_files/nomito_rep-seqs.qza

echo Creating no mitochrondrial table visualization file
## Visualize mitochondrial removal step results
qiime feature-table summarize \
 --i-table working_files/nomito_table.qza \
 --o-visualization visual_files/nomito_table.qzv \
 --m-sample-metadata-file $mapping\

#echo Creating visualization files for filtered OTU table, rep-seqs and taxonomy files

# check taxonomy file after mitochondrial removal step
qiime metadata tabulate \
 --m-input-file working_files/silva_extended_classification_taxonomy.qza \
 --o-visualization visual_files/nomito_taxonomy.qzv

echo "Building tree from filtered req-seq file which can take several hours"

## Build Phylogeny of merged Carib files
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences working_files/nomito_rep-seqs.qza \
  --p-n-threads $speed \
  --o-alignment extra/aligned_nomito_rep-seqs.qza \
  --o-masked-alignment extra/masked-alignment_rep-seqs.qza \
  --o-tree extra/unrooted-tree.qza \
  --o-rooted-tree working_files/nomito_rooted-tree.qza

echo "Output Taxonomy and Tree as tsv andd .nwk format"
qiime tools export \
 --input-path working_files/nomito_rooted-tree.qza \
 --output-path ./

qiime tools export \
 --input-path working_files/silva_extended_classification_taxonomy.qza \
 --output-path ./
 
qiime taxa barplot \
 --i-table working_files/nomito_table.qza \
 --i-taxonomy working_files/silva_extended_classification_taxonomy.qza \
 --m-metadata-file $mapping \
 --o-visualization visual_files/nomito_taxa-bar-plots.qzv


echo "Complete"
