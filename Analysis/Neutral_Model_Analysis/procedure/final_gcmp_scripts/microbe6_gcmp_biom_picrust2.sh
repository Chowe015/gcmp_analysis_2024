#!/bin/bash
#SBATCH --job-name=picrust
#SBATCH --nodes=1
#SBATCH --array=1-3
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mem=100G

#import data into fastp
echo "Picrust2 pipeline for GCMP Neutral Psuedo table Analysis..."
module load anaconda
conda activate picrust2


case $SLURM_ARRAY_TASK_ID in
    1)
        compartment="t_disp_m"
        ;;
    2)
        compartment="s_disp_m"
        ;;
    3)
        compartment="s_disp_t"
        ;;
esac

echo "Processing compartment: $compartment"

### Generate metagenome predictions
       picrust2_pipeline.py -s ./merged_sequence.fasta -i ./${compartment}/${compartment}_pseudo_feature_table_hdf5.biom -o ${compartment}/${compartment}_picrust2_out_pipeline -p $SLURM_NTASKS ;

            add_descriptions.py -i ${compartment}/${compartment}_picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o ${compartment}_picrust2_out_pipeline/pathway_out/pred_metagenome_unstrat_descrip.tsv.gz

            add_descriptions.py -i ${compartment}/${compartment}_picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o ${compartment}_picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

            add_descriptions.py -i ${compartment}/${compartment}_picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o ${compartment}_picrust2_out_pipeline/EC_metagenome_out/path_abun_unstrat_descrip.tsv.gz

echo "finished..."