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
        compartment="mucus"
        ;;
    2)
        compartment="tissue"
        ;;
    3)
        compartment="skeleton"
        ;;
esac

echo "Processing compartment: $compartment"

### Generate metagenome predictions
       picrust2_pipeline.py -s ./merged_sequence.fasta -i ./${compartment}/${compartment}_psudo_feature_table_hdf5.biom -o ${compartment}_picrust2_out_pipeline -p $SLURM_NTASKS


done

echo "finished... at $date..."
