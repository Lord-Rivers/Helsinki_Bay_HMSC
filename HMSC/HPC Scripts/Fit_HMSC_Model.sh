#!/bin/bash
#SBATCH --job-name=HMSC_model_fit_Tensorflow
#SBATCH  -M ukko
#SBATCH --partition=gpu 
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus-per-node=1
#SBATCH --constraint=v100 #Restrict this to only running on V100 GPU's because the A100's decrease convergence performence
#SBATCH --cpus-per-gpu=2
#SBATCH --time=30:00:00
#SBATCH --verbose
#SBATCH --array=0-3
#SBATCH --output=R-%x.%A_%a.out

module purge
module load Python/3.12.3-GCCcore-13.3.0
source $HOME/myenv/bin/activate

model_name=$NAME

Area=$AREA
samples=$SAMP
thin=$THIN
chains=4
divide=$samples*$thin; by=2; (( transient=(divide+by-1)/by ))
divide=3*$samples*$thin; by=2*100; ((verbose=(divide/by)))
#verbose=100
data_path=$(printf "%s/%s" $Area $model_name)
#screen_text=$(printf "#######\n Model Run starting\n Area: %s model_name: %s\n Samples: %.4d\n Thining: %.2d Transient steps: %.d\n verbose %.d\n######\n" $Area $model_name $samples $thin $transient $verbose)

input_path=$data_path/$(printf "INIT/HPC_INIT_samples_%.4d_thin_%.2d_chains_%.1d.rds" $samples $thin $chains)
output_path=$data_path/$(printf "Sampled/HPC_samples_%.4d_thin_%.2d_chain_%.1d.rds" $samples $thin $SLURM_ARRAY_TASK_ID)

printf "#####################\nModel Run starting\nArea: %s model_name: %s\nSamples: %.4d\nThining: %.2d Transient steps: %.d\nVerbose %.d\n#####################\n" $Area $model_name $samples $thin $transient $verbose
echo $data_path
echo $input_path
echo $output_path

srun python Gibs_sampling_tensor_flow.py --samples $samples --transient $transient --thin $thin --verbose $verbose --hmcthin 0 --input $input_path --output $output_path --chain $SLURM_ARRAY_TASK_ID