#!/bin/bash
# zsh -l

export chemeng185_raul_home_dir=/home/flores12
export chemeng185_pedro_home_dir=/home/pedrogm
export chemeng185_sela_home_dir=/home/selaber

export chemeng185_gdrive=/home/flores12/GoogleDrive/TeamDrives/chemeng_185b
export chemeng185_repo=/home/flores12/00_dropbox/01_norskov/00_git_repos/CHEMENG185B_DFT

export dropbox=/home/flores12/00_dropbox

# Exporting COMPENV (ase-espresso needs this set)
export COMPENV='rice'

#| - QE Configuration 
export PATH=/home/flores12/espresso/espresso/bin:$PATH

# Export location of ase-espresso python package ##########
export PYTHONPATH=/home/flores12/00_dropbox/01_norskov/00_git_repos:$PYTHONPATH

# Pseudo Potentials #######################################
export ESP_PSP_PATH=/home/flores12/espresso/pseudo

export SCRATCH=/farmshare/user_data/$USER

# Source openmpi module ###################################
module load openmpi/3.0.0

# /usr/share/lmos/lmos/libexec/lmod
# source /etc/profile.d/lmod.sh
#__|

# Python Configuration
eval "$(/home/flores12/anaconda3/bin/conda shell.bash hook)"
conda activate py27


# Sourcing my configuration files (Aliases and bash function)
source $chemeng185_raul_home_dir/00_dropbox/01_norskov/00_system_settings/00_system_specific/dropbox_paths.sh

source $chemeng185_raul_home_dir/00_dropbox/01_norskov/00_system_settings/aliases.sh
source $chemeng185_raul_home_dir/00_dropbox/01_norskov/00_system_settings/functions.sh


#| - Aliases
# alias chemeng185_repo_sync='bash $repo_185b/scripts/rclone_sync.sh'
alias chemeng185_repo_sync='bash $chemeng185_repo_185b/scripts/rclone_sync.sh'
alias source_bashrc='source $HOME/.bashrc'

alias jobs='python $sc/08_slurm_jobs/jobs.py'

# Custom sbatch command
sbatch_185b()
{
  sbatch --mail-user=$USER@stanford.edu $1
}

#__|


echo "sourced CHEMENG185B bashrc_shared.sh file"
