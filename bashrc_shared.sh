#!/bin/bash

export raul_home_dir=/home/flores12
export pedro_home_dir=/home/pedrogm
export sela_home_dir=/home/selaber

export gdrive_185b=/home/flores12/GoogleDrive/TeamDrives/chemeng_185b
export repo_185b=/home/flores12/00_dropbox/01_norskov/00_git_repos/CHEMENG185B_DFT

#| - QE Configuration 
# export PATH=/home/flores12/usr/q-e/bin:$PATH
# export PATH=/afs/ir/users/f/l/flores12/quantum_espresso/00/espresso/bin:$PATH
# export PATH=/afs/ir/users/f/l/flores12/quantum_espresso/00/espresso/bin:$PATH
export PATH=/home/flores12/espresso/espresso/bin:$PATH

# Export location of ase-espresso python package
export PYTHONPATH=/home/flores12/00_dropbox/01_norskov/00_git_repos/ase-espresso/espresso:$PYTHONPATH

# export SCRATCH=/farmshare/user_data/flores12
export SCRATCH=/farmshare/user_data/$USER

#__|

#| - Aliases
alias repo_185b_sync='bash $repo_185b/scripts/rclone_sync.sh'
alias source_bashrc='source $HOME/.bashrc'

# Custom sbatch command
sbatch_185b()
{
  sbatch --mail-user=$USER@stanford.edu $1
}

#__|



export ESP_PSP_PATH=/home/flores12/espresso/pseudo


# Source openmpi module
# /usr/share/lmos/lmos/libexec/lmod
# source /etc/profile.d/lmod.sh
module load openmpi/3.0.0




echo "sourced CHEMENG185B bashrc_shared.sh file"
