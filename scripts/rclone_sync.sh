echo ">>> Syncing CHEMENG185B folder"

export root_tmp=01_norskov/00_git_repos/CHEMENG185B_DFT

rclone --transfers=8 sync \
raul_dropbox:$root_tmp \
$dropbox/$root_tmp
