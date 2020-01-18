# TEMP

rc_185b_2()
#| - CHEMENG185B course files
{
echo ">>> Syncing CHEMENG185B folder"

# /home/raulf2012/Dropbox/01_norskov/00_git_repos/CHEMENG185B_DFT

export root_tmp=01_norskov/00_git_repos/CHEMENG185B_DFT

rclone --transfers=8 sync \
raul_dropbox:$root_tmp \
$dropbox/$root_tmp &
}
