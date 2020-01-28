export root_tmp=01_norskov

rclone sync \
raul_dropbox:$root_tmp/00_git_repos/ase-espresso \
$dropbox/$root_tmp/00_git_repos/ase-espresso \
--tpslimit 10 \
--transfers=4 \
--exclude .git/ \
--exclude __old__/ \
--exclude __temp__/ \
--exclude out_data/ \
--exclude out_plot/ \
--exclude pymatgen/ \
--exclude pymatgen_mine/ \
--verbose

