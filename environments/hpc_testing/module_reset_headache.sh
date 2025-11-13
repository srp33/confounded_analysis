module -t avail  &> before.txt 
module restore
bash /apps/spack/orc_scripts/expose_spack_modules.sh 
module -t avail  &> after.txt 

diff before.txt after.txt > diff.txt