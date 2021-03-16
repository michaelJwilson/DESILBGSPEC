export PETAL=9
export RR_TEMPLATE_DIR=/global/cscratch1/sd/mjwilson/DESILBGSPEC/templates/

srun -N 8 -n 256 -c 2 rrdesi_mpi -o /global/cscratch1/sd/mjwilson/DESILBGSPEC/80871/deep/redrock-$PETAL-80871-deep.h5 -z /global/cscratch1/sd/mjwilson/DESILBGSPEC/80871/deep/zbest-$PETAL-80871-deep.fits  /global/cscratch1/sd/mjwilson/DESILBGSPEC/80871/deep/coadd-$PETAL-80871-deep.fits
