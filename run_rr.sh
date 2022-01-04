# addrepo /global/homes/m/mjwilson/DESILBGSPEC/redrock

export PATH=/global/homes/m/mjwilson/DESILBGSPEC/redrock/bin/:$PATH
export PYTHONPATH=/global/homes/m/mjwilson/DESILBGSPEC/redrock/py/:$PYTHONPATH

export PETAL=8
export VERSION=v3.1 # 3.1 on 4th October 2021; 2.1 on June 29 2021.
export RR_TEMPLATE_DIR=/global/cscratch1/sd/mjwilson/DESILBGSPEC/templates/ 

# export TILE=80871
# export COADD=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$TILE/deep/coadd-$PETAL-$TILE-deep.fits 
# export RRH5=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$TILE/deep/redrock-$PETAL-$TILE-deep.h5
# export RRZ=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$TILE/deep/zbest-$PETAL-$TILE-deep.fits

# /global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/80869/20210408
# export TILE=80869 # HETDEX; cumulative night 20210408.
# export NIGHT=20210408

export TILE=80870   # HETDEX-2; cumulative night 20210513.
export NIGHT=20210513 

# export TILE=80862 # HETDEX-3.
# export TILE=80871 # COSMOS.
# export TILE=80872 # COSMOS-2

export COADD=/global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/$TILE/$NIGHT/coadd-$PETAL-$TILE-thru$NIGHT.fits
export RRH5=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$TILE/cumulative/$NIGHT/$VERSION/rrdetails-$PETAL-$TILE-thru$NIGHT.h5
export RRZ=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$TILE/cumulative/$NIGHT/$VERSION/redrock-$PETAL-$TILE-thru$NIGHT.fits

# export COADD=/global/cscratch1/sd/mjwilson/DESILBGSPEC/June21/$TILE/coadd-$TILE-$PETAL.fits
# export RRH5=/global/cscratch1/sd/mjwilson/DESILBGSPEC/June21/$TILE/$VERSION/redrock-$PETAL-$TILE.h5
# export RRZ=/global/cscratch1/sd/mjwilson/DESILBGSPEC/June21/$TILE/$VERSION/zbest-$PETAL-$TILE.h5

echo $COADD
echo $RRH5
echo $RRZ

# srun -N 16 -n 512 -c 2 rrdesi_mpi -i $COADD -o $RRH5 -d $RRZ
