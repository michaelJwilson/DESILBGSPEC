# interactive
#
# addrepo /global/homes/m/mjwilson/DESILBGSPEC/redrock/

export PETAL=4
export VERSION=v4.1 # 4.1 on Jan. 4th 2022; 3.1 on 4th October 2021; 2.1 on June 29 2021.

# [desi-data 5773]
export REDUX=/global/cfs/cdirs/desi/spectro/redux/

export FUJI_TEST3=f3/
export EVEREST=everest/

export RELEASE=$FUJI_TEST3
# export RELEASE=$FUJI_TEST3

export OUTDIR=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$RELEASE/$VERSION/

# Software.
export PATH=/global/homes/m/mjwilson/DESILBGSPEC/redrock/bin/:$PATH
export PYTHONPATH=/global/homes/m/mjwilson/DESILBGSPEC/redrock/py/:$PYTHONPATH

export RR_TEMPLATE_DIR=/global/cscratch1/sd/mjwilson/DESILBGSPEC/templates/

# Clauds: 80871 thru-20210430, 80872 thru-20210512.
#
#
# --  Input and output paths --
#
# Deprecated: export TILE=80871
# Deprecated: export COADD=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$TILE/deep/coadd-$PETAL-$TILE-deep.fits 
# Deprecated: export RRH5=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$TILE/deep/redrock-$PETAL-$TILE-deep.h5
# Deprecated: export RRZ=/global/cscratch1/sd/mjwilson/DESILBGSPEC/$TILE/deep/zbest-$PETAL-$TILE-deep.fits

mask_list=( 'DESILBG_TMG_FINAL'
            'DESILBG_BXU_FINAL'
            'DESILBG_G_FINAL')

# Clauds; cumulative night 20210430
export TILE=80871
export NIGHT=20210430

# Clauds-2; cumulative night 20210512
# export TILE=80872
# export NIGHT=20210512

# mask_list=( 'HETDEX_MAIN'
#             'HETDEX_HP')

# HETDEX; cumulative night 20210408. 
# export TILE=80869                                                                                                                              
# export NIGHT=20210408  
# e.g. /global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/80869/20210408

# HETDEX-2; cumulative night 20210513.
# export TILE=80870
# export NIGHT=20210513 

# HETDEX-3; cumulative night 20210416.
# export TILE=80862
# export NIGHT=20210416

# RRZ:   --outfile  output FITS file with best redshift per target.
# RRH5:  --details  output file for full redrock fit details.

export COADD=$REDUX/$RELEASE/tiles/cumulative/$TILE/$NIGHT/coadd-$PETAL-$TILE-thru$NIGHT.fits

export  SYM=$OUTDIR/tiles/cumulative/$TILE/$NIGHT/coadd-$PETAL-$TILE-thru$NIGHT.fits
export RRH5=$OUTDIR/tiles/cumulative/$TILE/$NIGHT/rrdetails-$PETAL-$TILE-thru$NIGHT.h5
export  RRZ=$OUTDIR/tiles/cumulative/$TILE/$NIGHT/redrock-$PETAL-$TILE-thru$NIGHT.fits

# Posterity:
#   export COADD=/global/cscratch1/sd/mjwilson/DESILBGSPEC/June21/$TILE/coadd-$TILE-$PETAL.fits
#   export RRH5=/global/cscratch1/sd/mjwilson/DESILBGSPEC/June21/$TILE/$VERSION/redrock-$PETAL-$TILE.h5
#   export RRZ=/global/cscratch1/sd/mjwilson/DESILBGSPEC/June21/$TILE/$VERSION/zbest-$PETAL-$TILE.h5

echo
echo 'Input:'
echo
echo $COADD
echo
echo 'Output:'
echo
echo $SYM
echo $RRH5
echo $RRZ
echo
echo
echo 'Linking'
echo
echo 'ln -s '$COADD' '$SYM
echo
echo 'Command:'
echo

for MASK in "${mask_list[@]}"
do
echo 'prospect_pages --spectra_files '$SYM' --zcat_files '$RRZ' --redrock_details_files '$RRH5' --outputdir '$OUTDIR'/tiles/cumulative/'$TILE'/'$NIGHT'/ --mask_type SV1_SCND_TARGET --targeting_mask '${MASK}' --template_dir /global/cscratch1/sd/mjwilson/DESILBGSPEC/templates/ --titlepage_prefix prospect_'${MASK}'_'$PETAL'_'$VERSION 
echo
done

for MASK in "${mask_list[@]}"
do
    prospect_pages --spectra_files $SYM \
                   --zcat_files $RRZ \
                   --redrock_details_files $RRH5 \
                   --outputdir $OUTDIR/tiles/cumulative/$TILE/$NIGHT/ \
                   --mask_type SV1_SCND_TARGET \
                   --targeting_mask ${MASK} \
                   --template_dir /global/cscratch1/sd/mjwilson/DESILBGSPEC/templates/ \
                   --titlepage_prefix prospect_${MASK}_$PETAL'_'$VERSION
done

