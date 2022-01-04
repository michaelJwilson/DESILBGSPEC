export PATH=/global/homes/m/mjwilson/DESILBGSPEC/redrock/bin/:$PATH
export PYTHONPATH=/global/homes/m/mjwilson/DESILBGSPEC/redrock/py/:$PYTHONPATH

mask_list=( 'HETDEX_MAIN'
            'HETDEX_HP')

# /global/cfs/cdirs/desi/spectro/redux

for MASK in "${mask_list[@]}"
do
    prospect_pages --datadir ${DESI_SPECTRO_REDUX}'/everest/tiles/cumulative' \
                   --dirtree_type cumulative \
                   --tiles 80870 \
                   --mask_type 'SV1_SCND_TARGET' --targeting_mask ${MASK} \
                   --outputdir /global/cscratch1/sd/mjwilson/DESILBGSPEC/80870/cumulative/20210513/v3.1/prospect/ \
                   --with_zcatalog --with_multiple_models \
                   --top_metadata TARGETID COADD_NUMEXP COADD_EXPTIME FIBER --titlepage_prefix ${MASK}'_everest'
done
