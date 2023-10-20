#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH -c 32
#SBATCH --clusters=garrawarla
#SBATCH --partition=gpuq
#SBATCH --account=mwaeor
#SBATCH --mem=200G
#SBATCH --tmp=200G
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --array=0-1

# salloc --nodes=1 --mem=200G --gres=gpu:1 --time=2:00:00 --clusters=garrawarla --partition=gpuq --account=mwaeor --tasks 1 --cpus-per-task=20 --tmp=200G

export srclist="/astro/mwaeor/dev/sdc3/catalog/sdc3_inner_lobes.fits" # <- YOUR SOURCELIST HERE
export visName="sub_lobes_points"                                      # <- NAME YOUR OUTPUT
export hypArgs="/astro/mwaeor/cjordan/sdc3/sdc3.toml --timesteps 0"   # <- single timestep, saves time.
# export hypArgs="/astro/mwaeor/cjordan/sdc3/sdc3.toml"
# this picks which range of channels to image out of the 901 channels.
export lo=$((SLURM_ARRAY_TASK_ID * 150))
export hi=$(((SLURM_ARRAY_TASK_ID+1) * 150))

module use /astro/mwaeor/software/modulefiles
module load hyperdrive/sdc3 wsclean/2.9 singularity
cd /nvmetmp; mkdir deleteme; cd deleteme
export sdc3store=/astro/mwaeor/dev/sdc3/store
export mergeName="106-196MHz_ts0720-0731"
# export mergeName="106-196MHz_ts0060-0071"
cp ${sdc3store}/uvfMerge/${mergeName}.uvfits .                        # <- unsub
# cp ${sdc3store}/mergeSub/${mergeName}/${mergeName}_sub.uvfits .       # <- sub lobes

# you can either subtract from a uvfits or convert it to ms

# subtract
export vis="106-196MHz_ts0060_${visName}.ms"
hyperdrive vis-sub $hypArgs \
    --data ${mergeName}.uvfits \
    --source-list $srclist \
    --outputs $vis

# convert
# export vis="${mergeName}_${visName}.ms"
# hyperdrive vis-convert $hypArgs \
#     --data "${mergeName}_sub.uvfits" \
#     --outputs $vis

# convert by doing a subtraction
# export vis="${mergeName}_${visName}.ms"
# hyperdrive vis-sub $hypArgs \
#     --data ${mergeName}_sub.uvfits \
#     --source-list $srclist \
#     --num-sources 0 \
#     --outputs $vis

# convert with casa
# export vis="${mergeName}_sub.ms"
# singularity exec \
#     --bind ${PWD}:/tmp --writable-tmpfs --pwd /tmp --home ${PWD} --cleanenv \
#     /pawsey/mwa/singularity/casa5/casa5.sif \
#     casa -c "importuvfits('${mergeName}_sub.uvfits', '${mergeName}_sub.ms');"

# this tricks wsclean into giving us its cleaning ssourcelist
singularity exec \
    --bind ${PWD}:/tmp --writable-tmpfs --pwd /tmp --home ${PWD} --cleanenv \
    /pawsey/mwa/singularity/casa5/casa5.sif \
    casa -c "tb.open('${vis}/POLARIZATION', nomodify=False); tb.putcell(rownr=0,columnname='CORR_TYPE',thevalue=[1]); tb.close(); print('done')"

# export lo=0
# export hi=150
export imgName="wsclean_${vis%%.ms}_ch${lo}-${hi}"
wsclean \
  -name $imgName \
  -reorder \
  -use-wgridder \
  -parallel-gridding 10 \
  -oversampling 4095 \
  -kernel-size 15 \
  -nwlayers 1000 \
  -grid-mode kb \
  -taper-edge 100 \
  -padding 2 \
  -size 8192 8192 \
  -scale 8asec \
  -weight uniform \
  -multiscale \
  -super-weight 4 \
  -auto-threshold 4 \
  -mgain 0.8 \
  -pol I \
  -channel-range $lo $hi \
  -channels-out 15 \
  -join-channels \
  -niter 1000 \
  -save-source-list \
  -minuv-l 400 \
  -fit-spectral-log-pol 2 \
  -multiscale-scales 1,2,3,4,5,6,7,8,9,10 \
  -no-update-model-required \
  ${vis}
  # -taper-gaussian 60 \


mkdir -p /astro/mwaeor/${USER}/sdc3/img-${visName}

(set -x; cp ${imgName}*-MFS-{psf,image,model,residual}.fits ${imgName}*.txt /astro/mwaeor/${USER}/sdc3/img-${visName})

printf '\a'

# zip -r ${vis}.zip ${vis}
# (set -x; cp ${vis}.zip /astro/mwaeor/${USER}/sdc3/img-${visName})

# echo ${imgName}*-MFS-dirty.fits
# echo /astro/mwaeor/${USER}/sdc3/img-${visName}

## beam correction

# salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=38 --tmp=500G
# module load singularity
# cd /nvmetmp

# export LST=0.0 # <- your LST here. it's easiest to just image timestep720 which is LST=0.0 (i think)
# # export FITS=/astro/mwaeor/dev/sdc3/img-sub_lobes/wsclean_106-196MHz_ts0720_sub_lobes_ch150-300-MFS-image.fits
# for FITS in \
#   "/astro/mwaeor/dev/sdc3/img-sub_lobes_trecs/wsclean_106-196MHz_ts0720_sub_lobes_trecs_ch0-150-MFS-image.fits" \
#   "/astro/mwaeor/dev/sdc3/img-sub_lobes_trecs/wsclean_106-196MHz_ts0720_sub_lobes_trecs_ch0-150-MFS-model.fits" \
#   "/astro/mwaeor/dev/sdc3/img-sub_lobes_trecs/wsclean_106-196MHz_ts0720_sub_lobes_trecs_ch150-300-MFS-image.fits" \
#   "/astro/mwaeor/dev/sdc3/img-sub_lobes_trecs/wsclean_106-196MHz_ts0720_sub_lobes_trecs_ch150-300-MFS-model.fits" \
# ; do
#   singularity exec -B $PWD --cleanenv --home /astro/mwaeor/$USER/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
#       /pawsey/mwa/mwaeor/dev/sdc3-pipeline/templates/beam_correct.py \
#       --lst 0.0 \
#       --fits $FITS
# done

