# Wizards of Oz AU SDC3 Pipeline

## Technical details

Please read the technical details PDF for a description of our approach [Technical Details](SDC3%20reproducibility%20-%20Wizards%20of%20Oz%203D.pdf)

This repository contains only the nextflow pipeline that ties together the rest of the software we used.

## Flows

- beam.nf - make a movie of the time-varying station beam
- transpose.nf
    - transpose the sdc3 data (901 files of 1440 timesteps) into files containing a specified number of channels and timesteps
    - subtract the specified sourcelist
    - grid and estimate the power spectra
    - produce images for diagnostics and source finding.

## Usage

- Ensure a `$profile` has been set up in `nextflow.config` for your cluster or local machine (see: `profiles.dug`, `profiles.garrawarla`).
- Ensure `nextflow` module version > 21 is available, or the `java` module is loaded, and the `nextflow` binary is in your path.
- configure `nextflow.config` to point to the directories containing the input files (e.g. sourcelists)

the basic syntax is

```bash
nextflow run ${flow}.nf -profile $profile -entry $entry
```

The following sections outline the key parameters for each flow, but for a full list of options,
see `nextflow.config`.

### Transposition

different workflows work best with data in different resolutions.

given a glob of uvfs, `--uvfGlob` with a single channel, multiple timestamps:
- split into timeblocks of `--splitTime`
- merge across spectral windows of `--spWidthHz`.

Optional timestep and channel filters can be provided for working on a subset
of the data.

Depending on flags, pass visibilities to the calibration and imaging stage.
- `--doCal` - enable calibration with `--calSrclist`
- `--imgUnsub` - produce images of the unsubtracted data
- `--imgModel` - produce images of the calibration model
- `--imgSub` - produce images of the subtracted data

```bash
nextflow run transpose.nf -entry uvfTranspose # -profle ... -w ...
```

### Imaging

produce images from a glob of visibility files, `--visGlob`

imaging jobs are split into multiple spectral windows by `--imgWidthHz`.
Multi frequency synthesis is done at a resolution determined by `--imgCleanWidthHz`

Optionally source find on the image with `--doAegean`

```bash
nextflow run transpose.nf -entry wscleanEntry # -profle ... -w ...
```

## Walkthrough

### sky modelling

The goal is to produce an initial point source sky model for calibration.

We do this by source finding on multiple images, and cross-matching across
channels to fit a spectral index.

We need to transpose the data into 90 * 1MHz visibility files containing the
central 144 timesteps because:
- the beam has time-varying projection effects which are less pronounced around lst=0
- imaging around lst=0 minimizes w-layers
- imaging multiple channels simultaneously lets us enforce a spectral fit on deconvolution components
- cross-matching complexity is reduced with less channels

we tune our imaging settings for optimial flux and position measurements:
- natural weighting
- minimum uv cutoff of 200l to reduce diffuse contributions
- tukey taper 1000l to smooth the sharp edge caused by uv cutoff
- gaussian taper near theoretical minimum psf width (from 8as at 106MHz to 4as at 196MHz)
- 2as pixel scale: minimum 2 pixels per beam
- maximum cleaning scale 100as
- enforce circular beam to avoid false ellipticity in fitted components
- auto mask and threshold at 5 sigma to avoid cleaning sidelobes

due to memory constraints we are limited to about 9000*9000 pixels, so we can only model 5 degrees
of the sky at a time. We split the sky into 5 degree tiles and image each tile separately.

we start by modelling unsubtracted images of the central field.

```bash
nextflow run transpose.nf -entry uvfTranspose -c skymodel.conf \
    --imgUnsub --doAegean --outModelSuffix "central" \
    # -profle ... -w ...
```

then we subtract the central skymodel and model the residual images in the outer patches.

```bash
nextflow run transpose.nf -entry uvfTranspose -c skymodel.conf --imgSub --doAegean \
    --imgUnsub --doAegean --outModelSuffix "central" \
    # -profle ... -w ...
```

merge the skymodels and iterate through the remaining patches