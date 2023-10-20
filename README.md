# Wizards of Oz AU SDC3 Pipeline

## Technical details

todo: link the google doc

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

### uv transpose

initial transpose implementation used casa measurement sets, this proved to be completely inadequate for the data volume, and better results were achieved using uvfits instead.

```bash
nextflow run transpose.nf -entry uvfTranspose # -profle ... -w ...
```