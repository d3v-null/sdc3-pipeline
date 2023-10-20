#!/usr/bin/env python

from collections import OrderedDict
from astropy.io import fits
import sys
from json import dump as json_dump
from argparse import ArgumentParser
import numpy as np

"""
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=38 --tmp=500G
cd /nvmetmp
cp /astro/mwaeor/dev/sdc3/Image/station_beam_time-varying.fits .
singularity exec -B \$PWD --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    /pawsey/mwa/mwaeor/dev/sdc3-pipeline/templates/fits_json.py \
    --fits station_beam_time-varying.fits \
    --json /astro/mwaeor/dev/sdc3/store/station_beam_time-varying.json
"""


def make_fits_axis_array(hdu, axis):
    count = hdu.header[f"NAXIS{axis}"]
    crval = hdu.header[f"CRVAL{axis}"]
    cdelt = hdu.header[f"CDELT{axis}"]
    crpix = hdu.header[f"CRPIX{axis}"]
    return cdelt * (np.arange(count) + (1 - crpix)) + crval


def get_parser():
    parser = ArgumentParser(
        description="get basic metadata from fits file in json format")

    parser.add_argument('--fits', help='source fits file')
    parser.add_argument('--json', help='Name of output json file')

    return parser

def main():
    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args([
            "--fits=${fits}",
            "--json=${json}",
        ]  # + shlex.split("${args}")
        )

    data = OrderedDict()
    with fits.open(args.fits) as hdus:
        for hdu in hdus:
            name = hdu.name
            data[name] = OrderedDict([('axes', [])])
            for axis in range(hdu.header['NAXIS']):
                axis_data = OrderedDict()
                for key in ['NAXIS', 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CROTA', 'CUNIT']:
                    if f"{key}{axis+1}" in hdu.header:
                        axis_data[key] = hdu.header[f"{key}{axis+1}"]
                if all([key in axis_data for key in ['NAXIS', 'CRVAL', 'CDELT', 'CRPIX']]):
                    axis_data['values'] = make_fits_axis_array(hdu, axis+1).tolist()
                data[name]['axes'].append(axis_data)
            for key in hdu.header:
                data[name][key] = hdu.header[key]

    with open(args.json, 'w') as json:
        json_dump(data, json, indent=4)


if __name__ == '__main__':
    main()