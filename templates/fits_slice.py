#!/usr/bin/env python

from collections import OrderedDict
from astropy.io import fits
import sys
from json import dump as json_dump
from argparse import ArgumentParser
import numpy as np
import shlex
import time
import os

"""
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=38 --tmp=500G
cd /nvmetmp
cp /astro/mwaeor/dev/sdc3/Image/station_beam_time-varying.fits .
singularity exec -B \$PWD --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    /pawsey/mwa/mwaeor/dev/sdc3-pipeline/templates/fits_slice.py \
    --fits station_beam_time-varying.fits \
    --axis UTC \
    --indices 0
"""

def make_fits_axis_array(hdu, axis):
    count = hdu.header[f"NAXIS{axis}"]
    crval = hdu.header[f"CRVAL{axis}"]
    cdelt = hdu.header[f"CDELT{axis}"]
    crpix = hdu.header[f"CRPIX{axis}"]
    return cdelt * (np.arange(count) + (1 - crpix)) + crval

# THIS IS WRONG
# def make_fits_axis_array(hdu, axis):
#     count = hdu.header[f"NAXIS{axis}"]
#     crval = hdu.header[f"CRVAL{axis}"]
#     cdelt = hdu.header[f"CDELT{axis}"]
#     crpix = hdu.header[f"CRPIX{axis}"]
#     return cdelt * (np.arange(count) - crpix) + crval


def get_parser():
    parser = ArgumentParser(
        description="split a fits file along the given axis")

    parser.add_argument('--fits', help='source fits file')
    parser.add_argument('--hdu', help='name of HDU containing axis to split', default="PRIMARY")
    parser.add_argument('--axis', help='CTYPE of axis to split')
    parser.add_argument('--indices', help='optional: only split these indices',
                        nargs='+', type=int, default=None)
    parser.add_argument('--digits', help='optional: digits used for padding',
                        type=int, default=None)


    return parser

def main():
    parser = get_parser()

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        # is being called directly from nextflow
        args = parser.parse_args(shlex.split("${args}"))

    print("env:")
    for k, v in os.environ.items():
        print(f' -> {k}: {v}')

    print(f"[{time.time()}] args: {args}")

    fits_base, fits_ext = os.path.splitext(args.fits)

    with fits.open(args.fits) as hdus:
        hdu = hdus[args.hdu]
        axis = None
        naxes = hdu.header['NAXIS']
        for a in range(naxes):
            if hdu.header.get(f"CTYPE{a+1}") == args.axis:
                axis = a
                axis_np = naxes - axis - 1
                print(f"found axis {a}(+1) -> np {axis_np}")
                break
        if not axis:
            raise UserWarning(f"Axis not found in header:{chr(0x0a)}{hdu.header}")

        values = np.array(make_fits_axis_array(hdu, axis + 1).tolist())
        header = hdu.header.copy()
        data = hdu.data.copy()
    header["NAXIS"] -= 1
    ax_len = len(values)
    digits = args.digits or len(str(ax_len))
    indices = args.indices or range(ax_len)

    for idx, value, slice_data in zip(indices, values[indices], np.rollaxis(data.take(indices, axis=axis_np), axis_np)):
        idx0 = str(idx).zfill(digits)
        slice_file = f"{fits_base}.{args.axis}{idx0}{fits_ext}"
        print(f"writing slice {idx0} to {slice_file} with value {value}")
        slice_hdus = fits.HDUList(hdus.copy())
        slice_hdu = slice_hdus[args.hdu]
        slice_hdu.data = slice_data
        print(f"shape data ({data.shape}) -> slice ({slice_data.shape})")
        for key in ['NAXIS', 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CROTA']:
            slice_hdu.header.pop(f"{key}{axis+1}", None)
        for a in range(axis+1, header['NAXIS']):
            for key in ['NAXIS', 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CROTA']:
                if f"{key}{a+1}" in slice_hdu.header:
                    # print(f"replacing {key}{a+1} ({slice_hdu.header.get(f'{key}{a+1}')}) with {key}{a} ({slice_hdu.header.get(f'{key}{a}')})")
                    slice_hdu.header[f"{key}{a}"] = slice_hdu.header.pop(f"{key}{a+1}")
        slice_hdu.header["NAXIS"] = header["NAXIS"] - 1

        slice_hdu.header[args.axis] = value
        with open(slice_file, 'wb') as slice_stream:
            slice_hdus.writeto(slice_stream, overwrite=True, output_verify='fix+warn')
            slice_stream.flush()
    with open("done.txt", "w") as f:
        f.write("done")
    print(f"[{time.time()}] os.path.curdir ({os.path.abspath(os.path.curdir)}):")
    for f in os.listdir("."):
        print(f" -> {f} {os.path.getsize(f)}")
    print(os.system(f'sync {fits_base}.*{fits_ext}'))
    print(f'[{time.time()}] sleeping to allow singularity to sync files...')
    time.sleep(5)
    print(f'[{time.time()}] done')


if __name__ == '__main__':
    main()