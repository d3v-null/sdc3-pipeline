#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// You can run this in garrawarla with
// module load java/17
// /astro/mwaeor/software/nextflow/nextflow-22.10.8-all run sourcefind.nf -w /astro/mwaeor/$USER/sdc3/work -profile garrawarla

include { bane } from './transpose.nf'
include { sourceFind } from './transpose.nf'

import groovy.transform.Synchronized
@Synchronized
def deepcopy(orig) {
    bos = new ByteArrayOutputStream()
    oos = new ObjectOutputStream(bos)
    oos.writeObject(orig); oos.flush()
    bin = new ByteArrayInputStream(bos.toByteArray())
    ois = new ObjectInputStream(bin)
    return ois.readObject()
}

import groovy.json.JsonSlurper
jslurp = new JsonSlurper()
def parseJson(path) {
    // TODO: fix nasty hack to deal with NaNs
    jslurp.parseText(path.getText().replaceAll(/-?(NaN|Infinity)/, '"$1"'))
}

process fitsJson {
    input:
        path(fits)

    output:
        path(json)

    storeDir "${params.storeDir}"
    tag "${fits.baseName}"
    label "python"
    // label "nvme"

    script:
    json = "${fits.baseName}.json"
    template "fits_json.py"
}

process fitsSlice {
    input:
        tuple path(fits), val(meta)

    output:
        tuple val(meta), path(out_glob)

    storeDir "${params.storeDir}/${meta.store}"

    tag "${fits.baseName}.${meta.axis}{${indices[0]}..${indices[-1]}}"

    label "python"
    label "nvme"
    label "mem_half"

    time { 1.hour }

    errorStrategy 'terminate'

    script:
    indices = meta.indices.collect { it.toString().padLeft(meta.digits, '0') }
    slice_glob = meta.axis + ( indices.size > 1 ? "{${indices.join(',')}}" : indices[0] )
    out_glob = "${fits.baseName}.${slice_glob}.${fits.extension}"
    args = (meta + [fits: fits]).findAll { k, v -> v != null && ['fits', 'hdu', 'axis', 'indices', 'digits'].contains(k) }
        .collect { k, v -> ["--${k}"] + (v instanceof List ? v : [v]) }
        .flatten()
        .join(' ')

    template "fits_slice.py"
}

process beamCorrect {
    input:
        tuple val(freq), val(meta), path(fits), path(beam)
    output:
        tuple val(freq), val(meta), path(corrected)

    storeDir "${params.storeDir}/${meta.store?:'mergeImg'}"

    tag "${fits.baseName}"

    label "python"

    time {5.min * task.attempt}

    script:
    corrected = "${fits.baseName}._pbeam.fits"
    args = "--fits ${fits} --beam ${beam} --threshold ${params.beamThreshold}"
    template "beam_correct.py"
}

workflow sliceUni {
    take:
        uni
    main:
        // slice the uniform image cube in frequency
        uni | fitsJson
        uni.combine(fitsJson.out).map { fits, json ->
                def data = parseJson(json)
                def prime_axes = data['PRIMARY']['axes']
                def freq_axis = prime_axes.find { it['CTYPE'] == 'FREQ' }
                def meta = [
                    axis: 'FREQ',
                    store: 'uni-slice-FREQ',
                    freqs: deepcopy(freq_axis['values']),
                ]
                def nfreqs = meta.freqs.size
                meta.digits = nfreqs.toString().length()
                meta.indices = 0..<nfreqs
                [fits, deepcopy(meta)]
            }
            | fitsSlice
    emit:
        fitsSlice.out.flatMap { meta, slices ->
            [meta.freqs, slices instanceof List ? slices : [slices]].transpose().collect { freq, slice ->
                def newMeta = [freq: freq]
                newMeta.centerFreqHz = meta.freq
                [freq, slice, deepcopy(meta) + newMeta]
            }
        }
}

workflow sliceBeam {
    take:
        beam
    main:
        // slice the beam image cube in frequency
        beam | fitsJson
        beam.combine(fitsJson.out).map { fits, json ->
                def data = parseJson(json)
                def prime_axes = data['PRIMARY']['axes']
                def freq_axis = prime_axes.find { it['CTYPE'] == 'FREQ' }
                def meta = [
                    axis: 'FREQ',
                    store: 'beam-slice-FREQ',
                    freqs: deepcopy(freq_axis['values']),
                ]
                def nfreqs = meta.freqs.size
                meta.digits = nfreqs.toString().length()
                meta.indices = 0..<nfreqs
                [fits, deepcopy(meta)]
            }
            | fitsSlice
    emit:
        fitsSlice.out.flatMap { meta, slices ->
            [meta.freqs, slices instanceof List ? slices : [slices]].transpose().collect { freq, slice ->
                def newMeta = [freq: freq]
                newMeta.centerFreqHz = meta.freq
                [freq, slice, deepcopy(meta) + newMeta]
            }
        }
}

workflow {
    channel.of(file(params.uniform_img_cube)) | sliceUni
    channel.of(file(params.station_beam)) | sliceBeam
    sliceUni.out.join(sliceBeam.out).map { freq, uniSlice, uniMeta , beamSlice, beamMeta ->
            // print("freq: ${freq} ${uniMeta} ${beamMeta}")
            uniMeta.centerFreqHz = freq
            [freq, uniMeta, uniSlice, beamSlice]
        }
        | beamCorrect
    beamCorrect.out.map { freq, meta, corrected ->
            ['', meta, corrected]
        }
        | bane


    beamCorrect.out.join(
            bane.out.map { _, meta, bkg, rms -> [meta.freq, bkg, rms]}
        )
        .map { imgName, meta, img, bkg, rms ->
            ['', meta, img, bkg, rms]
        }
        | sourceFind
}