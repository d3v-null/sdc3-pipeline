#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fitsJson {
    input:
        path(fits)

    output:
        path(json)

    storeDir "${params.outdir}"
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

    storeDir "${params.outdir}/${meta.store}"

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

process thumbnail {
    input:
        tuple path(img), val(meta)
    output:
        tuple path("${img.baseName}.png"), val(meta)

    storeDir "${params.outdir}/${meta.store}"

    tag "${img.baseName}"

    label "python"
    time 10.min

    script:
    def args = [
        fits: img,
        title: "${meta.time} ${meta.freq/1e6}MHz",
        output_name: "${img.baseName}.png",
        dpi: 50,
    ]
    // argstr = "${params.thumbnail_args}"
    argstr = (args + meta).findAll { k, v ->
            v != null && [
                'fits', 'output_name', 'title', 'dpi', 'limit', 'vmin_quantile',
                'vmax_quantile', 'transparent', 'symmetric', 'cmap'
            ].contains(k)
        }
        .collect { k, v -> ["--${k}"] + (v instanceof List ? v : [v]) }
        .flatten()
        .collect { "\"${it}\"" }
        .join(' ')

    template "/pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/thumbnail.py"
}

process ffmpeg {
    input:
        tuple val(name), path("??????????.png"), val(cachebust)

    output:
        tuple path("${name}.mp4"), path("${dot_cachebust}")

    storeDir "${params.outdir}"
    stageInMode "symlink"

    tag "${name}"

    label "ffmpeg"

    script:
    dot_cachebust = ".${name}.${cachebust}.cachebust"
    """
    #!/bin/bash -eux
    touch ${dot_cachebust}
    ffmpeg -y -framerate 60 \
        -pattern_type glob -i "??????????.png" \
        -vcodec libx264 \
        -vf "scale='min(3840,iw)':'min(2160,ih)':force_original_aspect_ratio=decrease,pad=ceil(iw/2)*2:ceil(ih/2)*2" \
        -pix_fmt yuv420p \
        "${name}.mp4"
    """
}

import groovy.json.JsonSlurper
jslurp = new JsonSlurper()
def parseJson(path) {
    // TODO: fix nasty hack to deal with NaNs
    jslurp.parseText(path.getText().replaceAll(/-?(NaN|Infinity)/, '"$1"'))
}

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

workflow sliceTimes {
    take:
        fitsRoot

    main:
        fitsRoot | fitsJson

        fitsRoot.combine(fitsJson.out).map{ fits, json ->
            def data = parseJson(json)
            def prime_axes = data['PRIMARY']['axes']
            def time_axis = prime_axes.find { it['CTYPE'] == 'UTC' }
            def freq_axis = prime_axes.find { it['CTYPE'] == 'FREQ' }
            def meta = [
                axis: 'UTC',
                store: 'slice-UTC',
                times: time_axis['values'],
                freqs: freq_axis['values'],
            ]
            def ntimes = meta.times.size
            meta.digits = ntimes.toString().length()
            meta.indices = 0..<ntimes
            // meta.indices = [0]
            [fits, meta]
        } | fitsSlice

    emit:
        sliceArgs = fitsSlice.out.flatMap { meta, slices ->
                [meta.times, slices instanceof List ? slices : [slices]].transpose().collect { time, slice ->
                    def nfreqs = meta.freqs.size
                    def newMeta = [
                        axis: 'FREQ',
                        store: "${meta.store}-FREQ",
                        time: time,
                        digits: nfreqs.toString().length(),
                        // indices: 0..<nfreqs
                        indices: [nfreqs-1]
                    ]
                    [slice, deepcopy(meta) + newMeta]
                }
            }
}

workflow sliceFreqs {
    take:
        // tuple of (slice, meta)
        sliceArgs

    main:
        sliceArgs | fitsSlice

    emit:
        sliceArgs = fitsSlice.out.flatMap { meta, slices ->
            [meta.freqs, slices instanceof List ? slices : [slices]].transpose().collect { freq, slice ->
                def newMeta = [freq: freq]
                [slice, deepcopy(meta) + newMeta]
            }
        }
}

workflow {
    print("processing ${params.time_varying_station_beam}")
    // def beam = file(params.time_varying_station_beam)

    channel.of(file(params.time_varying_station_beam))
        | sliceTimes
        | sliceFreqs
        | thumbnail

    thumbnail.out.map { img, meta -> [meta.freq, img] }
        .groupTuple(by: 0)
        .map { freq, files ->
            def name = "freq_${freq.round()}"
            def latest = files.collect { file -> file.lastModified() }.max()
            def cachebust = "${latest}_x" + sprintf("%04d", files.size())
            def sorted = files.collect { path -> file(deepcopy(path.toString())) }.sort()
            [name, sorted, cachebust]
        }
        | ffmpeg
}