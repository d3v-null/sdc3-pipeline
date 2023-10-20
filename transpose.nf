#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// You can run this in garrawarla with
// module load java/17
// /astro/mwaeor/software/nextflow/nextflow-22.10.8-all run transpose.nf -w /astro/mwaeor/$USER/sdc3/work -profile garrawarla -entry uvfTranspose

if (!params.doCal) {
    params.calSuffix = ""
}

import java.text.SimpleDateFormat
import java.util.LinkedHashMap
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

@Synchronized
def coerceList(x) {
    if (x instanceof List) {
        x
    } else {
        [x]
    }
}

// process cat {
//     input:
//     path(txt)

//     output:
//     path(out)

//     storeDir "s3://devtest2/."

//     shell:
//     out = 'out.txt'
//     """
//     #!/bin/bash -eux
//     cat ${txt} > ${out}
//     ls -al
//     """
// }

// workflow {
//     channel.from(file('s3://devtest/test.txt'))
//         | cat
//         | view { it -> it.text }
// }


import groovy.json.JsonSlurper
import groovy.json.StringEscapeUtils
import java.math.RoundingMode

jslurp = new JsonSlurper()
def parseJson(path) {
    // TODO: fix nasty hack to deal with NaNs
    jslurp.parseText(path.getText().replaceAll(/(NaN|-?Infinity)/, '"$1"'))
}

process uvMeta {
    input:
    tuple val(baseName), path(vis)

    output:
    tuple val(baseName), path(uvmeta)

    storeDir "${params.storeDir}/uvMeta"

    tag "${baseName}"

    label "python"
    time 15.minute

    script:
    uvmeta = "${baseName}.json".toString()
    template "uvmeta.py"
}

process msMeta {
    input:
    tuple val(baseName), path(ms)

    output:
    tuple val(baseName), path(meta)

    storeDir "${params.storeDir}/msMeta"

    tag "${baseName}"

    label "casa"
    time 15.minute

    script:
    meta = "${baseName}.json".toString()
    """
    #!/usr/bin/env casa --nogui --nologger --log2term -c
    import json
    import numpy as np
    tb.open('${ms}')
    times = tb.query('', columns='distinct(TIME)').getcol('TIME')
    tb.close()
    ms.open('${ms}')
    spwinfo = ms.getspectralwindowinfo()['0']
    print('spwinfo0', spwinfo)
    start_freq = spwinfo['Chan1Freq']
    chan_width = spwinfo['ChanWidth']
    num_chan = spwinfo['NumChan']
    ms.close()
    with open('${meta}', 'w') as f:
        json.dump({
            'times': [{'mjd_sec': time} for time in times.tolist()],
            'freqs': np.linspace(start_freq, start_freq + chan_width * num_chan, num_chan).tolist()
        }, f)
    """
}

process uvSplit {
    input:
    tuple val(baseName), val(splitName), path(uvfits), val(times)

    output:
    tuple val(baseName), val(splitName), path(split)

    storeDir "${params.storeDir}/uvSplit"

    tag "${baseName}_${splitName}"

    label "python"
    label "nvme"
    time 15.minute

    script:
    args = "--times ${times.join(' ')}"
    split = "${baseName}_${splitName}.uvfits"
    template "uvsplit.py"
}

process msSplit {
    input:
    tuple val(baseName), val(splitName), path(ms), val(times)

    output:
    tuple val(baseName), val(splitName), path(split)

    storeDir "${params.storeDir}/msSplit"

    tag "${baseName}_${splitName}"

    label "casa"
    time 15.minute

    script:
    split = "${baseName}_${splitName}.ms"
    """
    #!/usr/bin/env casa --nogui --nologger --log2term -c
    tb.open('${ms}')
    tb.query('', columns='distinct(TIME)', sortlist='TIME', name='alltimes.tab')\
        .query('rownr() in [${times.join(",")}]', name='seltimes.tab', style='python')
    tb.query('TIME in [select TIME from seltimes.tab]', style='python')\
        .copy('${split}', deep=True, valuecopy=True, norows=False)
    """

//     """
//     #!/bin/bash -eux

    // times = np.sort(tb.query('', columns='distinct(TIME)').getcol('TIME'))
// cat > split.py << EoF
// tb.open('${ms}')
// tb.taql('select distinct(TIME) from ${ms}').copy('alltimes.tab')
// tb.open('alltimes.tab')
// tb.selectrows([${times.join(',')}]).copy('seltimes.tab')
// tb.open('${ms}')
// tb.taql('select * from ${ms} where TIME in [select TIME from seltimes.tab]').copy('${split}')
// EoF
//     # singularity exec \
//     #     --bind \$PWD:/tmp --writable-tmpfs --pwd /tmp --home \$PWD --cleanenv \
//     #     docker://d3vnull0/casa:latest \
//     casa -c split.py
//     """
}

process uvfMerge {
    input:
    tuple val(mergeName), val(meta), path(uvfs)

    output:
    tuple val(mergeName), val(meta), path(merge)

    storeDir "${params.storeDir}/uvfMerge"

    tag "${mergeName}"
    label "hyperdrive"

    stageInMode "symlink"

    label "cpu_half"
    label "mem_half"
    label "nvme"

    time 1.hour

    script:
    uvfs = coerceList(uvfs)
    merge = "${mergeName}.uvfits"
    """
    #!/bin/bash -eux
    RUST_BACKTRACE=full sdc3_vis_convert -vvv ${uvfs.join(' ')} --time ${meta.ts.join(' ')} --output ${merge}
    """
}
// process uvMerge {
//     input:
//     tuple val(mergeName), path(uvfits)

//     output:
//     tuple val(mergeName), path(merge)

//     storeDir "${params.storeDir}/uvMerge"

//     tag "${mergeName}"

//     module "hyperdrive/dev"

//     label "cpu_half"
//     label "mem_half"
//     label "nvme"

//     time 1.hour

//     script:
//     uvfits = coerceList(uvfits)
//     merge = "${mergeName}.uvfits"
//     """
//     #!/bin/bash -eux
//     RUST_BACKTRACE=full sdc3_vis_convert -vvv ${uvfits.join(' ')} --output ${merge}
//     """
// }

process msMerge {
    input:
    tuple val(mergeName), path(ms)

    output:
    tuple val(mergeName), path(merge)

    storeDir "${params.storeDir}/msMerge"

    tag "${mergeName}"

    label "nvme"
    label "mem_half"
    label "hyperdrive"

    time 4.hour

    script:
    merge = "merge_${mergeName}.ms"
    """
    #!/bin/bash -eux
    sdc3_vis_convert ${ms} --output ${merge}
    """
    // label "casa"
    // time 6.hour

    // script:
    // merge = "merge_${mergeName}.ms"
    // vis = "[" + ms.collect { "'${it}'" }.sort().join(",") + "]"
    // """
    // #!/usr/bin/env casa --nogui --nologger --log2term -c
    // # get the reference frequency and channel width from the first MS
    // vis=${vis}
    // for v in vis:
    //     tb.open(v)
    //     print(v, tb.getcell('DATA', 0))
    //     tb.close()
    // ms.open(vis[0])
    // spwinfo = ms.getspectralwindowinfo()['0']
    // start_freq = qa.tos(qa.quantity(spwinfo['Chan1Freq'], 'Hz'))
    // chan_width = qa.tos(qa.quantity(spwinfo['ChanWidth'], 'Hz'))
    // frame = spwinfo['Frame']
    // num_chan = len(vis)
    // print('spwinfo0', spwinfo)
    // ms.close()
    // virtualconcat(vis=vis, concatvis='v${merge}')
    // mstransform(vis='v${merge}', outputvis='${merge}',
    //     regridms=True, mode='frequency', combinespws=True,
    //     start=start_freq, width=chan_width, nchan=num_chan,
    //     interpolation='linear', outframe=frame, datacolumn='all')
    // tb.open('${merge}')
    // print('merged', tb.getcell('DATA', 0))
    // tb.close()
    // ms.open('${merge}')
    // print('merged spwinfo', ms.getspectralwindowinfo())
    // ms.close()
    // """
}

process hypCal {
    input:
    tuple val(mergeName), val(meta), path(data)

    output:
    tuple val(mergeName), val(meta), path(soln), path(model)

    storeDir "${params.storeDir}/mergeCal/${mergeName}${params.calSuffix}"

    tag "${mergeName}"

    label "mem_half"
    label "cpu_half"
    label "gpu_nvme"
    label "hyperdrive"

    time 45.min

    script:
    model = "model_${mergeName}.ms"
    soln = "soln_${mergeName}.fits"
    // beam_args="--beam ska_airy"
    // array_pos_args="--array-position 116.764 -26.825 0"
    // --source-list ${params.sourcelist}
    """
    hyperdrive di-cal -vvv \
        "${params.hypToml}" \
        ${params.hyp_dical_args} \
        --model-filenames "${model}" \
        --data "${data}" \
        """ + (
            params.subSrclist ? "--source-list ${params.subSrclist}" : ""
        ) + """ \
        --outputs "${soln}"
    """
}

process plotSols {
    input:
    tuple val(mergeName), val(meta), path(soln)
    output:
    tuple val(mergeName), val(meta), path(plotGlob)

    storeDir "${params.storeDir}/mergeCal/${mergeName}"

    tag "${mergeName}"

    label "hyperdrive"

    script:
    plotGlob = "${soln.baseName}_{phases,amps}*.png"
    """
    hyperdrive solutions-plot ${params.hyp_sols_plot_args} ${soln} \
        --max-amp 2
    """
}

// quick image for visibilities
workflow wscleanEntry {
    viss = channel.fromPath("/astro/mwaeor/cjordan/sdc3/106-196MHz_ts0720-0731_sub.uvfits", type: 'any').map { [it.baseName, it] }.take(1)
    viss.flatMap { baseName, vis ->
            def (_, tbName, subName) = baseName.split('_')
            ((params.startFreqHz)..(params.endFreqHz)).by(params.inWidthHz).withIndex()
                .collate((params.imgWidthHz / params.inWidthHz).intValue(), false)
                .findAll { it[-1][-1] != it[0][-1] }
                // .findAll { it[-1][0] > 195e6 }
                // .take(1)
                .collect { fi_chunk ->
                    startFreqHz = fi_chunk[0][0]
                    startFreqIdx = fi_chunk[0][-1]
                    endFreqHz = fi_chunk[-1][0]
                    endFreqIdx = fi_chunk[-1][-1]
                    spwName = getSpwName([startFreqHz, endFreqHz])
                    wscArgs = "${params.wscleanCommon} ${params.wscleanMsw} ${params.wscleanMswExtra}"
                    wscArgs += " -channel-range ${startFreqIdx} ${endFreqIdx}"
                    channels = Math.ceil(params.imgWidthHz / params.imgCleanWidthHz).intValue()
                    if (channels > 1) {
                        wscArgs += " -channels-out ${channels} -join-channels"
                    }
                    meta = [spwName:spwName, startFreqIdx:startFreqIdx, endFreqIdx:endFreqIdx]
                    ["${spwName}_${tbName}_${subName}_msw", meta, vis, wscArgs]
                }
        }
        .take(1)
        | wsclean
        | view()

    wsclean.out.flatMap { imgName, files ->
            def (spwName, tbName, subName, imgType) = imgName.split('_')
            files.findAll { it.baseName =~ /.*MFS-image/ }.collect { fits ->
                ["${spwName}_${subName}_${imgType}", fits]
            }
        }
}

process hypSub {
    input:
    tuple val(mergeName), val(meta), path(data)
    output:
    tuple val(mergeName), val(meta), path(subVis)

    storeDir "${params.storeDir}/mergeSub/${mergeName}"

    tag "${mergeName}"
    label "hyperdrive"
    label "cpu_half"
    label "mem_half"
    label "gpu"

    time 2.hour

    script:
    subVis = "${mergeName}_sub.uvfits"
    """
    ${params.hyperdrive} vis-sub "${params.hypToml}" \
        --data ${data} \
        """ + (
            params.subSrclist ? "--source-list ${params.subSrclist}" : ""
        ) + """ \
        --outputs "${subVis}" \
    """
        // --output-vis-time-average 120s
}

// process hypIonoSub {
//     input:
//     tuple val(mergeName), path(vis)
//     output:
//     tuple val(mergeName), path(subVis), path(json)

//     storeDir "${params.storeDir}/mergePeel/${mergeName}"

//     tag "${mergeName}"
//     label "hyperdrive"
//     label "cpu_half"
//     label "gpu"

//     time 2.hour

//     script:
//     subVis = "${vis.baseName}_ionosub.uvfits"
//     """
//     ${params.hyperdrive} peel \
//         --data "${metafits}" "${vis}" \
//         --source-list "${subSrclist}" \
//         --iono-sub 1000 \
//         --sub 8000 \
//         --outputs "${subVis}" "${json}"
//     """
// }

process wsclean {
    input:
    tuple val(imgName), val(meta), path(vis), val(args_)
    output:
    tuple val(imgName), val(meta), path("wsclean_${imgName}*-{image,model,residual,sources,dirty,psf}.{fits,txt}")

    storeDir "${params.storeDir}/mergeImg/${imgName}"

    tag "${imgName}"
    label "wsclean_hyp_casa"
    label "cpu_full"
    label "mem_full"
    label "nvme"
    time 12.hour

    script:
    vis_ms = vis.collect {"${it.baseName}.ms"}
    vis = vis.collect()
    args = deepcopy(args_)
    if (params.wscleanSaveSrcs) {
        args += " -save-source-list"
    }
    channels = Math.ceil(params.imgWidthHz / params.imgCleanWidthHz).intValue()
    if (meta.startFreqIdx != null && meta.endFreqIdx != null) {
        args += " -channel-range ${meta.startFreqIdx} ${meta.endFreqIdx}"
        channels = Math.min(channels, meta.endFreqIdx - meta.startFreqIdx + 1)
    }
    if (channels > 1) {
        args += " -channels-out ${channels} -join-channels"
        if (params.wscleanSaveSrcs) {
            args += " -fit-spectral-log-pol 2"
        }
    }
    """
    #!/bin/bash -eux
    """ + (
            // convert any uvfits to ms
            [vis, vis_ms].transpose().collect { uv, ms ->
                (uv.extension == 'uvfits' ? \
                """hyperdrive vis-convert --data "${uv}" --outputs "${ms}" -- "${params.hypToml}" """ : "")
            }.join("\n")
    ) + """
    """ + (params.wscleanSaveSrcs?(
        // this is a really really really awful hack to make wsclean
        // think that xx is I so that it can export a srclist
        vis_ms.collect {
            """${params.casa} -c "tb.open('${it}/POLARIZATION', nomodify=False); tb.putcell(rownr=0,columnname='CORR_TYPE',thevalue=[1]); tb.close()" """
        }.join("\n")
    ):"") + """
    ${params.wsclean} \
        -name wsclean_${imgName} \
        -size ${params.imgSize} ${params.imgSize} \
        -scale ${params.imgScale}asec \
        -abs-mem ${0.9 * task.memory.toGiga()} \
        ${args} \
        ${vis_ms.join(" ")}
    """
}

process pbCorrect {
    input:
    tuple val(imgName), val(meta), path(fits)
    output:
    tuple val(imgName), val(meta), path(corrected)

    storeDir "${params.storeDir}/mergeImg/${imgName}"

    tag "${fits.baseName}"

    label "python"

    time {5.min * task.attempt}

    script:
    corrected = "${fits.baseName}._pbeam.fits"
    args = "--fits ${fits} --lst ${meta.lst_rad}"
    template "beam_correct.py"
}

process bane {
    input:
    tuple val(imgName), val(meta), path(fits)
    output:
    tuple val(imgName), val(meta), path(bkg), path(rms)

    storeDir "${params.storeDir}/mergeImg/${imgName}"

    tag "${imgName}.${fits.baseName}"

    label "bane"
    label "cpu_full"
    label "mem_full"

    time {5.min * task.attempt}

    script:
    bkg = "${fits.baseName}_bkg.fits"
    rms = "${fits.baseName}_rms.fits"
    """
    #!/bin/bash
    BANE ${fits} --cores=${task.cpus} --debug
    """
}

process sourceFind {
    input:
    tuple val(imgName), val(meta), path(img), path(bkg), path(rms), path(psf)
    output:
    tuple val(imgName), val(meta), \
        path("aegean_${img.baseName}_comp.fits"), \
        path("aegean_${img.baseName}_comp.reg") // TODO: why does aegean add this crap on the end?

    storeDir "${params.storeDir}/mergeImg/${imgName}"

    tag "${imgName}.${img.baseName}"

    label "aegean"
    label "cpu_half"
    label "mem_full"

    time 10.minute

    script:
    bkg = "${img.baseName}_bkg.fits"
    rms = "${img.baseName}_rms.fits"
    """
    #!/bin/bash
    // BANE ${img} --cores ${task.cpus}
    aegean ${img} \
        --hdu 0 --slice 0 \
        --noise=${rms} --background=${bkg} --cores=${task.cpus} \
        --table aegean_${img.baseName}.fits,aegean_${img.baseName}.reg
    ls -alt
    """
    // high_res = 0.00240241718342004
    // low_res = 0.00266511900788808

    // # ##now run aegean on the low res image, using the higher res table as a prior
    // aegean hyperdrive_things/sub_ts0720-0731_c000-150_minuv-l-200_7200_niter100000-image._pbeam.fits \
    //     --table aegean_sub_ts0720-0731_c000-150_minuv-l-200_7200_niter100000-image._pbeam.fits,aegean_sub_ts0720-0731_c000-150_minuv-l-200_7200_niter100000-image._pbeam.reg \
    //     --cores=20  --hdu 0 --slice 0  \
    //     --noise=low_res_rms.fits --background=low_res_bkg.fits \
    //     --priorized 1 --ratio 1.11 \
    //     --input=aegean_sub_ts0720-0731_c150-300_minuv-l-200_7200_niter100000-image._pbeam_comp_cropped.fits \
    //     --region mimas_reg.mim
}

process imgQuantiles {
    input:
    tuple val(imgName), path(fits)
    output:
    tuple val(imgName), path(csv)

    storeDir "${params.storeDir}/mergeImg/${imgName}"

    tag "${imgName}.${fits.baseName}"

    label "python"

    time {5.min * task.attempt}

    script:
    csv = "quantile_${fits.baseName}.csv"
    template "img_meta.py"
}

process chipskaGrid {
    input:
    tuple val(gridName), val(meta), path(viss)
    output:
    tuple val(gridName), val(meta), path("{crosspower,residpower,residpowerimag,totpower,flagpower,fg_num,outputweights}_xx_${bias_mode}.iter.${gridName}_{${meta.spwNames.join(',')}}.dat")

    storeDir "${params.storeDir}/grid/${gridName}"
    tag "${gridName}"
    label "chipska"
    label "cpu_half"
    label "mem_half"
    label "gpu"

    time 2.hour

    script:
    period = 2500 // to simulate 1000h from 4h of 10s snapshots
    eorband = 1
    eorfield = 0
    maxu = meta.maxu?:600
    nbins = meta.nbins?:80
    nchans = meta.nchans?:150
    bias_mode = 1
    freq_idx_start = 0
    bandwidth = 0
    """
    export DATADIR="\$PWD"
    export OUTPUTDIR="\$PWD/"
    export OMP_NUM_THREADS=${task.cpus}
    """ + (
        // gives {bv,bvdiff,noise,noisediff,weights}_freq*_xx.${gridName}.dat
        viss.collect { vis ->
            """
            /chips/gridvisdska "${vis}" ${meta.gpsTime?:0} ${gridName} ${eorband} -f ${eorfield}
            """
        }.join("\n")
    ) + """
    """ + (
        [meta.startFreqsHz, meta.spwNames].transpose().collect { startFreqHz, spwName ->
            """
            /chips/prepare_ska ${gridName} ${nchans} ${freq_idx_start} 'xx' ${gridName}_${spwName} ${eorband} -p ${period} -c ${meta.chanWidthHz} -n ${startFreqHz}
            /chips/lssa_fg_ska ${gridName}_${spwName} ${nchans} ${nbins} 'xx' ${maxu} ${gridName}_${spwName} ${bias_mode} ${eorband} ${bandwidth}
            """
        }.join("\n")
    ) + """
    ls -al
    """
}

process chipsPlot {
    input:
        tuple val(group), val(meta), path(grid)
    output:
        tuple val(group), val(meta), path("chips${dims}D_${pol}_${suffix}.png")

    storeDir "${params.storeDir}/gridPlot/${group}"

    tag "${group}.${ptype}${endSuffix}"

    label "python"

    time 15.minute

    script:
    pol = meta.pol?:"both"
    if (pol == "both") {
        pol = "xx+yy"
    }
    ptype = meta.ptype?:""
    dims = ptype[0]?:2
    endSuffix = ""
    if (ptype =~ /.*comp/) {
        endSuffix = "_comparison"
    } else if (ptype =~ /.*diff/) {
        endSuffix = "_diff"
    } else if (ptype =~ /.*ratio/) {
        endSuffix = "_ratio"
    } else if (ptype[0] == "2") {
        endSuffix = "_crosspower"
    }
    suffix = ""
    if ((meta.tags?:[])[1]) {
        suffix = "_${meta.tags[1]}${suffix}"
    }
    if ((meta.tags?:[])[0]) {
        suffix = "${meta.tags[0]}${suffix}"
    } else {
        suffix = "${meta.ext}${suffix}"
    }
    basedir = "./"
    args = [
        title: meta.title,
        // file group
        basedir: "./",
        chips_tag: meta.ext,
        chips_tag_one: (meta.tags?:[])[0],
        chips_tag_two: (meta.tags?:[])[1],
        chips_tag_three: (meta.tags?:[])[2],
        // plot group
        plot_type: ptype,
        polarisation: meta.pol,
        min_power: meta.min_power,
        max_power: meta.max_power,
        // plot group 2D
        colourscale: meta.colourscale,
        max_neg_power: meta.max_neg_power,
        min_neg_power: meta.min_neg_power,
        // plot group 1D
        wedge_factor: meta.wedge_factor,
        low_k_edge: meta.low_k_edge,
        high_k_edge: meta.high_k_edge,
        num_k_edges: meta.nbins,
        kperp_max: meta.kperp_max,
        kperp_min: meta.kperp_min,
        kparra_min: meta.kparra_min,
        plot_wedge_cut_2D: meta.plot_wedge_cut_2D,
        chips_tag_one_label: (meta.labels?:[])[0],
        chips_tag_two_label: (meta.labels?:[])[1],
        chips_tag_three_label: (meta.labels?:[])[2],
        // chips group
        lowerfreq: meta.startFreqHz,
        chan_width: meta.chanWidthHz,
        umax: meta.maxu,
        // density_correction: meta.density_correction,
        // omega_matter: meta.omega_matter,
        // omega_baryon: meta.omega_baryon,
        // omega_lambda: meta.omega_lambda,
        // hubble: meta.hubble,

    ].findAll { _, v -> v != null }
        .collect { k, v -> """--${k} "${v}" """ }
        .join(" ")
    if (meta.plot_delta) {
        args += " --plot_delta"
        endSuffix += "_delta"
    }
    if (meta.plot_wedge_cut) {
        args += " --plot_wedge_cut_2D"
    }
    if (meta.nchans) {
        args += " --N_chan ${meta.nchans}"
    }
    if (meta.kperp) {
        args += " --K_perp ${meta.kperp}"
    }

    meta.endSuffix = endSuffix
    suffix = "${suffix}${endSuffix}"

    template "jline_plotchips.py"
}

process ffmpeg {
    input:
        tuple val(name), path("??????????.png"), val(cachebust)

    output:
        tuple path("${name}.mp4"), path("${dot_cachebust}")

    storeDir "${params.storeDir}/video"

    stageInMode "symlink"

    tag "${name}"
    label "ffmpeg"
    module = 'singularity'
    container = "${params.ffmpeg_sif}"

    script:
    dot_cachebust = ".${name}.${cachebust}.cachebust"
    """
    #!/bin/bash -eux
    touch ${dot_cachebust}
    ffmpeg -y -framerate 5 \
        -pattern_type glob -i "??????????.png" \
        -vcodec libx264 \
        -vf "scale='min(3840,iw)':'min(2160,ih)':force_original_aspect_ratio=decrease,pad=ceil(iw/2)*2:ceil(ih/2)*2" \
        -pix_fmt yuv420p \
        "${name}.mp4"
    """
}

@Synchronized
def getSpwName(freqsHz) {
    freqsMHz = freqsHz.collect { sprintf("%03d", (it / 1e6).intValue()) }.sort(false)
    spwName = [freqsMHz[0], freqsMHz[-1]].unique().join("-") + "MHz"
    // print("${freqsMHz[0]}, ${freqsMHz[-1]} -> ${spwName}")
    spwName
}

@Synchronized
def getTimeblockName(timesteps) {
    tsIdxs = timesteps.collect { sprintf("%04d", it.intValue()) }.sort()
    tbName = "ts" + [tsIdxs[0], tsIdxs[-1]].unique().join("-")
    // print("${tsIdxs[0]}, ${tsIdxs[-1]} -> ${tbName}")
    tbName
}

workflow uvfTranspose {
    uvfs = channel.fromPath(params.uvfGlob).map { [it.baseName, it] }
    uvfs | uvMeta
    baseMeta = uvMeta.out.map { baseName, metaJson ->
            meta = parseJson(metaJson)
            times = meta.times?:[]
            freqs = meta.freqs?:[]
            def summary = [
                times: meta.times?:[],
                freqs: meta.freqs?:[]
            ]
            [baseName, summary]
        }
    merges = baseMeta.join(uvfs)
        .flatMap { baseName, meta, uvf ->
            // get spwIndex if spWidthHz is set
            spwIdx = -1
            firstFreqHz = meta.freqs[0]
            if ((params.spWidthHz?:0) > 0) {
                spwStartHz = params.startFreqHz?:0
                spwIdx = ((firstFreqHz - spwStartHz) / params.spWidthHz)
                    .setScale(0, RoundingMode.DOWN).intValue()
            }
            times = meta.times?:[]
            times.withIndex()
                .collate(params.splitTime)
                .collect { time_chunk ->
                    lst_rad = time_chunk[(time_chunk.size/2).intValue()][0].lst_rad
                    time_idxs = time_chunk.collect { it[1] }
                    // print("tb ${tb} lst_rad ${lst_rad}")
                    [spwIdx, time_idxs, lst_rad, firstFreqHz, uvf]
                }
        }
        // deleteme: filter times and channels
        .filter { _spwIdx, ts, lst_rad, fq, _uvf ->
            tb = deepcopy(getTimeblockName(ts))
            spw = deepcopy(getSpwName([fq]))
            filtTb = params.filterTimeblocks?:[tb]
            filtSpw = params.filterSpws?:[spw]
            // print("filter tb ${tb} spw ${spw} filtTb ${filtTb} filtSpw ${filtSpw}")
            (!params.doFilter) || (\
                (filtTb).contains(tb) \
                && (filtSpw).contains(spw)
            )
        }
        .groupTuple(by: 0..2)
        .map { _spwIdx, ts, lst_rad, freqs, mergeUvfs ->
            meta = [
                splitName: getTimeblockName(ts),
                freqs: freqs.target.sort(false) as ArrayList,
                ts: ts.sort(false) as ArrayList,
                lst_rad: lst_rad,
                spwName: getSpwName(freqs)
            ]
            mergeName = "${meta.spwName}_${meta.splitName}"
            [mergeName, deepcopy(meta), mergeUvfs.target as ArrayList]
        }
    merges
        // .take(1) // <- deleteme
        .view { mergeName, _meta, mergeUvfs ->
            "${mergeName} (${mergeUvfs.size}): ${mergeUvfs[0]} -> ${mergeUvfs[-1]}"}
        | uvfMerge
        // | mergeQA
        // | (mergeQA & hypSub)
    // channel.empty() | hypSub

    // calibrate
    if (params.doCal) {
        uvfMerge.out.map { name, meta, vis ->
                [name, deepcopy(meta), vis]
            } | hypCal
    } else {
        channel.empty() | hypCal
    }

    // plot calibration solutions
    hypCal.out.map { name, meta, soln, _model -> [ name, deepcopy(meta), soln ] }
        | plotSols

    // make movies of calibration solutions
    solFrames = plotSols.out.flatMap { name, meta, pngs ->
            coerceList(pngs).collect { png ->
                if (png.baseName =~ /.*_phases.*/) {
                    ["${meta.spwName}_phases", png]
                } else if (png.baseName =~ /.*_amps.*/) {
                    ["${meta.spwName}_amps", png]
                } else {
                    [meta.spwName, png]
                }
            }
        }

    // subtract
    if (params.doSub) {
        if (params.doCal) {
            uvfMerge.out.join(hypCal.out)
                .map { mergeName, meta, uvf, _meta, soln, _model ->
                    ["${mergeName}${params.calSuffix}", deepcopy(meta), [uvf, soln] ]
                }
                | hypSub
        } else {
            uvfMerge.out.map { mergeName, meta, uvf ->
                    ["${mergeName}", deepcopy(meta), [uvf] ]
                }
                | hypSub
        }
    } else {
        channel.empty() | hypSub
    }

    // grid subtracted
    if (params.doGrid) {
        hypSub.out.map { mergeName, meta, vis ->
                def tbName = getTimeblockName(meta.ts)
                def gpsTime = params.firstGpsTime + (meta.ts[0]?:0) * params.intTime
                def chanWidthHz = params.inWidthHz
                def chansPerGrid = (params.gridWidthHz / chanWidthHz).intValue()
                def (startFreqsHz, spwNames) = deepcopy(meta.freqs).unique(false).sort(false).withIndex()
                    .collate(chansPerGrid, false)
                    .collect{ fi_chunk ->
                        startFreqHz = fi_chunk[0][0]
                        endFreqHz = fi_chunk[-1][0]
                        spwName = getSpwName([startFreqHz, endFreqHz])
                        // spwName = sprintf("%03d", (startFreqHz/1e6).intValue())
                        [startFreqHz, deepcopy(spwName)]
                    }
                    .transpose()
                def newMeta = [
                    gpsTime: gpsTime,
                    chanWidthHz: chanWidthHz,
                    nchans: chansPerGrid,
                    startFreqsHz: startFreqsHz,
                    spwNames: spwNames,
                    maxu: params.gridMaxU,
                    nbins: params.gridNKBins,
                    tbName: tbName,
                ]
                ["${tbName}${params.calSuffix}${params.subSuffix}", newMeta, vis]
            }
            | chipskaGrid
    } else {
        channel.empty() | chipskaGrid
    }

    // plot power spectra
    if (params.plotGrids) {
        def chipsGrids = chipskaGrid.out.flatMap { gridName, meta, grids ->
                [meta.spwNames, meta.startFreqsHz].transpose().collect { spwName, startFreqHz ->
                    spwGrids = grids.findAll { grid -> grid.baseName =~ /.*${gridName}_${spwName}.*/ }
                    newMeta = [
                        startFreqHz: startFreqHz,
                        spwName: spwName,
                    ]
                    ["${gridName}_${spwName}", deepcopy(meta) + newMeta, spwGrids]
                }
            }
            .filter { gridName, meta, grids ->
                // filter out grids with no data
                if (spwGrids.size == 0) {
                    print("no grids found for spwName=${gridName} in ${grids}")
                }
                grids.size > 0
            }

        if (params.plotSingle1DDelta) {
            gridPlotsSingle1DDelta = chipsGrids.map { gridName, meta, grid ->
                    def spwName = meta.spwNames[0]
                    def startFreqHz = meta.startFreqsHz[0]
                    def newMeta = [
                        ext: gridName,
                        ptype: "1D",
                        pol: "xx",
                        plot_delta: true,
                        min_power: 1e-3,
                        max_power: 1e+15,
                        title: "crosspower\\n${gridName}",
                        plot_name: 'chips1d',
                        tags: [],
                    ]
                    [gridName, deepcopy(meta) + newMeta, grid]
                }
        } else {
            gridPlotsSingle1DDelta = channel.empty()
        }

        if (params.plotSingle2D) {
            gridPlotsSingle2D = chipsGrids.map { gridName, meta, grid ->
                    def spwName = meta.spwNames[0]
                    def startFreqHz = meta.startFreqsHz[0]
                    def newMeta = [
                        ext: gridName,
                        ptype: "2D",
                        pol: "xx",
                        min_power: 1e-4,
                        max_power: 1e+2,
                        title: "crosspower\\n${gridName}",
                        plot_name: 'chips2d',
                        tags: [],
                    ]
                    [gridName, deepcopy(meta) + newMeta, grid]
                }
        } else {
            gridPlotsSingle2D = channel.empty()
        }

        gridPlotsSingle1DDelta.mix(gridPlotsSingle2D) | chipsPlot
    } else {
        channel.empty() | chipsPlot
    }

    // make movies of grid plots
    chipsFrames = chipsPlot.out.flatMap { name, meta, pngs ->
            coerceList(pngs).collect { png ->
                if (png.baseName =~ /.*${meta.ext}.*/) {
                    tokens = png.baseName.split("${meta.ext}")
                    prefix = tokens[0]
                    suffix = tokens[-1]
                    ["${prefix}${meta.spwName}${suffix}", png]
                } else {
                    throw new Exception("meta.ext=${meta.ext} not found in png.baseName=${png.baseName}")
                    [meta.spwName, png]
                }
            }
        }


    // image
    if (params.imgUnsub) {
        mergeImgs = uvfMerge.out.flatMap {mergeName, meta, vis ->
                deepcopy(meta.freqs).unique(false).sort(false).withIndex()
                    .collate((params.imgWidthHz / params.inWidthHz).intValue(), false)
                    // .take(1) // deleteme
                    .collect { fi_chunk ->
                        startFreqHz = fi_chunk[0][0]
                        startFreqIdx = fi_chunk[0][-1]
                        endFreqHz = fi_chunk[-1][0]
                        endFreqIdx = fi_chunk[-1][-1]
                        spwName = getSpwName([startFreqHz, endFreqHz])
                        wscArgs = "${params.wscleanCommon} ${params.wscleanMsw} ${params.wscleanMswExtra}"
                        newMeta = [
                            spwName:spwName,
                            startFreqIdx:startFreqIdx,
                            endFreqIdx:endFreqIdx
                        ]
                        ["${spwName}_${meta.splitName}_msw", deepcopy(meta) + newMeta, vis, wscArgs]
                    }
            }
    } else {
        mergeImgs = channel.empty()
    }
    if (params.imgModel) {
        modelImgs = hypCal.out.flatMap {mergeName, meta, _soln, model ->
                deepcopy(meta.freqs).unique(false).sort(false).withIndex()
                    .collate((params.imgWidthHz / params.inWidthHz).intValue(), false)
                    // .take(1) // deleteme
                    .collect { fi_chunk ->
                        startFreqHz = fi_chunk[0][0]
                        startFreqIdx = fi_chunk[0][-1]
                        endFreqHz = fi_chunk[-1][0]
                        endFreqIdx = fi_chunk[-1][-1]
                        spwName = getSpwName([startFreqHz, endFreqHz])
                        wscArgs = "${params.wscleanCommon} ${params.wscleanMsw} ${params.wscleanMswExtra}"
                        newMeta = [
                            spwName:spwName,
                            startFreqIdx:startFreqIdx,
                            endFreqIdx:endFreqIdx
                        ]
                        ["${spwName}_${meta.splitName}_model_msw", deepcopy(meta) + newMeta, model, wscArgs]
                    }
            }
    } else {
        modelImgs = channel.empty()
    }
    if (params.imgSub) {
        subImgs = hypSub.out.flatMap {mergeName, meta, vis ->
                deepcopy(meta.freqs).unique(false).sort(false).withIndex()
                    .collate((params.imgWidthHz / params.inWidthHz).intValue(), false)
                    // .take(1) // deleteme
                    // .findAll { it[-1][-1] != it[0][-1] }
                    .collect { fi_chunk ->
                        startFreqHz = fi_chunk[0][0]
                        startFreqIdx = fi_chunk[0][-1]
                        endFreqHz = fi_chunk[-1][0]
                        endFreqIdx = fi_chunk[-1][-1]
                        spwName = getSpwName([startFreqHz, endFreqHz])
                        wscArgs = "${params.wscleanCommon} ${params.wscleanMsw} ${params.wscleanMswExtra}"
                        newMeta = [
                            spwName:spwName,
                            startFreqIdx:startFreqIdx,
                            endFreqIdx:endFreqIdx
                        ]
                        ["${spwName}_${meta.splitName}${params.calSuffix?:''}${params.subSuffix?:''}_msw", deepcopy(meta) + newMeta, vis, wscArgs]
                    }
            }
    } else {
        subImgs = channel.empty()
    }

    mergeImgs.mix(subImgs).mix(modelImgs)
        | wsclean

    // break out image products into images and psfs
    imgPsfs = wsclean.out.map { imgName, meta, fitss ->
            imgs = fitss.findAll { it.name =~ /.*-image.fits/ }
            psfs = fitss.findAll { it.name =~ /.*-psf.fits/ }
            if (imgs.size != 1 && psfs.size != 1) {
                imgs = fitss.findAll { it.name =~ /.*-MFS-image.fits/ }
                psfs = fitss.findAll { it.name =~ /.*-MFS-psf.fits/ }
            }
            assert imgs.size == 1 && psfs.size == 1
            [imgName, meta, imgs[0], psfs[0]]
        }
    imgPsfs.map { imgName, meta, img, psf -> [imgName, meta, img] }
        | pbCorrect
        | bane

    pbCorrect.out
        .join(bane.out.map { imgName, meta, bkg, rms -> [imgName, bkg, rms] })
        .join(imgPsfs.map { imgName, meta, img, psf -> [imgName, psf] })
        | sourceFind

    solFrames.mix(chipsFrames)
        .groupTuple()
        .map { name, pngs_ ->
            pngs = coerceList(pngs_)
            def latest = pngs.collect { it.lastModified() }.max()
            def cachebust = "${latest}_x" + sprintf("%04d", pngs.size())
            [name, pngs.sort(), cachebust]
        }
        | ffmpeg

    // wsclean.out.flatMap { imgName, files ->
    //         def (spwName, tbName, subName, imgType) = imgName.split('_')
    //         files.findAll { it.baseName =~ /.*MFS-image/ }.collect { fits ->
    //             ["${spwName}_${subName}_${imgType}", fits]
    //         }
    //     }
}

// workflow msTranspose {
//     uv = channel.fromPath(params.msGlob, type: 'dir').map { [it.baseName, it] }
//     uv | msMeta
//     baseMeta = msMeta.out.map { baseName, metaJson ->
//             meta = parseJson(metaJson)
//             times = meta.times?:[]
//             freqs = meta.freqs?:[]
//             def summary = [
//                 times: meta.times?:[],
//                 freqs: meta.freqs?:[]
//             ]
//             [baseName, summary]
//         }
//     baseMeta.join(uv)
//         .flatMap { baseName, meta, uv ->
//             times = meta.times?:[]
//             (0..<times.size).collate(params.splitTime)
//                 .withIndex()
//                 .take(1) // <- deleteme
//                 .collect { ts, i ->
//                     splitName = sprintf("%04d", i)
//                     [baseName, splitName, uv, ts]
//                 }
//         }
//         | msSplit
//     // now split bandwidth into chunks of params.spWidthHz
//     merges = baseMeta.flatMap { _, meta -> [meta.freqs?:[]] }
//         .collect()
//         .map { freqs -> freqs.min() }
//         .combine( baseMeta.join(msSplit.out) )
//         .map { startFreqHz, baseName, meta, splitName, uv ->
//             firstFreqHz = meta.freqs[0]
//             spwIdx = 0
//             if (params.spWidthHz > 0) {
//                 spwIdx = ((firstFreqHz - startFreqHz) / params.spWidthHz)
//                     .setScale(0, RoundingMode.DOWN).intValue()
//             }
//             [spwIdx, splitName, firstFreqHz, uv]
//         }
//         .groupTuple(by: 0..1).map { _, splitName, freqs, uv ->
//             freqs.sort()
//             spwStartMHz = sprintf("%03d", (freqs[0] / 1e6).intValue())
//             spwEndMHz = sprintf("%03d", (freqs[-1] / 1e6).intValue())
//             spwName = spwStartMHz == spwEndMHz ? "${spwStartMHz}MHz" : "${spwStartMHz}-${spwEndMHz}MHz"
//             ["${spwName}_${splitName}", uv.sort()]
//         }
//         | view { chunkName, vis -> "${chunkName} (${vis.target.size}): ${vis.target[0]} -> ${vis.target[-1]}"}
//     merges
//         | msMerge
//         | hypCal
//         | plotSols
//     mergeCounts = merges.map { mergeName, vis -> [mergeName, vis.target.size] }
//     msMerge.out.join(mergeCounts)
//         | wsclean
// }
