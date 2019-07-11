#!/usr/bin/env nextflow

Channel
    .fromPath(params.input)
    .splitCsv(header:true, sep:'\t', quote:'"')
    .map{ row-> tuple(row.sample_id, file(row.read_1), file(row.read_2), ( row.read_l == "" ? null : file(row.read_l))) }
    .set { samples_fastqc_ch; }


def summary = [:]
summary['Pipeline Name']  = 'trim-filter-unicycler'
summary['Input']          = params.input
summary['Output dir']     = params.outdir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile


summary.each{ k, v -> println "${k}: ${v}" }


/*
 * FastQC Pre-Trim
 */
process fastqc_pre_trim {
    tag "$sampleId"
    cpus 8
    conda '/home/dfornika/miniconda3/envs/fastqc-0.11.8'

    input:
    set sampleId, file(read_1), file(read_2), file(read_l) from samples_fastqc_ch

    output:
	file "*_R{1,2}_*_fastqc.zip" into illumina_fastqc_results_ch
        file "*MinION*fastqc.zip" optional true into minion_fastqc_results_ch

    script:
    """
    fastqc -q -t 8 *.fastq.gz
    """
}

/*
 * MultiQC Illumina
 */
process multiqc_illumina_pre_trim {
    cpus 8
    conda '/home/dfornika/miniconda3/envs/multiqc-1.7'
    publishDir "${params.outdir}", mode: 'copy', pattern: "multiqc_report_illumina.html"

    input:
	file '*_fastqc.zip' from illumina_fastqc_results_ch.collect()


    script:
    """
    multiqc -n multiqc_report_illumina.html .
    """
}

/*
 * MultiQC MinION
 */
process multiqc_minion_pre_trim {
    cpus 8
    conda '/home/dfornika/miniconda3/envs/multiqc-1.7'
    publishDir "${params.outdir}", mode: 'copy', pattern: "multiqc_report_minion.html"

    input:
	file '*_fastqc.zip' from minion_fastqc_results_ch.collect()


    script:
    """
    multiqc -n multiqc_report_minion.html .
    """
}

/*
 * trim with trim_galore, run fastqc on trimmed files.

process trim_illumina {
    tag "$sampleId"
    cpus 4
    memory '1 GB'
    conda 'trim-galore=0.6.3'
    input:
        set val(sampleId), file(read_1), file(read_2), file(read_l) from samples_trim_ch
    when:
    false
    output:
        set val("$sampleId"), file("*_val_1.fq"), file("*_val_2.fq") into samples_trimmed_ch
        file "*fastqc*"
    script:
    """
    trim_galore --paired $read1 $read2
    """
}
*/




