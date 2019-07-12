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
 * FastQC (Pre-Trimming)
 */
process fastqc_pre_trim {
    tag "$sample_id"
    cpus 8
    conda '/home/dfornika/miniconda3/envs/fastqc-0.11.8'

    input:
    set sample_id, file(read_1), file(read_2), file(read_l) from samples_fastqc_ch

    output:
    file "*_R{1,2}_*_fastqc.zip" into illumina_fastqc_results_ch
    file "*MinION*fastqc.zip" optional true into minion_fastqc_results_ch
    set sample_id, file(read_1), file(read_2) into trim_illumina_ch
    set sample_id, file(read_l) into trim_minion_ch
    
    script:
    """
    fastqc -q -t 8 *.fastq.gz
    """
}

/*
 * MultiQC Illumina (Pre-Trimming)
 */
process multiqc_illumina_pre_trim {
    cpus 8
    conda '/home/dfornika/miniconda3/envs/multiqc-1.7'
    publishDir "${params.outdir}", mode: 'copy', pattern: "multiqc_report_illumina.html"

    input:
    file '*_fastqc.zip' from illumina_fastqc_results_ch.collect()

    output:
    file 'multiqc_report_illumina.html'
    
    script:
    """
    multiqc -n multiqc_report_illumina.html .
    """
}

/*
 * MultiQC MinION (Pre-Trimming)
 */
process multiqc_minion_pre_trim {
    cpus 8
    conda '/home/dfornika/miniconda3/envs/multiqc-1.7'
    publishDir "${params.outdir}", mode: 'copy', pattern: "multiqc_report_minion.html"

    input:
    file '*_fastqc.zip' from minion_fastqc_results_ch.collect()

    output:
    file 'multiqc_report_minion.html'
    
    script:
    """
    multiqc -n multiqc_report_minion.html .
    """
}

/*
 * Trim with illumina reads with trim_galore.
 */
process trim_illumina {
    tag "$sample_id"
    cpus 4
    memory '2 GB'
    conda '/home/dfornika/miniconda3/envs/trim-galore-0.6.3'
    input:
    set sample_id, file(read_1), file(read_2) from trim_illumina_ch

    output:
    set file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into illumina_trimmed_ch

    script:
    """
    trim_galore --cores 4 --basename $sample_id --paired $read_1 $read_2
    """
}

/*
 * Trim minion reads with porechop.
 */
process trim_minion {
    tag "$sample_id"
    cpus 16
    memory '2 GB'
    conda '/home/dfornika/miniconda3/envs/porechop-0.2.3'

    input:
    set sample_id, file(read_l) from trim_minion_ch

    when:
    !(read_l.name =~ /^input.\d/) 
    
    output:
    file("*.trim.fastq.gz") into minion_trimmed_ch

    script:
    """
    porechop --threads 16  -i $read_l -o ${sample_id}.MinION.trim.fastq.gz
    """
}

/*
 * FastQC (Post-Trimming)
 */
process fastqc_post_trim {
    tag "$sample_id"
    cpus 8
    conda '/home/dfornika/miniconda3/envs/fastqc-0.11.8'

    input:
    set file(read_1), file(read_2) from illumina_trimmed_ch
    file read_l from minion_trimmed_ch

    output:
    file "*_R{1,2}_*_fastqc.zip" into illumina_trimmed_fastqc_results_ch
    file "*MinION*fastqc.zip" into minion_trimmed_fastqc_results_ch

    script:
    """
    fastqc -q -t 8 *.fastq.gz *.fq.gz
    """
}

/*
 * MultiQC Illumina (Post-Trimming)
 */
process multiqc_illumina_post_trim {
    cpus 8
    conda '/home/dfornika/miniconda3/envs/multiqc-1.7'
    publishDir "${params.outdir}", mode: 'copy', pattern: "multiqc_report_illumina_post_trimmming.html"

    input:
    file '*_fastqc.zip' from illumina_trimmed_fastqc_results_ch.collect()

    output:
    file 'multiqc_report_illumina_post_trimming.html'

    script:
    """
    multiqc -n multiqc_report_illumina_post_trimming.html .
    """
}

/*
 * MultiQC MinION (Post-Trimming)
 */
process multiqc_minion_post_trim {
    cpus 8
    conda '/home/dfornika/miniconda3/envs/multiqc-1.7'
    publishDir "${params.outdir}", mode: 'copy', pattern: "multiqc_report_minion_post_trimming.html"

    input:
    file '*_fastqc.zip' from minion_trimmed_fastqc_results_ch.collect()

    output:
    file 'multiqc_report_minion_post_trimming.html'

    script:
    """
    multiqc -n multiqc_report_minion_post_trimming.html .
    """
}
