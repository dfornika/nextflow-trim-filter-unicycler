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
    set sample_id, file(read_l) into filt_minion_ch
    
    script:
    """
    fastqc -q -t 8 *.fastq.gz
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
    set file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into illumina_trimmed_fastqc_ch
    set sample_id, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into illumina_trimmed_unicycler_ch

    script:
    """
    trim_galore --cores 4 --basename $sample_id --paired $read_1 $read_2
    """
}


/*
 * Filter minion reads with filtlong.
 */
process filtlong_minion {
    tag "$sample_id"
    cpus 4
    memory '2 GB'
    conda '/home/dfornika/miniconda3/envs/filtlong-0.2.0'

    input:
    set sample_id, file(read_l) from filt_minion_ch

    when:
    !(read_l.name =~ /^input.\d/) 
    
    output:
    set sample_id, file("*.filt.fastq.gz") into trim_minion_ch
    
    script:
    """
    filtlong  \
    --min_length 1000 \
    --keep_percent 90 \
    --length_weight 0.5 \
    --mean_q_weight 10 \
    --target_bases 500000000 \
    $read_l | gzip > ${sample_id}.MinION.filt.fastq.gz
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
    file("*.trim.fastq.gz") into minion_trimmed_filtered_fastqc_ch
    set sample_id, file("*.trim.fastq.gz") into minion_trimmed_filtered_unicycler_ch
    script:
    """
    porechop --threads 16  -i $read_l -o ${sample_id}.MinION.filt.trim.fastq.gz
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
    set file(read_1), file(read_2) from illumina_trimmed_fastqc_ch
    file read_l from minion_trimmed_filtered_fastqc_ch

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
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.html"

    input:
    file '*_fastqc.zip' from illumina_fastqc_results_ch.mix(illumina_trimmed_fastqc_results_ch).collect()

    output:
    file '*.html'

    script:
    """
    multiqc -n multiqc_report_illumina.html .
    """
}

/*
 * MultiQC MinION (Post-Trimming)
 */
process multiqc_minion_post_trim {
    cpus 8
    conda '/home/dfornika/miniconda3/envs/multiqc-1.7'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.html"

    input:
    file '*_fastqc.zip' from minion_fastqc_results_ch.mix(minion_trimmed_fastqc_results_ch).collect()

    output:
    file '*.html'

    script:
    """
    multiqc -n multiqc_report_minion.html .
    """
}

illumina_trimmed_unicycler_ch

/*
 * Assemble with Unicycler
 */
process unicycler_assemble {
    cpus 16
    conda '/home/dfornika/miniconda3/envs/unicycler-0.4.7'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_assembly.fasta"

    input:
    set sample_id, file(read_1), file(read_2), file(read_l) from illumina_trimmed_unicycler_ch.join(minion_trimmed_filtered_unicycler_ch, remainder: true)

    output:
    file '*_assembly.fasta'

    script:
    if( read_l.name =~ /^input.\d/ )
        """
        unicycler \
        --threads 16 \
        -1 $read_1 \
        -2 $read_2 \
        -o .
        cp assembly.fasta ${sample_id}_assembly.fasta
        """
    else
	"""
        unicycler \
        --threads 16 \
        -1 $read_1 \
        -2 $read_2 \
        -l $read_l \
        -o .
        cp assembly.fasta ${sample_id}_assembly.fasta
        """
}


