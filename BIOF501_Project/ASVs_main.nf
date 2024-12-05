nextflow.enable.dsl=2

// Define parameters  
params.threads = 4 // Number of threads to use for parallel processing
params.input = './raw_fastq' // Directory containg the raw FASTQ files for the pipeline
params.outdir = './results' // Directory where the output files will be saved 

// Channel for FASTQ file listing
Channel
    .fromFilePairs("${params.input}/*_{1,2}.fastq", flat: true) // Matches the paired FASTQ files
    .set { reads_ch } // Sets the channel variable for the workflow

// Process 1: FastP QC
process fastpQC {
    tag { sample_id } // Tags the fastQC process with the corresponding sample ID for logging

    input:
    tuple val(sample_id), path(read1), path(read2) // Takes input paired FASTQ files with sample ID

    output:
    tuple val(sample_id),
          path("${sample_id}_1.fastp.fastq"), // Cleaned forward read
          path("${sample_id}_2.fastp.fastq") // Cleaned reverse read

    publishDir "${params.outdir}/fastpQC", mode: 'copy' // Saves the output to its own directory within "results"

    script:
    """
    fastp -i ${read1} -I ${read2} \\
          -o ${sample_id}_1.fastp.fastq \\
          -O ${sample_id}_2.fastp.fastq \\
          -w ${params.threads}
    """
}

// Process 2: Merge paired-end reads
process mergeReads {
    tag { sample_id } // Tags the process with the sample ID for logging

    input:
    tuple val(sample_id), path(R1), path(R2) // Input cleaned forward and reverse reads

    output:
    tuple val(sample_id), path("${sample_id}.merged.fastq") // Merged output file

    publishDir "${params.outdir}/mergeReads", mode: 'copy' // Saves the output to its own directory within "results"


    script:
    """
    vsearch --fastq_mergepairs ${R1} --reverse ${R2} \\
            --fastqout ${sample_id}.merged.fastq \\
            --threads ${params.threads}
    """
}

// Process 3: Filter merged reads
process filterReads {
    tag { sample_id } // Tags the process with the sample ID for logging

    input:
    tuple val(sample_id), path(merged_fastq) // Input merged FASTQ file

    output:
    tuple val(sample_id), path("${sample_id}.filtered.fasta") // Filtered FASTA file

    publishDir "${params.outdir}/filterReads", mode: 'copy' // Saves the output to its own directory within "results"
 
    script:
    """
    vsearch --fastx_filter ${merged_fastq} \\
            --fastq_maxee 1.0 \\
            --fastaout ${sample_id}.filtered.fasta
    """
}

// Process 4: Concatenate reads
process concatenateReads {
    input:
    path(filtered_fastas) // Input filtered FASTA files

    output:
    path("concat.fasta") // Output concatenated FASTA file

    publishDir "${params.outdir}/concatenateReads", mode: 'copy' // Saves the output to its own directory within "results"

    script:
    """
    cat ${filtered_fastas} > concat.fasta
    """
}

// Process 5: Dereplicate sequences
process dereplicateSequences {
    input:
    path("concat.fasta") // Input concatenated FASTA file

    output:
    path("derep.fasta") // Output dereplicated FASTA file

    publishDir "${params.outdir}/dereplicateSequences", mode: 'copy' // Saves the output to its own directory within "results"

    script:
    """
    vsearch --derep_fulllength concat.fasta \\
            --output derep.fasta \\
            --sizeout --threads ${params.threads}
    """
}

// Process 6: ASV Denoising
process denoiseASV {
    input:
    path("derep.fasta") // Input dereplicated FASTA file

    output:
    path("centroids.fasta") // Output ASV centroids

    publishDir "${params.outdir}/denoiseASV", mode: 'copy' // Saves the output to its own directory within "results"

    script:
    """
    vsearch --cluster_unoise derep.fasta \\
            --centroids centroids.fasta \\
            --threads ${params.threads}
    """
}

// Process 7: Chimera checking
process chimeraCheck {
    input:
    path("centroids.fasta") // Input ASV centroids file

    output: 
    path("nochimeras.fasta") // Output non-chimeric ASVs

    publishDir "${params.outdir}/chimeraCheck", mode: 'copy' // Saves the output to its own directory within "results"

    script:
    """
    vsearch --uchime3_denovo centroids.fasta \\
            --nonchimeras nochimeras.fasta \\
            --threads ${params.threads}
    """
}

// Process 8: Remove singletons
process removeSingletons {
    input:
    path("nochimeras.fasta") // Input non-chimeric ASVs

    output:
    path("ASVs.fasta") // Output ASVs with no singletons

    publishDir "${params.outdir}/removeSingletons", mode: 'copy' // Saves the output to its own directory within "results"

    script:
    """
    vsearch --sortbysize nochimeras.fasta \\
            --output ASVs.fasta \\
            --minsize 2 --threads ${params.threads}
    """
}

// Process 9: Create ASV count matrix
process createCountMatrix {
    input:
    path("concat.fasta") // Input concatenated FASTA
    path("ASVs.fasta") // Input ASVs file 

    output:
    path("ASV_counts.tsv") // Output ASV count matrix

    publishDir "${params.outdir}/createCountMatrix", mode: 'copy' // Saves the output to its own directory within "results"

    script:
    """
    vsearch --usearch_global concat.fasta \\
            --db ASVs.fasta \\
            --id 0.99 \\
            --otutabout ASV_counts.tsv \\
            --threads ${params.threads}
    """
}

// Workflow definition
workflow {
    // Step 1: FastP QC
    cleaned_reads_ch = fastpQC(reads_ch)

    // Step 2: Merge paired-end reads
    merged_reads_ch = mergeReads(cleaned_reads_ch)

    // Step 3: Filter merged reads
    filtered_reads_ch = filterReads(merged_reads_ch)

    // Step 4: Concatenate reads
    concat_fasta = concatenateReads(filtered_reads_ch.map { it[1] })

    // Step 5: Dereplicate sequences
    derep_fasta = dereplicateSequences(concat_fasta)

    // Step 6: ASV Denoising
    denoise_fasta = denoiseASV(derep_fasta)

    // Step 7: Chimera checking
    nochimeras_fa = chimeraCheck(denoise_fasta)

    // Step 8: Remove singletons
    asvs_fa = removeSingletons(nochimeras_fa)

    // Step 9: Create ASV count matrix
    asv_counts = createCountMatrix(concat_fasta, asvs_fa)
}
