nextflow.enable.dsl=2

params.input = "raw_fastq"         // Input directory
params.output = "vsearch_output"  // Output directory
params.threads = 16               // Number of threads

// Flags to skip steps
params.skip_fastp = false
params.skip_merge = false
params.skip_filter = false
params.skip_concat = false
params.skip_derep = false
params.skip_denoise = false
params.skip_chimera = false
params.skip_swarm = false
params.skip_singleton = false
params.skip_count = false


// Step 1: FastP QC
process fastpQC {
    input:
    path input_dir

    output:
    path "${params.output}/fastp"

    script:
    """
    mkdir -p ${params.output}/fastp
    for R1 in ${input_dir}/*_R1.fastq.gz; do
        SAMPLE=\$(basename \${R1} _R1.fastq.gz)
        R2="${input_dir}/\${SAMPLE}_R2.fastq.gz"

        if [[ -f "\${R2}" ]]; then
            fastp -i \${R1} -I \${R2} \
                  -o ${params.output}/fastp/\${SAMPLE}_R1.fastq.gz \
                  -O ${params.output}/fastp/\${SAMPLE}_R2.fastq.gz \
                  -j ${params.output}/fastp/\${SAMPLE}_report.json \
                  -h ${params.output}/fastp/\${SAMPLE}_report.html \
                  -w ${params.threads}
        fi
    done
    """
}

// Step 2: Merge paired-end reads
process mergeReads {
    input:
    path "${params.output}/fastp"

    output:
    path "${params.output}/merged"

    script:
    """
    mkdir -p ${params.output}/merged
    for R1 in ${params.output}/fastp/*_R1.fastq.gz; do
        SAMPLE=\$(basename \${R1} _R1.fastq.gz)
        R2="${params.output}/fastp/\${SAMPLE}_R2.fastq.gz"

        if [[ -f "\${R2}" ]]; then
            vsearch --fastq_mergepairs \${R1} --reverse \${R2} \
                    --fastqout ${params.output}/merged/\${SAMPLE}.merged.fastq \
                    --threads ${params.threads}
        fi
    done
    """
}

// Step 3: Filter merged reads
process filterReads {
    input:
    path "${params.output}/merged"

    output:
    path "${params.output}/filtered"

    script:
    """
    mkdir -p ${params.output}/filtered
    for M in ${params.output}/merged/*.merged.fastq; do
        SAMPLE=\$(basename \${M} .merged.fastq)
        vsearch --fastx_filter \${M} \
                --fastq_maxee 1.0 \
                --fastq_minlen 200 \
                --fastaout ${params.output}/filtered/\${SAMPLE}.filtered.fasta
    done
    """
}

// Step 4: Concatenate reads
process concatenateReads {
    input:
    path "${params.output}/filtered"

    output:
    path "${params.output}/concat"

    script:
    """
    mkdir -p ${params.output}/concat
    rm -f ${params.output}/concat/concat.fasta
    for F in ${params.output}/filtered/*.filtered.fasta; do
        SAMPLE=\$(basename \${F} .filtered.fasta)
        HEADER=\$(grep '>' \${F} | head -n 1)
        NEW_HEADER=">\${SAMPLE}"
        sed "s/>.*/\${NEW_HEADER}/" \${F} >> ${params.output}/concat/concat.fasta
    done
    """
}

// Step 5: Dereplicate sequences
process dereplicateSequences {
    input:
    path "${params.output}/concat"

    output:
    path "${params.output}/derep"

    script:
    """
    mkdir -p ${params.output}/derep
    vsearch --derep_fulllength ${params.output}/concat/concat.fasta \
            --output ${params.output}/derep/derep.fasta \
            --sizeout --threads ${params.threads}
    """
}

// Step 6: ASV Denoising
process denoiseASV {
    input:
    path "${params.output}/derep"

    output:
    path "${params.output}/denoise"

    script:
    """
    mkdir -p ${params.output}/denoise
    vsearch --cluster_unoise ${params.output}/derep/derep.fasta \
            --centroids ${params.output}/denoise/centroids.fasta \
            --threads ${params.threads}
    """
}

// Step 7: Chimera checking
process chimeraCheck {
    input:
    path "${params.output}/denoise"

    output:
    path "${params.output}/nochimeras"

    script:
    """
    mkdir -p ${params.output}/nochimeras
    vsearch --uchime3_denovo ${params.output}/denoise/centroids.fasta \
            --nonchimeras ${params.output}/nochimeras/nochimeras.fasta \
            --threads ${params.threads}
    """
}

// Step 8: Swarm Clustering
process swarmClustering {
    input:
    path "${params.output}/nochimeras"

    output:
    path "${params.output}/swarm"

    script:
    """
    mkdir -p ${params.output}/swarm
    swarm -d 1 -f -z -t ${params.threads} \
          -o ${params.output}/swarm/swarms.txt \
          ${params.output}/nochimeras/nochimeras.fasta
    """
}

// Step 9: Remove singletons
process removeSingletons {
    input:
    path "${params.output}/nochimeras"

    output:
    path "${params.output}/ASVs"

    script:
    """
    mkdir -p ${params.output}/ASVs
    vsearch --sortbysize ${params.output}/nochimeras/nochimeras.fasta \
            --output ${params.output}/ASVs/ASVs.fasta \
            --minsize 2 --threads ${params.threads}
    """
}

// Step 10: Create ASV count matrix
process createCountMatrix {
    input:
    path "${params.output}/ASVs"

    output:
    path "${params.output}/ASVs"

    script:
    """
    vsearch --usearch_global ${params.output}/concat/concat.fasta \
            --db ${params.output}/ASVs/ASVs.fasta \
            --otutabout ${params.output}/ASVs/ASV_counts.tsv \
            --threads ${params.threads}
    """
}
workflow {
    if (!params.skip_fastp) {
        fastpQC(input_dir: params.input)
    }
    if (!params.skip_merge) {
        mergeReads(input_dir: "${params.output}/fastp")
    }
    if (!params.skip_filter) {
        filterReads(input_dir: "${params.output}/merged")
    }
    if (!params.skip_concat) {
        concatenateReads(input_dir: "${params.output}/filtered")
    }
    if (!params.skip_derep) {
        dereplicateSequences(input_file: "${params.output}/concat/concat.fasta")
    }
    if (!params.skip_denoise) {
        denoiseASV(input_file: "${params.output}/derep/derep.fasta")
    }
    if (!params.skip_chimera) {
        chimeraCheck(input_file: "${params.output}/denoise/centroids.fasta")
    }
    if (!params.skip_swarm) {
        swarmClustering(input_file: "${params.output}/nochimeras/nochimeras.fasta")
    }
    if (!params.skip_singleton) {
        removeSingletons(input_file: "${params.output}/nochimeras/nochimeras.fasta")
    }
    if (!params.skip_count) {
        createCountMatrix(input_dir: "${params.output}/concat", db_file: "${params.output}/ASVs/ASVs.fasta")
    }
}
