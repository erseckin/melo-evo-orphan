nextflow.enable.dsl=2

// Reuse your existing modules
include { indexation_genome_STAR ; mapping_STAR } from './module/mapping/main.nf'
include { rsem_index ; rsem_calculate ; matrix_entete_rsem ;
          matrix_colonne_rsem_TPM ; matrix_rsem as matrix_rsem_TPM ;
          matrix_colonne_rsem_FPKM ; matrix_rsem as matrix_rsem_FPKM } from './module/count/main.nf'

workflow {
    /* --- Inputs --- */
    // genome + annotation as staged files (portable across nodes/partitions)
    genome_ch = Channel.fromPath(params.genome)
    anno_ch   = Channel.fromPath(params.annotation).collect()

    // SINGLE-END reads: e.g. "/path/*.fastq.gz"
    reads_se = Channel
        .fromPath(params.reads)
        .map { fq -> tuple(fq.baseName, fq) }   // (sample_id, fastq)

    /* --- STAR index + mapping --- */
    indexation_genome_STAR(genome_ch)
    mapping_STAR(indexation_genome_STAR.out.index.collect(), reads_se, anno_ch)

    /* --- RSEM index + quantify (single-end, from STAR transcriptome BAM) --- */
    rsem_index(genome_ch, Channel.fromPath(params.annotation), indexation_genome_STAR.out.index.collect())
    rsem_calculate(mapping_STAR.out.transcri_bam, rsem_index.out.rsem_index.collect())

    /* --- Build TPM / FPKM matrices only --- */
    name_matrix_rsem_T = Channel.fromPath('matrice_rsem_TPM.txt')
    name_matrix_rsem_F = Channel.fromPath('matrice_rsem_FPKM.txt')

    // header
    matrix_entete_rsem(rsem_calculate.out.rsem_results)
    rsem_header = matrix_entete_rsem.out.entete

    // TPM matrix
    matrix_colonne_rsem_TPM(rsem_calculate.out.rsem_results)
    tpm_cols = matrix_colonne_rsem_TPM.out.matrix_colonne_tmp
                 .toSortedList( { a, b -> a[0] <=> b[0] } )
                 .flatten()
                 .filter { it =~ /tmp/ }
                 .collect()
    matrix_rsem_TPM(tpm_cols, rsem_header, name_matrix_rsem_T)

    // FPKM matrix
    matrix_colonne_rsem_FPKM(rsem_calculate.out.rsem_results)
    fpkm_cols = matrix_colonne_rsem_FPKM.out.matrix_colonne_tmp
                  .toSortedList( { a, b -> a[0] <=> b[0] } )
                  .flatten()
                  .filter { it =~ /tmp/ }
                  .collect()
    matrix_rsem_FPKM(fpkm_cols, rsem_header, name_matrix_rsem_F)
}

