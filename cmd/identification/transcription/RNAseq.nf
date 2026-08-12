nextflow.enable.dsl=2


include { indexation_genome_STAR ; mapping_STAR } from './module/mapping/main.nf'
include { index_samtools } from './module/samtools/main.nf'
include { bigwig } from './module/bigwig/main.nf'
include { feature_counts ; matrix_entete_feature_counts ; matrix_colonne_feature_counts ; matrix_feature_counts } from './module/count/main.nf'
include { rsem_index ; rsem_calculate ; matrix_entete_rsem ; matrix_colonne_rsem_expected_count ; matrix_rsem_expected_count ; matrix_rsem ; matrix_colonne_rsem_TPM ; matrix_rsem as matrix_rsem_TPM ; matrix_rsem as matrix_rsem_FPKM ; matrix_colonne_rsem_FPKM } from './module/count/main.nf'
 



workflow {
	toto=Channel.fromPath(params.genome)
	/*toto.view()*/
	
	/*  ---Mapping--- */
	indexation_genome_STAR(Channel.fromPath(params.genome))
	indexation_genome_STAR.out.index.view()
	reads=Channel.fromFilePairs(params.reads)
	/*reads.view()*/
	anno=Channel.fromPath(params.annotation).collect()
	/*anno.view()*/
	mapping_STAR(indexation_genome_STAR.out.index.collect(),reads,anno)
	/*mapping_STAR.out.bam.view()*/
	
	/*  ---samtools--- */
	index_samtools(mapping_STAR.out.bam)
	
	/*  ---bigwig--- */
	bigwig(mapping_STAR.out.bam, index_samtools.out.index_bai_ch.collect())
	
	/*  ---feature_counts--- */
	feature_counts(mapping_STAR.out.bam,params.annotation)
	
	/*  ---entete de la matrice_feature_counts--- */
	matrix_entete_feature_counts(feature_counts.out.feature_counts_ch.first())
	matrix_entete_feature_counts=matrix_entete_feature_counts.out.entete_feature_counts_ch
	
	/*  --- colonnes de la matrice_feature_counts---  */
	matrix_colonne_feature_counts(feature_counts.out.feature_counts_ch)
	matrix_colonne_feature_counts=matrix_colonne_feature_counts.out.matrix_colonne_feature_counts_ch.toSortedList( { a, b -> a[0] <=> b[0]} ).flatten().filter { it =~ /tmp/ }.collect()
	matrix_colonne_feature_counts.view()
	
	/*  --- matrice feature_counts--- */
	matrix_feature_counts(matrix_colonne_feature_counts,matrix_entete_feature_counts)
	
	/*   ---rsem---  */
	rsem_index(Channel.fromPath(params.genome),Channel.fromPath(params.annotation),indexation_genome_STAR.out.index.collect())
    /*rsem_index.out.view()*/
    rsem_calculate(mapping_STAR.out.transcri_bam,rsem_index.out.rsem_index.collect())
    /*mapping_STAR.out.transcri_bam.view()*/
    /*rsem_calculate.out.rsem_results.view()*/
    
    /*   ---nom de fichier de sortie des matrices de rsem---  */
    name_matrix_rsem_exe=Channel.fromPath('matrice_rsem_expected_count.txt')
    name_matrix_rsem_T=Channel.fromPath('matrice_rsem_TPM.txt')
    name_matrix_rsem_F=Channel.fromPath('matrice_rsem_FPKM.txt')
    
    /*entete de la matrice*/
    matrix_entete_rsem(rsem_calculate.out.rsem_results)
    matrix_entete_rsem=matrix_entete_rsem.out.entete
    
    /*expected_count*/
    matrix_colonne_rsem_expected_count(rsem_calculate.out.rsem_results)
   	matrix_colonne_rsem_expected_count=matrix_colonne_rsem_expected_count.out.matrix_colonne_tmp.toSortedList( { a, b -> a[0] <=> b[0]} ).flatten().filter { it =~ /tmp/ }.collect()
	matrix_rsem(matrix_colonne_rsem_expected_count,matrix_entete_rsem,name_matrix_rsem_exe)
	
	/*TPM*/
	matrix_colonne_rsem_TPM(rsem_calculate.out.rsem_results)
	matrix_colonne_rsem_T=matrix_colonne_rsem_TPM.out.matrix_colonne_tmp.toSortedList( { a, b -> a[0] <=> b[0]} ).flatten().filter { it =~ /tmp/ }.collect()
	matrix_rsem_TPM(matrix_colonne_rsem_T,matrix_entete_rsem,name_matrix_rsem_T)
	
	/*FPKM*/
	matrix_colonne_rsem_FPKM(rsem_calculate.out.rsem_results)
	matrix_colonne_rsem_F=matrix_colonne_rsem_FPKM.out.matrix_colonne_tmp.toSortedList( { a, b -> a[0] <=> b[0]} ).flatten().filter { it =~ /tmp/ }.collect()
	matrix_rsem_FPKM(matrix_colonne_rsem_F,matrix_entete_rsem,name_matrix_rsem_F)
}


