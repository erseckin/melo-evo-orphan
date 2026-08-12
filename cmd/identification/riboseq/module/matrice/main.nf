nextflow.enable.dsl=2

process matrix_colonne {

	input:
		tuple val (pair_ID), path (file)
	
	output:
		tuple val(pair_ID), path ("${pair_ID}.tmp") , emit: col
	
	script:
		
	"""
		cut -f7 $file > $pair_ID".tmp"
	"""
}

/*
 * voir si cela peut faire partie d'un process revoir la synthase
 */

matrix_colonne_ch
  .toSortedList( { a, b -> a[0] <=> b[0]} )
  .flatten()
  .filter { it =~ /tmp/ }
  .set { matrice_colonne_triee_ch }



process matrix {
	publishDir "$params.outdir/results/matrice/", mode: 'copy'
	
	input:
		path (file) from matrice_colonne_triee_ch.collect()
		path (entete) from entete_matrice_ch
	
	output:
		file ("*.txt") into matrice_channel
	
	script:

	"""
		paste $entete $file > matrice_feature_counts.txt
		sed -i '1d' matrice_feature_counts.txt
		sed -i 's/Aligned.sortedByCoord.out.bam//g' matrice_feature_counts.txt
	"""
}


process matrix_entete_rsem {
	tag "matrix_entete_rsem"
	label "matrix_entete_rsem"
			
	input:
		tuple val (pair_ID), path (file) from rsem_calculate1_ch.first()
			
	output:
		path "entete_rsem.txt" into (entete_rsem1_ch , entete_rsem2_ch , entete_rsem3_ch)
	
	script:
		
	"""	
		cut -f1 $file > entete
		mv entete entete_rsem.txt
	"""
}

process matrix_colonne_rsem_expected_count {

	input:
		tuple val (pair_ID), path (file) from rsem_calculate2_ch
	
	output:
		tuple val(pair_ID), path ("${pair_ID}.tmp") into matrix_colonne_rsem1_ch
	
	script:
		
	"""
		cut -f5 $file > $pair_ID".tmp"
		sed -i '1d' $pair_ID".tmp"
		sed -i '1i$pair_ID' $pair_ID".tmp"
		
	"""
}

matrix_colonne_rsem1_ch
  .toSortedList( { a, b -> a[0] <=> b[0]} )
  .flatten()
  .filter { it =~ /tmp/ }
  .set { matrix_rsem1_colonne_triee_ch }



process matrix_rsem_expected_count {
	publishDir "$params.outdir/results/matrice/", mode: 'copy'
	
	input:
		path (file) from matrix_rsem1_colonne_triee_ch.collect()
		path (entete) from entete_rsem1_ch
	
	output:
		file ("*.txt") 
	
	script:

	"""
		paste $entete $file > matrice_rsem_expected_count.txt
	"""
}

process matrix_colonne_rsem_TPM {

	input:
		tuple val (pair_ID), path (file) from rsem_calculate3_ch
	
	output:
		tuple val(pair_ID), path ("${pair_ID}.tmp") into matrix_colonne_rsem_TPM_ch
	
	script:
		
	"""
		cut -f6 $file > $pair_ID".tmp"
		sed -i '1d' $pair_ID".tmp"
		sed -i '1i$pair_ID' $pair_ID".tmp"
		
	"""
}

matrix_colonne_rsem_TPM_ch
  .toSortedList( { a, b -> a[0] <=> b[0]} )
  .flatten()
  .filter { it =~ /tmp/ }
  .set { matrix_rsem_TPM_colonne_triee_ch }



process matrix_rsem_TPM {
	publishDir "$params.outdir/results/matrice/", mode: 'copy'
	
	input:
		path (file) from matrix_rsem_TPM_colonne_triee_ch.collect()
		path (entete) from entete_rsem2_ch
	
	output:
		file ("*.txt") 
	
	script:

	"""
		paste $entete $file > matrice_rsem_TPM.txt
	"""
}

process matrix_colonne_rsem_FPKM {

	input:
		tuple val (pair_ID), path (file) from rsem_calculate4_ch
	
	output:
		tuple val(pair_ID), path ("${pair_ID}.tmp") into matrix_colonne_rsem_FPKM_ch
	
	script:
		
	"""
		cut -f7 $file > $pair_ID".tmp"
		sed -i '1d' $pair_ID".tmp"
		sed -i '1i$pair_ID' $pair_ID".tmp"
		
	"""
}

matrix_colonne_rsem_FPKM_ch
  .toSortedList( { a, b -> a[0] <=> b[0]} )
  .flatten()
  .filter { it =~ /tmp/ }
  .set { matrix_rsem_FPKM_colonne_triee_ch }



process matrix_rsem_FPKM {
	publishDir "$params.outdir/results/matrice/", mode: 'copy'
	
	input:
		path (file) from matrix_rsem_FPKM_colonne_triee_ch.collect()
		path (entete) from entete_rsem3_ch
	
	output:
		file ("*.txt") 
	
	script:

	"""
		paste $entete $file > matrice_rsem_FPKM.txt
	"""
}
