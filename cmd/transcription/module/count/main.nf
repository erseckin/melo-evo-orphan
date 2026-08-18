nextflow.enable.dsl=2

process feature_counts {
	label "feature_counts"
	tag "${feature_counts}"

	publishDir "$params.outdir/results/feature_counts/", pattern:"*.txt", mode: 'copy', saveAS: { filename -> "${pair_ID}.txt"}
	publishDir "$params.outdir/results/feature_counts_log/", pattern:"*.txt.summary", mode: 'copy', saveAS: { filename -> "${pair_ID}.txt.summary"}
	
	input:
		tuple val(pair_ID), path (file)
		val annotation
	
	output:
		tuple val (pair_ID), path ("${pair_ID}_count.txt"), emit : feature_counts_ch
		file "*"
	
	script:	
	
	"""	
		featureCounts -a $annotation -o $pair_ID"_count.txt" $params.option_feature_counts $file
		
	"""
}

process matrix_entete_feature_counts {
	tag "matrix_entete_feature_counts"
	label "matrix_entete_feature_counts"
			
	input:
		tuple val (pair_ID), path (file)
			
	output:
		path "entete.feature.counts.txt", emit : entete_feature_counts_ch
	
	script:
		
	"""	
		cut -f1 $file > entete
		mv entete entete.feature.counts.txt
	"""
}

process matrix_colonne_feature_counts {

	input:
		tuple val (pair_ID), path (file)
	
	output:
		tuple val(pair_ID), path ("${pair_ID}.tmp"), emit : matrix_colonne_feature_counts_ch
	
	script:
		
	"""
		cut -f7 $file > $pair_ID".tmp"
	"""
}

process matrix_feature_counts {
	publishDir "$params.outdir/results/matrice/", mode: 'copy'
	
	input:
		path (file)
		path (entete)
	
	output:
		path ("matrice_feature_counts.txt"), emit : matrice_feature_counts_ch
	
	script:

	"""
		paste $entete $file > matrice_feature_counts.txt
		sed -i '1d' matrice_feature_counts.txt
		sed -i 's/Aligned.sortedByCoord.out.bam//g' matrice_feature_counts.txt
	"""
}



process rsem_index {
	label "rsem_index"
	tag "${rsem_index}"
	
	input:
		path(genome)
		path(annotation)
		path(index_star)
			
	output:
	
		path ("genome_rsem*")  ,emit : rsem_index
			
	script:
	"""	
		rsem-prepare-reference --gff3 $annotation --star $genome genome_rsem
		
	"""
} 

process rsem_calculate {
	label "rsem_calculate"
	tag "${rsem_calculate}"

	publishDir "$params.outdir/results/rsem_count/", mode: 'copy', pattern:"*.genes.results"
	publishDir "$params.outdir/results/rsem_count/", mode: 'copy', pattern:"*.isoforms.results"
			
	input:


		tuple val(pair_ID), path (file)
		path(rsem_index)
	
	output:
		tuple val  (pair_ID), path ("*.genes.results"), emit : rsem_results
		file "*.isoforms.results"
		
	script:
	"""	
		rsem-calculate-expression $params.option_rsem $file genome_rsem $pair_ID
		
	"""
} 
 

process matrix_entete_rsem {
	tag "${matrix_entete_rsem}"
	label "matrix_entete_rsem"
			
	input:
		tuple val (pair_ID), path (file)
			
	output:
		path "entete_rsem.txt" , emit : entete
	
	script:
		
	"""	
		cut -f1 $file > entete
		mv entete entete_rsem.txt
	"""
}

process matrix_colonne_rsem_expected_count {

	input:
		tuple val (pair_ID), path (file)
	
	output:
		tuple val(pair_ID), path ("${pair_ID}.tmp"),  emit: matrix_colonne_tmp
	
	script:
		
	"""
		cut -f5 $file > $pair_ID".tmp"
		sed -i '1d' $pair_ID".tmp"
		sed -i '1i$pair_ID' $pair_ID".tmp"
		
	"""
}



process matrix_rsem_expected_count {

	publishDir "$params.outdir/results/matrice/", mode: 'copy'
	
	input:
		
		path (file)
		path (entete)
	output:
		file ("*.txt"), emit : matrix_rsem_expected_count_ch
	
	script:

	"""
		paste $entete $file > matrice_rsem_expected_count.txt
	"""
}





process matrix_rsem {

	publishDir "$params.outdir/results/matrice/", mode: 'copy'
	
	input:
		
		path (file)
		path (entete)
		path (name)
		
	output:
		path ("$name"), emit : matrice_ch
	
	script:

	"""
		paste $entete $file > $name
	"""
}

process matrix_colonne_rsem_TPM {

	input:
		tuple val (pair_ID), path (file)
	
	output:
		tuple val(pair_ID), path ("${pair_ID}.tmp"), emit : matrix_colonne_tmp
	
	script:
		
	"""
		cut -f6 $file > $pair_ID".tmp"
		sed -i '1d' $pair_ID".tmp"
		sed -i '1i$pair_ID' $pair_ID".tmp"
		
	"""
}


process matrix_colonne_rsem_FPKM {

	input:
		tuple val (pair_ID), path (file)
	
	output:
		tuple val(pair_ID), path ("${pair_ID}.tmp") , emit : matrix_colonne_tmp
	
	script:
		
	"""
		cut -f7 $file > $pair_ID".tmp"
		sed -i '1d' $pair_ID".tmp"
		sed -i '1i$pair_ID' $pair_ID".tmp"
		
	"""
}
