nextflow.enable.dsl = 2

process EXPAND_AMBIGUOUS_SEQUENCES {
	tag "$meta.uid"
	conda "python=3.7"

	input:
	tuple val(meta), val(seq)

	output:
	tuple val(meta), path("*.acgt.fasta")

	script:
	"""
	#!/usr/bin/env python

	from itertools import product

	iupac_codes = dict(
	    A="A", C="C", G="G", T="T",
	    R="AG", Y="CT", S="GC", W="AT", K="GT", M="AC",
		B="CGT", D="AGT", H="ACT", V="ACG", N="ACGT"
    )

	choices = [iupac_codes[nucl] for nucl in "$seq"]

	with open("${meta.uid}.acgt.fasta", "w") as writer:
	    for i, comb in enumerate(product(*choices)):
	        seq = "".join(comb)
	        writer.write(f">${meta.uid}.{i}\\n{seq}\\n")
	"""
}

process MAKE_BLAST_DB {
	container "quay://biocontainers/blast=2.12.0--hf3cf87c_4"
	
	input:
	path fasta

	output:
	path "*"

	script:
	"""
	makeblastdb -dbtype nucl -in $fasta -out ref
	"""
}

process BLASTN {
    tag "$meta.uid"
	cpus 5
	container "quay://biocontainers/blast=2.12.0--hf3cf87c_4"
	publishDir "$params.outdir/blastn", mode: "copy"

    input:
	tuple val(meta), path(query)
	path db

    output:
	tuple val(meta), path("*.tsv")

    script:
	extra_fields = meta.collect{key, value -> key;}.join(" ")
	extra_values = meta.collect{key, value -> value;}.join('\t')
    """
	cat <(
	  echo $params.blast_fields $extra_fields | tr " " "\t") <(
	  blastn -num_threads $task.cpus $params.blast_args \\
		-perc_identity $params.min_pident \\
		-query $query -db ref \\
		-outfmt "6 $params.blast_fields" \\
	-out - | awk '{print \$0"\t$extra_values"}' \\
    ) > blast.${meta.uid}.tsv
    """
}

process SUMMARIZE_BLAST {
	tag "summary"
	conda "pandas=1.4"
	publishDir "$params.outdir", mode: "copy"

	input:
	path blast_file

	output:
	path "summary.csv"

	script:
	"""
	#!/usr/bin/env python

	import pandas as pd

	blast = pd.read_csv("$blast_file", sep="\\t")
	counts = blast.groupby(["sacc", "primer"], as_index=False).agg(dict(uid="nunique", gene="first"))
	counts = counts.groupby(["sacc", "gene"]).uid.agg(lambda cts: sum(ct>1 for ct in cts))
	counts = counts.unstack(fill_value=0)
	counts["nb_genes"] = (counts > 0).sum(axis=1)
	counts = counts[counts.nb_genes > 0].sort_values(by="nb_genes", ascending=False)

	counts.to_csv("summary.csv")
	"""
}

workflow AMBIGUOUS_BLASTN{
	take:
	primers // one fasta file with all primers. Sequence title needs to be "{primer_name} {fwd,rev}"
	contigs 
	
	main:
	db = MAKE_BLAST_DB(contigs)
	primers = primers.splitFasta( record: [header: true, seqString: true ] )
		.map{
			def (primer_name, orient, name) = it.header.tokenize();
			def (gene, count) = primer_name.contains("_") ? primer_name.tokenize("_") : [primer_name, "0"];
			[[uid: "${primer_name}_${orient}", primer: primer_name, orient: orient, gene: gene, raw_name: name], it.seqString]
		}
	primers_acgt = EXPAND_AMBIGUOUS_SEQUENCES(primers)
	blast = BLASTN(primers_acgt, db)
	blast_all = blast.map{it[1]}.collectFile(
		name: "blast_results.tsv", storeDir: params.outdir,
		skip: 1, keepHeader: true
	)
	SUMMARIZE_BLAST(blast_all)
}

workflow {
	primers = Channel.fromPath(params.primers)
	contigs = file(params.contigs, checkIfExists: true)

	AMBIGUOUS_BLASTN(primers, contigs)
}
