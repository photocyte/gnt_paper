nextflow.enable.dsl=2

process prokka_do {
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    input:
      path fasta
    output:
      path("${fasta}_PROKKA")
    echo true
    shell:
      '''
      prokka !{fasta} --compliant --centre MooreLab --outdir !{fasta}_PROKKA --locustag !{fasta} --prefix !{fasta}_PROKKA
      '''
}

process prokka_do_with_pep {
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    input:
      path fasta
      path peps
    output:
      path("${fasta}_PROKKA")
    echo true
    shell:
      '''
      prokka !{fasta} --proteins !{peps} --compliant --centre MooreLab --outdir !{fasta}_PROKKA --locustag !{fasta} --prefix !{fasta}_PROKKA
      '''
}


workflow {
 prokka_do(Channel.fromPath(params.fasta))
}
