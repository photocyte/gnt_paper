nextflow.enable.dsl=2

include { fasterq_dump } from './modules/sra-tools/sra-tools.nf'
include { rnaspades_PE; rnaspades_S; rnaspades } from './modules/spades/spades.nf'
include { prokka_do; prokka_do_with_pep } from './modules/prokka/prokka.nf'
include { tblastn } from './modules/blast/blast.nf'

process determine_library_layout {
input:
    path fastqs
output:
    tuple(env(LIBRARY_LAYOUT), path(fastqs))
shell:
'''
if [ `ls -1 ./*.fastq.gz | wc -l` -eq 1 ]
then
LIBRARY_LAYOUT="SINGLE"
elif [ `ls -1 ./*.fastq.gz | wc -l` -eq 2 ]
then
LIBRARY_LAYOUT="PAIRED"
elif [ `ls -1 ./*.fastq.gz | wc -l` -eq 3 ]
then
LIBRARY_LAYOUT="TRIPLE"
else
LIBRARY_LAYOUT="UNKNOWN"
fi
'''
}

process fasta_add_prefix {
 conda "seqkit"
 input:
  val prefix
  path fasta
 output:
  path "${prefix}_spades.fa"
 shell:
 '''
 seqkit replace -p "^" -r "!{prefix}_" !{fasta} > !{prefix}_spades.fa
 '''
}

process seqkit_blast_filter {
    conda "seqkit"
    input:
     path(blast_results)
     path(fasta)
    output:
     path "results/${fasta}"
shell:
    '''
    mkdir results
    cat !{blast_results} | awk '$3 > 75' | cut -f 2 | sort | uniq > filter_by.txt
    seqkit grep -f filter_by.txt !{fasta} > results/!{fasta}
    '''
}

workflow rnaAssemble {
take: 
 sraId
main:
 fasterq_dump(sraId)
 determine_library_layout(fasterq_dump.out)

 rnaspades(determine_library_layout.out,"no") //"no" = strand specific or not
 
 fasta_add_prefix(sraId,rnaspades.out.scaffolds)
 emit:
  fasta_add_prefix.out
}

workflow filterByBlast {
take:
 target_pep_tuple
main:
 target = target_pep_tuple.map{it[0]}
 pep = target_pep_tuple.map{it[1]}
    
 tblastn(target,pep)
    
 //out_target_tuple = tblastn.out.combine(target.first())
 //seqkit_tblastn = out_target_tuple.map{it[0]}
 //seqkit_target = out_target_tuple.map{it[1]}
    
 seqkit_blast_filter(tblastn.out.blast_out,tblastn.out.target)
emit:
 seqkit_blast_filter.out   
}

workflow {
target_list = params.targets.split(";").toList() // Groovy split, plus have to convert to "java.util.Collection"
targets = Channel.fromList( target_list )
peps = Channel.fromPath(params.peps)
rnaAssemble(targets)
filterByBlast(rnaAssemble.out.combine(peps))
prokka_tup = filterByBlast.out.combine(peps)
prokka_do_with_pep(prokka_tup.map{it[0]},prokka_tup.map{it[1]})
}
