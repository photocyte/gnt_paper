nextflow.enable.dsl=2

process makeblastdb_nucl {
input:
 path target
output:
 tuple path("${target}"), path("${target}*")
shell:
'''
makeblastdb -in !{target} -dbtype nucl
'''
}

process blastn {
cpus 4
    input:
     path queries
     path dbs
    output:
     path "blast.out"
    shell:
    '''
   blastn -task blastn -num_threads !{task.cpus} -db !{target[0]} -query !{queries} -evalue 1e-20 -out blast.out -outfmt 6
'''
}

process tblastn {
    cpus 4
    input:
     path target
     path queries
    output:
     path "blast.out",emit:"blast_out"
     path target,emit:"target"
    shell:
    '''
    makeblastdb -version
    tblastn -version 
    makeblastdb -in !{target} -dbtype nucl
    tblastn -num_threads !{task.cpus} -db !{target} -query !{queries} -evalue 1e-20 -out blast.out -outfmt 6

    if [ ! -s blast.out ]
    then
    rm -f blast.out
    fi
    '''
}

workflow {
 query = Channel.fromPath(query)
 target = Channel.fromPath(target)
 tblastn(query,target)
}
