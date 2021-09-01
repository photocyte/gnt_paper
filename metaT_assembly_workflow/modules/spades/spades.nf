nextflow.enable.dsl=2

process spades_meta_PE_nano {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
executor "pbs"
queue "home"
time '3d'
clusterOptions "-V -A moore-group -m n"
cpus 16
memory '725 GB'
input:
 tuple path(R1_reads),path(R2_reads)
 path(nanopore)
tag {"spades ${R1_reads}"}
shell:
'''
spades.py -t !{task.cpus} -m !{task.memory.toGiga()} --meta -k21 --pe-1 1 !{R1_reads} --pe-2 1 !{R2_reads} --nanopore !{nanopore} -o spades_out
'''
}

process spades_meta_PE {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
executor "pbs"
queue "home"
clusterOptions "-V -A moore-group -m n"
time '3d'
cpus 16
memory '725 GB'
input:
 tuple path(R1_reads),path(R2_reads)
tag {"spades ${R1_reads}"}
shell:
'''
spades.py -t !{task.cpus} -m !{task.memory.toGiga()} --meta -k21 --pe-1 1 !{R1_reads} --pe-2 1 !{R2_reads} --nanopore !{nanopore} -o spades_out
'''
}

process rnaspades_PE {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
cache 'lenient'
scratch 'ram-disk'
maxForks 1
executor "local"
queue "home"
clusterOptions "-V -A moore-group -m n"
time '3d'
cpus 16
memory '360 GB'
input:
 tuple path(R1_reads),path(R2_reads)
output:
 path("spades_out/transcripts.fasta"),emit:"scaffolds"
 path("spades_out/transcripts.paths"),emit:"paths"
 path("spades_out/assembly_graph_with_scaffolds.gfa"),emit:"gfa"
tag {"spades ${R1_reads}"}
shell:
'''
spades.py --rna -t !{task.cpus} -m !{task.memory.toGiga()} --pe-1 1 !{R1_reads} --pe-2 1 !{R2_reads} -o spades_out
'''
}

process rnaspades_PE_S {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
cache 'lenient'
scratch 'ram-disk'
maxForks 1
executor "local"
queue "home"
clusterOptions "-V -A moore-group -m n"
time '3d'
cpus 16
memory '360 GB'
input:
 tuple path(R1_reads),path(R2_reads),path(U_reads)
output:
 path("spades_out/transcripts.fasta"),emit:"scaffolds"
 path("spades_out/transcripts.paths"),emit:"paths"
 path("spades_out/assembly_graph_with_scaffolds.gfa"),emit:"gfa"
tag {"spades ${R1_reads}"}
shell:
'''
spades.py --rna -t !{task.cpus} -m !{task.memory.toGiga()} --pe-1 1 !{R1_reads} --pe-2 1 !{R2_reads} --pe-s 1 !{U_reads} -o spades_out
'''
}

process rnaspades_S {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
cache 'lenient'
scratch 'ram-disk'
maxForks 1
executor "local"
queue "home"
clusterOptions "-V -A moore-group -m n"
time '3d'
cpus 16
memory '360 GB'
input:
 path(U_reads)
output:
 path("spades_out/transcripts.fasta"),emit:"scaffolds"
 path("spades_out/transcripts.paths"),emit:"paths"
 path("spades_out/assembly_graph_with_scaffolds.gfa"),emit:"gfa"
tag {"spades ${R1_reads}"}
shell:
'''
spades.py --rna -t !{task.cpus} -m !{task.memory.toGiga()} --pe-s 1 !{U_reads} -o spades_out
'''
}

process rnaspades {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
cache 'lenient'
scratch 'ram-disk'
maxForks 2
executor "local"
queue "home"
clusterOptions "-V -A moore-group -m n"
time '3d'
cpus 8
memory '360 GB'
input:
 tuple val(LIBRARY_LAYOUT),path(reads)
 val(STRAND_SPECIFIC)
output:
 path("spades_out/transcripts.fasta"),emit:"scaffolds"
 path("spades_out/transcripts.paths"),emit:"paths"
 path("spades_out/assembly_graph_with_scaffolds.gfa"),emit:"gfa"
tag {"spades-${LIBRARY_LAYOUT}-${STRAND_SPECIFIC}SS-${reads[0]}"}
shell:
'''
if [ "!{LIBRARY_LAYOUT}" == "SINGLE" ]
then
READ_STR="--pe-s 1 !{reads[0]}"
elif [ "!{LIBRARY_LAYOUT}" == "PAIRED" ]
then
READ_STR="--pe-1 1 !{reads[0]} --pe-2 1 !{reads[1]}"
elif [ "!{LIBRARY_LAYOUT}" == "TRIPLE" ]
then
READ_STR="--pe-1 1 !{reads[0]} --pe-2 1 !{reads[1]} --pe-s 1 !{reads[2]}"
else
echo "Unknown library layout."
exit
fi

if [ "!{STRAND_SPECIFIC}" == "yes" ]
then
SS_STR="--ss fr"
elif [ "!{STRAND_SPECIFIC}" == "YES" ]
then
SS_STR="--ss fr"
else
SS_STR=""
fi

spades.py --rna -t !{task.cpus} -m !{task.memory.toGiga()} ${READ_STR} ${SS_STR} -o spades_out
'''
}

process spades_small {
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
//scratch 'ram-disk'
maxForks 1
executor "local"
cpus 4
memory '50 GB'
input:
 tuple path(R1_reads),path(R2_reads),path(U_reads)
output:
 path("spades_out/scaffolds.fasta"),emit:"scaffolds"
 path("spades_out/scaffolds.paths"),emit:"paths"
 path("assembly_graph_with_scaffolds.gfa"),emit:"gfa"
tag {"spades ${R1_reads}"}
shell:
'''
spades.py --only-assembler --cov-cutoff off -k 21 -t !{task.cpus} -m !{task.memory.toGiga()} --pe-1 1 !{R1_reads} --pe-2 1 !{R2_reads} --pe-s 1 !{U_reads} -o spades_out
'''
}


workflow {

forward = Channel.fromPath(params.forward)
reverse = Channel.fromPath(params.reverse)
pairedReads = forward.combine(reverse)
rnaspades_PE(pairedReads)

}
