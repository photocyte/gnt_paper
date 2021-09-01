nextflow.enable.dsl=2

process dump_and_bowtie {
  conda 'sra-tools=2.10.8 bowtie2=2.4.2 samtools=1.10 pv=1.6.6' 
  storeDir "results"
  errorStrategy "retry"
  maxRetries 50
  maxForks 1
  scratch 'ram-disk'
  cpus 6
  input: 
   tuple val(SRA), val(index_name), path(index_files)
  output:
   path "${SRA}-${index_name}.aligned_reads.wrg.bam"
  tag "${SRA}"
shell:
'''
    ##If you don't set this, by default it will cache huge files in your home directory. Gross!
    mkdir -p sra_cache
    vdb-config -s "/repository/user/main/public/root=$(pwd)/sra_cache"

    ## Pull down the reads
    fasterq-dump --threads !{task.cpus} --print-read-nr --split-files --temp ./ !{SRA}

    ## Delete source fastqs if too small.  Some sort of problem occured?
    find . -name '*.fastq' -type 'f' -size -160k -delete
        
    ## Concatenate the fastqs via a named pipe
    ## Use pipe viewer(pv) for some simple progress statistics.
    du -hs ./*.fastq
    FILESIZE=`du -scb ./*.fastq | grep total | cut -f 1`
    
    ## Run bowtie2
    cat *.fastq | pv -s "$FILESIZE" --interval 600 | bowtie2 --no-unal --threads !{task.cpus} --very-sensitive-local -x !{index_name} -q -U - --rg-id !{SRA} --rg "SM:!{SRA}" | samtools view -b -h -F 4 | samtools sort -@ 4 -T !{SRA}-!{index_name} -o !{SRA}-!{index_name}.aligned_reads.wrg.bam

    ## Cleanup the fastqs
    ls -1 | grep ".fastq" | xargs -P 1 -n 30 rm -f
'''
}

process bowtie2_index {
conda "bowtie2=2.4.2"
//module 'bowtie2:samtools'
input:
 path FNA
output:
 tuple val("${FNA}"), path("${FNA}*")
script:
"""
bowtie2-build ${FNA} ${FNA} 
"""
}

workflow {
 FASTAs = Channel.fromPath(params.fasta)
 bowtie2_index(FASTAs)

 SRAs = Channel.fromPath(params.sra_ids_file).splitText().map{it -> it.trim()}
 
 bowtie2_cmds = SRAs.combine(bowtie2_index.out)
 dump_and_bowtie(bowtie2_cmds)

}
