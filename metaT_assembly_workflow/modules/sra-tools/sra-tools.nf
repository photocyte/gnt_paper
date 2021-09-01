nextflow.enable.dsl=2

process fasterq_dump {
    conda 'sra-tools>2.10.8 pigz' 
    cpus 4 
    publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
    storeDir "results/${task.process}"
    maxForks 1
    scratch 'ram-disk'
    input:
      val 'sraid'
    output:
      path("${sraid}*.fastq.gz") //It gets complicated dealing with single vs paired libraries.
    tag "${sraid}"
    shell:
      '''
    ##See https://stackoverflow.com/questions/7642743/how-to-generate-random-numbers-in-the-busybox-shell
    ##Works on busybox, like the quay.io docker containers typically are
      date
      echo "Sleeping for random amount of time to avoid SRA rate limiting"
      sleep 17
      date
      echo "Done sleeping."
      ##If you don't set this, by default it will cache huge files in your home directory. Gross!
      mkdir -p sra_cache
      vdb-config -s "/repository/user/main/public/root=$(pwd)/sra_cache"

      echo "trying to prefetch the SRA accession..."
      prefetch !{sraid} --max-size 420000000000

      fasterq-dump --threads !{task.cpus} --print-read-nr --split-files --temp ./ !{sraid}
 
      for f in ./*.fastq
      do
      pigz --processes !{task.cpus} $f
      done
      wait
      '''
}

workflow {
    sraid = Channel.from(params.sraid)
    fasterq_dump(sraid)
}
