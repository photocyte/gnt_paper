#cat ./SraRunTable-10.txt | tail -n +2 | grep -P "" | ruby -rcsv -e 'CSV.foreach(ARGV.shift) {|row| puts row[0]}' /dev/stdin | sort | uniq > target_SRAs.txt

##If running into the SIGBUS issues, then should do this.
if [ ! -d .nextflow ]
then
THETMP=$(mktemp -d /dev/shm/TRF_XXXXXXXX)
echo $THETMP
ln -s ${THETMP} .nextflow
fi

nextflow run meta_T.nf --sra_ids_file target_SRAs.txt --fasta ./target_genes/targets.fa -resume -with-trace -with-report
