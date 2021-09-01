##Assuming the files are unzipped to the contained directory structure
#unzip 

##Fetching only those read alignments with >90% identity
find . -name "*.m8" -type f -exec bash -c '[[ $(wc -l < "$1") -gt 2 ]] && echo "$1"' _ '{}' \; | xargs -n 1 awk '$3>90' | gzip > 90.m8.gz

##Fetching only those read alignments with >40% identity
find . -name "*.m8" -type f -exec bash -c '[[ $(wc -l < "$1") -gt 2 ]] && echo "$1"' _ '{}' \; | xargs -n 1 awk '$3>40' | gzip > 40.m8.gz


##Fetching those with >90% identity, and e-value <1e-10
#zcat 90.m8.gz | awk '$10<1e-10' | cut -f 1 | cut -f 1 -d "." | sort | uniq -c | sort -nr

##Just fetching those with >90% identity
#zcat 90.m8.gz | cut -f 1 | cut -f 1 -d "." | sort | uniq -c | sort -nr

##New way to do it, that doesn't require having all the files unzipped!
unzip -p results.zip "*.m8" | awk '$3>40' | gzip > 40.m8.gz 
