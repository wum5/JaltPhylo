#columns for blast output:
#1-qseqid 2-qlen 3-sseqid 4-slen 5-frames 6-pident 7-nident 8-length 
#9-mismatch 10-gapopen 11-qstart 12-qend 13-sstart 14-send 15-evalue 16-bitscore

cat *.blastn >all.rawblastn
awk -F '\t' '{if($1!=$3 && $6>0.8 && ($12-$11)/$2>0.7) print $1,$3,$15}' all.rawblastn >all.evalue
sed -i 's/  / /g' all.evalue
awk '{if($3==0.0)print $1,$2,180; else print $1,$2,-log($3)/log(10)}' all.evalue >all.evalue-log