# move to target folder
cd /N/dc2/projects/jaltomt/Phylogenomics/phylogeny/gene_tre
FILENAME='double_origin.txt'

for file in *$FILENAME*; do 
sed -i 's/Solyc/S.lycopersicum/g' $file
sed -i 's/JA0701/J.repandidentata/g' $file
sed -i 's/JA0456/J.procumbens/g' $file
sed -i 's/JA0694/J.darcyana/g' $file
sed -i 's/JA0450/J.auriculata/g' $file
sed -i 's/JA0798/J.quipuscoae/g' $file
sed -i 's/JA0711/J.calliantha/g' $file
sed -i 's/JA0723/J.yungayensis/g' $file
sed -i 's/JA0608/J.biflora/g' $file
sed -i 's/JA0702/J.sinuosa/g' $file
sed -i 's/JA0726/J.aijana/g' $file
sed -i 's/JA0432/J.umbellata/g' $file
sed -i 's/JA0010/J.grandibaccata/g' $file
sed -i 's/JA0719/J.dendroidea/g' $file
sed -i 's/JA0816/J.incahuasina/g' $file
done

#for file in *$FILENAME*; do 
#sed -i 's/S.lycopersicum/Solycp/g' $file
#sed -i 's/J.repandidentata/JA0701/g' $file
#sed -i 's/J.procumbens/JA0456/g' $file
#sed -i 's/J.darcyana/JA0694/g' $file
#sed -i 's/J.auriculata/JA0450/g' $file
#sed -i 's/J.quipuscoae/JA0798/g' $file
#sed -i 's/J.calliantha/JA0711/g' $file
#sed -i 's/J.yungayensis/JA0723/g' $file
#sed -i 's/J.biflora/JA0608/g' $file
#sed -i 's/J.sinuosa/JA0702/g' $file
#sed -i 's/J.aijana/JA0726/g' $file
#sed -i 's/J.umbellata/JA0432/g' $file
#sed -i 's/J.grandibaccata/JA0010/g' $file
#sed -i 's/J.dendroidea/JA0719/g' $file
#sed -i 's/J.incahuasina/JA0816/g' $file
#done

