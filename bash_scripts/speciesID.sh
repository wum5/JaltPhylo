# move to target folder
#cd /N/dc2/projects/jaltomt/phylogeny/concatenate/besttree
cd /N/dc2/projects/jaltomt/introgression
#cd /N/dc2/projects/jaltomt/
FILENAME=new

for file in *$FILENAME*; do 
sed -i 's/JA0010/J. grandibaccata/g' $file
sed -i 's/JA0432/J. umbellata/g' $file
sed -i 's/JA0450/J. auriculata/g' $file
sed -i 's/JA0456/J. procumbens/g' $file
sed -i 's/JA0608/J. biflora/g' $file
sed -i 's/JA0694/J. darcyana/g' $file
sed -i 's/JA0701/J. repandidentata/g' $file
sed -i 's/JA0702/J. sinuosa/g' $file
sed -i 's/JA0711/J. calliantha/g' $file
sed -i 's/JA0719/J. dendroidea/g' $file
sed -i 's/JA0723/J. yungayensis/g' $file
sed -i 's/JA0726/J. aijana/g' $file
sed -i 's/JA0798/J. quipuscoae/g' $file
sed -i 's/JA0816/J. incahausina/g' $file
sed -i 's/Solyc/S. lycopersicum/g' $file
done

