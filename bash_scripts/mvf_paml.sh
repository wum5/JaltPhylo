#!/bin/bash

#PBS -N paml_jaltBrch
#PBS -l nodes=1:ppn=1,walltime=48:00:00
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module unload python
module load python/3.3.0
module load paml

CD=/N/dc2/projects/jaltomt/Phylogenomics/paml
OD=/N/dc2/projects/jaltomt/Phylogenomics/paml/paml_new
cd /N/dc2/projects/jaltomt/Softwares/mvftools-dev-master



## branch-site test on ancestral branch of Clade including red, orange and green lineages
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_noCap_codon --out $OD/Clade_out --allelegroups RED:JA0450,JA0432,JA0608,JA0702,JA0719,JA0723,JA0726,JA0816,JA0711,JA0798 OTHER:JA0456,JA0701,JA0694,Solyc --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_Clade --pamltmp PAMLtemp_Clade --startcontig 0 --endcontig 6772 --target JA0450 JA0432 JA0608 JA0702 JA0719 JA0723 JA0726 JA0816 JA0711 JA0798 --targetspec 10 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree

## branch-site test on ancestral branch of orange lineages
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_noCap_codon --out $OD/Orange_out --allelegroups RED:JA0432,JA0608,JA0702,JA0719,JA0723,JA0726,JA0816 OTHER:JA0456,JA0701,JA0694,JA0450,JA0711,JA0798,Solyc --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_Orange --pamltmp PAMLtemp_Orange --startcontig 0 --endcontig 6772 --target JA0432 JA0608 JA0702 JA0719 JA0723 JA0726 JA0816 --targetspec 7 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree

## branch-site test on ancestral branch of black lineages
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_noCap_codon --out $OD/Black_out --allelegroups RED:JA0456,JA0701,JA0694 OTHER:JA0450,JA0432,JA0608,JA0702,JA0719,JA0723,JA0726,JA0816,JA0711,JA0798,Solyc --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_Black --pamltmp PAMLtemp_Black --startcontig 0 --endcontig 6772 --target JA0456 JA0701 JA0694 --targetspec 3 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree

## branch-site test on ancestral branch of green lineages
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_noCap_codon --out $OD/Green_out --allelegroups RED:JA0711,JA0798 OTHER:JA0456,JA0701,JA0694,JA0450,JA0432,JA0608,JA0702,JA0719,JA0723,JA0726,JA0816,Solyc --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_Green --pamltmp PAMLtemp_Green --startcontig 0 --endcontig 6772 --target JA0711 JA0798 --targetspec 2 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree

## branch-site test on ancestral branch of red lineage
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_noCap_codon --out $OD/Red_out --allelegroups RED:JA0450 OTHER:JA0456,JA0701,JA0694,JA0432,JA0608,JA0702,JA0719,JA0723,JA0726,JA0816,JA0711,JA0798,Solyc --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_Red --pamltmp PAMLtemp_Red --startcontig 0 --endcontig 6772 --target JA0450 --targetspec 1 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree

## clade test on Clade including orange and green lineages
#python3.3 mvf_analyze_codon2.py GroupUniqueAlleleWindow --mvf $CD/Jalt_noCap_codon --out $OD/CladeOG_out --allelegroups RED:JA0432,JA0608,JA0702,JA0719,JA0723,JA0726,JA0816,JA0711,JA0798 OTHER:JA0456,JA0701,JA0694,JA0450,Solyc --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_CladeOG --pamltmp PAMLtemp_CladeOG --startcontig 0 --endcontig 6772 --target JA0432 JA0608 JA0702 JA0719 JA0723 JA0726 JA0816 JA0711 JA0798 --targetspec 9 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree


## branch-site test on Jaltomata ancestral branch
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_wtCap_codon --out $OD/JaltBrch_out --allelegroups RED:JA0456,JA0701,JA0694,JA0450,JA0432,JA0608,JA0702,JA0719,JA0723,JA0726,JA0816,JA0711,JA0798 OTHER:Solyc,Capana --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 SOL:Solyc CAP:Capana --branchlrt $OD/Geneoutput_JaltBrch --pamltmp PAMLtemp_JaltBrch --startcontig 0 --endcontig 4250 --target JA0456 JA0701 JA0694 JA0450 JA0432 JA0608 JA0702 JA0719 JA0723 JA0726 JA0816 JA0711 JA0798 --targetspec 13 --raxmlpath raxmlHPC --outgroup Capana --allsampletree

## clade test on Jaltomata ancestral branch
#python3.3 mvf_analyze_codon2.py GroupUniqueAlleleWindow --mvf $CD/Jalt_wtCap_codon --out $OD/JaltClade_out --allelegroups RED:JA0456,JA0701,JA0694,JA0450,JA0432,JA0608,JA0702,JA0719,JA0723,JA0726,JA0816,JA0711,JA0798 OTHER:Solyc,Capana --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 SOL:Solyc CAP:Capana --branchlrt $OD/Geneoutput_JaltClade --pamltmp PAMLtemp_JaltClade --startcontig 0 --endcontig 4250 --target JA0456 JA0701 JA0694 JA0450 JA0432 JA0608 JA0702 JA0719 JA0723 JA0726 JA0816 JA0711 JA0798 --targetspec 13 --raxmlpath raxmlHPC --outgroup Capana --allsampletree



## PhyloGWAS on the derived floral traits in Jaltomata (nectar)
python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_noSolyc_codon --out $OD/Jalt_nectar --allelegroups RED:JA0432,JA0608,JA0719,JA0726,JA0816,JA0711,JA0798 OTHER:JA0456,JA0701,JA0694,JA0450,JA0723,JA0702 --windowsize -1 --uselabels --speciesgroups PRO:JA0456 REP:JA0701 DAR:JA0694 AUR:JA0450 UMB:JA0432 BIF:JA0608 SIN:JA0702 DEN:JA0719 YUN:JA0723 AIJ:JA0726 INC:JA0816 CAL:JA0711 QUI:JA0798 --branchlrt $OD/Geneoutput_nectar --pamltmp PAMLtemp_nectar --startcontig 0 --endcontig 0 --target JA0432 JA0608 JA0719 JA0726 JA0816 JA0711 JA0798 --targetspec 8 --raxmlpath raxmlHPC --allsampletree

