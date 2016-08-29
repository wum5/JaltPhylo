#!/bin/bash

#PBS -N Jal-paml-Green_sp
#PBS -l nodes=1:ppn=1,walltime=32:00:00
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module unload python
module load python/3.3.0
module load paml

CD=/N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog/MVF_PAML/withCap
OD=/N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog/MVF_PAML/noCap
cd /N/dc2/projects/jaltomt/software/mvftools-dev-master


#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_ortho_codon --out $CD/JaltBrch_out --allelegroups RED:JA0711,JA0798,JA0723,JA0432,JA0702,JA0816,JA0719,JA0694,JA0014,JA0701 OTHER:Solyc,Capana --windowsize -1 --uselabels --speciesgroups UMB:JA0432 DAR:JA0694 REP:JA0701,JA0014 SIN:JA0702 CAL:JA0711 DEN:JA0719 YUN:JA0723 INC:JA0816 QUI:JA0798 SOL:Solyc CAP:Capana --branchlrt $CD/Geneoutput_JaltBrch --pamltmp PAMLtemp_JaltBrch --startcontig 0 --endcontig 3906 --target JA0711 JA0798 JA0723 JA0432 JA0702 JA0816 JA0719 JA0694 JA0014 JA0701 --targetspec 10 --raxmlpath raxmlHPC --outgroup Capana --allsampletree
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_ortho_codon --out $CD/SolnBrch_out --allelegroups RED:Solyc OTHER:JA0711,JA0798,JA0723,JA0432,JA0702,JA0816,JA0719,JA0694,JA0014,JA0701,Capana --windowsize -1 --uselabels --speciesgroups UMB:JA0432 DAR:JA0694 REP:JA0701,JA0014 SIN:JA0702 CAL:JA0711 DEN:JA0719 YUN:JA0723 INC:JA0816 QUI:JA0798 SOL:Solyc CAP:Capana --branchlrt $CD/Geneoutput_SolnBrch --pamltmp PAMLtemp_SolnBrch --startcontig 0 --endcontig 3906 --target Solyc --targetspec 1 --raxmlpath raxmlHPC --outgroup Capana --allsampletree

#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_ortho_codon --out $CD/Clade2_out --allelegroups RED:JA0711,JA0798,JA0723,JA0432,JA0702,JA0816,JA0719 OTHER:JA0694,JA0014,JA0701,Solyc,Capana --windowsize -1 --uselabels --speciesgroups UMB:JA0432 DAR:JA0694 REP:JA0701,JA0014 SIN:JA0702 CAL:JA0711 DEN:JA0719 YUN:JA0723 INC:JA0816 QUI:JA0798 SOL:Solyc CAP:Capana --branchlrt $CD/Geneoutput_Clade2 --pamltmp PAMLtemp_Clade2_sp --startcontig 0 --endcontig 3906 --target JA0711 JA0798 JA0723 JA0432 JA0702 JA0816 JA0719 --targetspec 7 --raxmlpath raxmlHPC --outgroup Capana --allsampletree
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_ortho_codon --out $CD/blackBrch_out --allelegroups RED:JA0694,JA0014,JA0701 OTHER:JA0711,JA0798,JA0723,JA0432,JA0702,JA0816,JA0719,Solyc,Capana --windowsize -1 --uselabels --speciesgroups UMB:JA0432 DAR:JA0694 REP:JA0701,JA0014 SIN:JA0702 CAL:JA0711 DEN:JA0719 YUN:JA0723 INC:JA0816 QUI:JA0798 SOL:Solyc CAP:Capana --branchlrt $CD/Geneoutput_blackBrch --pamltmp PAMLtemp_blackBrch_sp --startcontig 0 --endcontig 3906 --target JA0694 JA0014 JA0701 --targetspec 3 --raxmlpath raxmlHPC --outgroup Capana --allsampletree
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $CD/Jalt_ortho_codon --out $CD/greenBrch_out --allelegroups RED:JA0711,JA0798 OTHER:JA0694,JA0014,JA0701,JA0723,JA0432,JA0702,JA0816,JA0719,Solyc,Capana --windowsize -1 --uselabels --speciesgroups UMB:JA0432 DAR:JA0694 REP:JA0701,JA0014 SIN:JA0702 CAL:JA0711 DEN:JA0719 YUN:JA0723 INC:JA0816 QUI:JA0798 SOL:Solyc CAP:Capana --branchlrt $CD/Geneoutput_greenBrch --pamltmp PAMLtemp_greenBrch_sp --startcontig 0 --endcontig 3906 --target JA0711 JA0798 --targetspec 2 --raxmlpath raxmlHPC --outgroup Capana --allsampletree

#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $OD/Jalt_ortho_codon --out $OD/Clade2_out --allelegroups RED:JA0711,JA0798,JA0723,JA0432,JA0702,JA0816,JA0719 OTHER:JA0694,JA0014,JA0701,Solyc --windowsize -1 --uselabels --speciesgroups UMB:JA0432 DAR:JA0694 REP:JA0701,JA0014 SIN:JA0702 CAL:JA0711 DEN:JA0719 YUN:JA0723 INC:JA0816 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_Clade2 --pamltmp PAMLtemp_Clade2 --startcontig 0 --endcontig 6223 --target JA0711 JA0798 JA0723 JA0432 JA0702 JA0816 JA0719 --targetspec 7 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $OD/Jalt_ortho_codon --out $OD/blackBrch_out --allelegroups RED:JA0694,JA0014,JA0701 OTHER:JA0711,JA0798,JA0723,JA0432,JA0702,JA0816,JA0719,Solyc --windowsize -1 --uselabels --speciesgroups UMB:JA0432 DAR:JA0694 REP:JA0701,JA0014 SIN:JA0702 CAL:JA0711 DEN:JA0719 YUN:JA0723 INC:JA0816 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_blackBrch --pamltmp PAMLtemp_blackBrch --startcontig 0 --endcontig 6223 --target JA0694 JA0014 JA0701 --targetspec 3 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree
#python3.3 mvf_analyze_codon.py GroupUniqueAlleleWindow --mvf $OD/Jalt_ortho_codon --out $OD/greenBrch_out --allelegroups RED:JA0711,JA0798 OTHER:JA0694,JA0014,JA0701,JA0723,JA0432,JA0702,JA0816,JA0719,Solyc --windowsize -1 --uselabels --speciesgroups UMB:JA0432 DAR:JA0694 REP:JA0701,JA0014 SIN:JA0702 CAL:JA0711 DEN:JA0719 YUN:JA0723 INC:JA0816 QUI:JA0798 SOL:Solyc --branchlrt $OD/Geneoutput_greenBrch --pamltmp PAMLtemp_greenBrch --startcontig 0 --endcontig 6223 --target JA0711 JA0798 --targetspec 2 --raxmlpath raxmlHPC --outgroup Solyc --allsampletree
