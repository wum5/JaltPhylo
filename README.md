# Jaltomata Phylogenomics

## Table of Contents
* [Overview](#overview)
* [Contributors](#contributors)
* [Raw Data Processing](#raw-data-processing)
* [Homolog Inference](#homolog-inference)
* [Ortholog Inference](#ortholog-inference)
* [Alignment construction](#alignment-construction)
* [Phylogeny Construction](#phylogeny-construction)
* [Introgression Analysis](#introgression-analysis)
* [Ancestral Segregating Allele Analysis](#ancestral-segregating-allele-analysis)
* [Adaptive Evolution Analysis](#adaptive-evolution-analysis)

## Overview
* Raw scripts/Pipeline for the "Jaltomato Phylogenomics" Project.
* Some scripts were written by Ya Yang for her study (https://bitbucket.org/yangya/phylogenomic_dataset_construction).
* Scripts associated with MVF-format data processsing can be found in mvftools (https://github.com/jbpease/mvftools).
* Still in updating!

## Contributors 
* Meng Wu
* https://github.com/wum5/JaltPhylo

## Raw Data Processing
##### Trim low-quality reads and the first 15-bp of reads due to non-random hexamer primers
```
qsub trim.sh
qsub clip5end.sh
qsub FastaQC.sh
```
##### Build transcript assembly and predict CDS using Transdecoder
```
qsub trinity.sh
qsub transdecoder.sh
```
##### Rename CDS files and reduce redundancy
```
for file in *_dir; do cp $file/longest_orfs.cds outDIR/$file'.cds'; done
python fix_names_from_transdecoder.py <DIR> <DIR>
cat *_NR.fa *_RP.fa > *_cds.fa
qsub cd-hit-est.sh
```

## Homolog Inference
##### Make all-by-all blast and Infer putative homolog groups using similarity
```
qsub blastn.sh
cat *blastn > all.rawblast
python blast_to_mcl.py all.rawblast <hit_fraction_cutoff>
mcl all.rawblast.hit-frac0.4.minusLogEvalue --abc -te 5 -tf 'gq(10)' -I 2.5 -o hit-frac0.4_I2.5_e10
python write_fasta_files_from_mcl.py <fasta files> <mcl_outfile> <minimal_ingroup_taxa> <outDIR>
```
##### make initial alignments and then cut long internal branch
```
qsub mafft.sh
qsub phyutility.sh
qsub fasttree.sh
python cut_long_branches_iter.py <inDIR> <outDIR>
```
##### refine the final clusters
```
qsub mafft.sh
qsub phyutility.sh
qsub raxml.sh
```
##### Cut long internal branches, trim spurious tips and mask monophyletic/paraphyletic tips of the same taxon
```
python cut_long_internal_branches.py <inDIR> <internal_branch_length_cutoff> <minimal_taxa> <outDIR>
python trim_tips.py <treDIR> <outDIR> <relative_cutoff> <absolute_cutoff1> <absolute_cutoff2>
python mask_tips_by_taxonID_transcripts.py <treDIR> <aln-clnDIR> <outDIR>
```

## Ortholog Inference
##### Paralogy pruning to infer orthologs
```
python prune_paralogs_MI.py <homologDIR> <tree_ending> <relative_tip_cutoff> <absolute_tip_cutoff> <minimal_taxa> <outDIR>
python write_ortholog_fasta_files.py <fasta file with all seqs> <ortholog tree DIR> <outDIR> <MIN_TAXA>
```
##### Rename the sequence files based on Tomato Gene Model and add Capsella orthologous sequences
```
python cluster_gene_ID.py <inDIR> <treDIR> <outDIR>
python CapsellaOrtholog.py <inDIR> Tomato_Capsella.txt Capsicum.annuum.L_Zunla-1_v2.0_CDS.fa <outDIR>
```

## Alignment Construction
##### Run Guidance to make sequence alignments
```
python directory_subpackage.py <inDIR> <num_subdir> .fa
qsub guidance.sh
```
##### Re-run Guidance on unprocessed sequences
```
for file in Solyc*; do cp $file/MSA.PRANK.Without_low_SP_Col.With_Names outDIR/$file; done
python find_unprocessed_files.py <processedDIR> <originalDIR> <unprocessedDIR>
```
##### Post-alignment treatment_1, remove Capsella sequences and delete gaps or missing bases
```
qsub mask_bySW.sh
python orf_aln_process.py <inDIR> <outDIR> -s Capana -d 15
```
##### Calculate pair-wise genetic distance
```
python3.3 fasta2mvf.py --fasta alignments_Dir/* --out genes_mvf --contigbyfile --overwrite
python3.3 mvf_analyze_dna.py --mvf genes_mvf --out genetic_dist PairwiseDistanceWindow
```

## Phylogeny Construction
#### Concatenated tree and Consensus tree using RAxML
```
qsub raxml_concatenate.sh
module load phylip; consense
raxmlHPC -L MRE -z genetrees.tre -m GTRCAT -n T1
```
#### Colascence tree by ASTRAL
```
python astral_prepare.py inDIR tre_file bs_file
qsub astral.sh
```
#### Gene tree analysis with BUCKy
```
python seqformat_converter.py <inDIR> <outDIR> .phy .nex
python mrbayes_prepare.py <inDIR>
for file in *.nex; do mkdir "${file%.*nex}"; mv $file "${file%.*nex}"; done
qsub mrbayes.sh
qsub bucky.sh
```
#### Visualize phylogenetic tree
```
rstrip phylo_construct.R
```

## Introgression Analysis
##### Run ABBA using MVF
```
python3.3 fasta2mvf.py --fasta <concatenated_fasta> --out transcriptome --overwrite
python ABBA_trio.py
qsub trios.sh
```
##### Infer pairwise species-specific/common ABBA-BABA sites 
```
python ABBA_parse.py -mvf MVF_FILE -test pairwise
sh speciesID.sh
```
##### Infer direction of introgression by using D-foil test (example)
```
python3.3 mvf_analyze_dna.py --mvf transcriptome --out SIN_CAL_DAR_PRO --samples JA0702 JA0711 JA0694 JA0456 Solyc --windowsize 6201996 PatternCount
python dfoil.py --out myfile --infile SIN_CAL_DAR_PRO —pvalue 0.00001
```

## Ancestral Segregating Allele Analysis
##### Mapping reads to tomato reference genome and call SNPs 
```
sh mapping.sh
sh snp_call.sh
python mvf_join.py --mvf SL2.50ch00.mvf SL2.50ch01.mvf SL2.50ch02.mvf SL2.50ch03.mvf SL2.50ch04.mvf SL2.50ch05.mvf SL2.50ch06.mvf SL2.50ch07.mvf SL2.50ch08.mvf SL2.50ch09.mvf SL2.50ch10.mvf SL2.50ch11.mvf SL2.50ch12.mvf --out combined.mvf
```
##### Count ancestral segregating alleles
```
python ancestral_variation.py -i comibined.mvf -t species_hetero
python ancestral_variation.py -i comibined.mvf -t shared_hetero
python ancestral_variation.py -i comibined.mvf -t shared_snp
```

## Adaptive Evolution Analysis
##### Separate alignments with/without Capana and remove JA0010 from alignments
```
python orf_aln_process.py -i <inDIR> -o <outDIR> -s JA0010
grep -lir 'Capana' ./ | xargs mv -t <outDIR>
python seqformat_converter.py <inDIR> <outDIR> .fa .phy
sh edit_phy2.sh
```
##### Post-alignment treatment_2
```
python codemlScript.py <outDIR> <codeml_build> <treeFile>
qsub paml.sh
find */rub -empty -type f
python SWAMP.py -i <inDIR> -b <branchcodes.txt> -t 5 -w 15 -m 50
```
##### Remove all gaps and missing bases before PAML
```
for file in Solyc*; do cp inDIR/*masked.phy outDIR; done
python orf_aln_process.py -i <inDIR> -o <outDIR> -s seqname -d 14
```
##### Run PAML using MVF
```
python3.3 fasta2mvf.py --fasta inDIR/* --out outDIR/Jalt_ortho_dna --contigbyfile --overwrite
python3.3 mvf_translate.py --mvf Jalt_ortho_dna --out Jalt_ortho_codon
qsub mvf_paml.sh
python CombinedPAML.py <NS_out> <Geneoutput> GeneFunction.txt > PAML_final.txt
```
