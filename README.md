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
* [Adaptive Evolution Analysis](#adaptive-evolution-analysis)

## Overview
* Raw scripts/Pipeline for the "Jaltomato Phylogenomics" Project.
* Some scripts were written by Ya Yang for her study (https://bitbucket.org/yangya/phylogenomic_dataset_construction)
* Still in updating!

## Contributors 
* Meng Wu
* https://github.com/wum5/JaltPhylo

## Raw Data Processing
##### Trim low-quality reads using shear.py
```
qsub trim.sh
```
##### Remove the first 15-bp of reads due to non-random hexamer primers
```
qsub clip5end.sh
```
##### Check read quality using FastaQC
```
qsub FastaQC.sh
```
##### Build transcript assembly
```
qsub trinity.sh
```
##### Predict CDS using Transdecoder
```
qsub transdecoder.sh
```
##### Rename CDS files
```
for file in *_dir; do cp $file/longest_orfs.cds outDIR/$file'.cds'; done
python fix_names_from_transdecoder.py <DIR> <DIR>
```
##### Reduce redundancy
```
cat *_NR.fa *_RP.fa > *.cds.fa
qsub cd-hit-est.sh
```

## Homolog Inference
##### Make all-by-all blast
```
qsub blastn.sh
```
##### Inferring putative homolog groups using similarity
```
cat *blastn >all.rawblast
python blast_to_mcl.py all.rawblast <hit_fraction_cutoff>
mcl all.rawblast.hit-frac0.4.minusLogEvalue --abc -te 5 -tf 'gq(10)' -I 2.5 -o hit-frac0.4_I2.5_e10
python write_fasta_files_from_mcl.py <fasta files> <mcl_outfile> <minimal_ingroup_taxa> <outDIR>
```
##### make initial alignments
```
qsub mafft.sh
qsub phyutility.sh
qsub fasttree.sh
```
##### Cut long internal branch
```
python cut_long_branches_iter.py <inDIR> <outDIR>
```
##### refine the final clusters
```
qsub mafft.sh
qsub phyutility.sh
qsub raxml.sh
```
##### Cut long internal branches
```
python cut_long_internal_branches.py <inDIR> <internal_branch_length_cutoff> <minimal_taxa> <outDIR>
```
##### Trim spurious tips
```
python trim_tips.py <treDIR> <outDIR> <relative_cutoff> <absolute_cutoff1> <absolute_cutoff2>
```
##### Mask monophyletic and paraphyletic tips that belongs to the same taxon
```
python mask_tips_by_taxonID_transcripts.py <treDIR> <aln-clnDIR> <outDIR>
```

## Ortholog Inference
##### Paralogy pruning to infer orthologs
```
python prune_paralogs_MI.py <homologDIR> <tree_file_ending> <relative_long_tip_cutoff> <absolute_long_tip_cutoff> <minimal_taxa> <outDIR>
```
##### Write sequence files from ortholog trees
```
python write_ortholog_fasta_files.py <fasta file with all seqs> <ortholog tree DIR> <outDIR> <MIN_TAXA>
```
##### Rename the sequence files based on Tomato Gene Model rather than Cluster ID
```
python cluster_gene_ID.py <inDIR> <treDIR> <outDIR>
```
##### Add Capsella-Tomato 1-to-1 orthologous sequence into sequence files
```
python CapsellaOrtholog.py <inDIR> Tomato_Capsella.txt Capsicum.annuum.L_Zunla-1_v2.0_CDS.fa <outDIR>
```

## Alignment Construction
##### Run Guidance to make sequence alignments
```
python directory_subpackage.py <inDIR> <num_subdir> .fa
python guidance_process.py <outDIR> <num_subdir> <thread> <hour>
for file in *sh; do qsub $file; done
```
##### Re-run Guidance on unprocessed sequences
```
for file in Solyc*; do cp $file/MSA.PRANK.Without_low_SP_Col.With_Names outDIR/$file; done
python find_unprocessed_files.py <processedDIR> <originalDIR> <unprocessedDIR>
```
##### Post-alignment treatment_1 (before constructing phylogeny)
```
qsub mask_bySW.sh
```
##### Remove Capana sequences and delete gaps or missing bases from alignments
```
python orf_aln_process.py <inDIR> <outDIR>
```

## Phylogeny Construction
#### Concatenated tree
```
python ConcatSeq.py <inDIR> concat_withCap.fa
python seqformat_converter.py <fastaDIR> <phylipDIR> .fa
qsub raxml_concatenate.sh
```

## Introgression Analysis
##### Run ABBA using MVF
```
python ConcatSeq.py <inDIR> Jalt_concat_dna.fa
python3.3 fasta2mvf.py --fasta Jalt_concat_dna.fa --out jalt_concat_dna --overwrite
sh trios.sh
```

## Adaptive Evolution Analysis
##### Post-alignment treatment_2 (before PAML)
```
python seqformat_converter.py <inDIR> <outDIR> .fa
sh edit_phy2.sh
python codemlScript.py <outDIR> <codeml_build> <treeFile>
qsub paml.sh
python SWAMP.py -i <inDIR> -b <branchcodes.txt> -t 5 -w 15 -m 50
for file in Solyc*; do cp inDIR/*masked.phy outDIR; done
for file in Solyc*; do sed -i 's/N/-/g' $file; done
python orf_aln_process.py <inDIR> <outDIR> seqname
```
##### Run PAML using MVF
```
python3.3 fasta2mvf.py --fasta inDIR/* --out outDIR/Jalt_ortho_dna --contigbyfile --overwrite
python3.3 mvf_translate.py --mvf Jalt_ortho_dna --out Jalt_ortho_codon
qsub mvf_paml.sh
python CombinedPAML.py Clade2_out Geneoutput_Clade2 GeneFunction.txt > Clade2_final.txt
```
