# Jaltomata Phylogenomics

## Description
* Raw scripts/Pipeline for the "Jaltomato Phylogenomics" Project.
* Some scripts were written by Ya Yang for her study (https://bitbucket.org/yangya/phylogenomic_dataset_construction)
* Still in updating!

## Authors: 
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
##### Delete empty files
```
find . -size 0 -delete
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
##### Filter clusters with specific species
```
python species_in_clusters.py <inDIR> <outDIR>
```
##### Write sequence files from ortholog trees
```
python write_ortholog_fasta_files.py <fasta file with all seqs> <ortholog tree DIR> <outDIR> <MIN_TAXA>
```
##### Rename the sequence files based on Tomato Gene Model rather than Cluster ID
```
python cluster_gene_ID.py <inDIR> > Cluster2Gene.txt
python SeqRename.py <inDIR> <outDIR> Cluster2Gene.txt
```
##### Add Capsella-Tomato 1-to-1 orthologous sequence into sequence files
```
python CapsellaOrtholog.py <inDIR> Tomato_Capsella.txt Capsicum.annuum.L_Zunla-1_v2.0_CDS.fa <outDIR>
```

## Alignment Construction and Quality Check
##### Run Guidance to make sequence alignments
```
python directory_subpackage.py <inDIR> <num_subdir> .fa
python guidance.py <outDIR> <num_subdir>
for file in *sh; do qsub $file; done
```
##### Re-run Guidance on unprocessed sequences
```
for file in Solyc*; do cp $file/MSA.PRANK.Without_low_SP_Col.With_Names outDIR/$file; done
python find_unprocessed_files.py <processedDIR> <originalDIR> <unprocessedDIR>
```
##### Post-alignment treatment_1 (before constructing phylogeny)
```
python OrfBoundary.py <inDIR> <outDIR>
qsub mask_bySW.sh
```
##### Concatenate alignments
```
python seqformat_converter.py <phylipDIR> <fastaDIR> .phy
python ConcatSeq.py <inDIR> concat_withCap.fa
python seqformat_converter.py <fastaDIR> <phylipDIR> .fa
```

## Adaptive Evolution Analysis
##### Sliding Window before PAML
```
python seqformat_converter.py <inDIR> <outDIR> .fa
sh edit_phy2.sh
python codemlScript.py <outDIR> <codeml_build> <treeFile>
qsub paml.sh
python SWAMP.py -i <inDIR> -b <branchcodes.txt> -t 5 -w 15 -m 50
for file in Solyc*; do cp inDIR/*masked.phy outDIR; done
for file in Solyc*; do sed -i 's/N/-/g' $file; done
python DeleteSites.py <inDIR> <outDIR>
python StopCodon.py <inDIR> .fa
```
##### Run PAML using MVF
```
python3.3 fasta2mvf.py --fasta inDIR/* --out outDIR/Jalt_ortho_dna --contigbyfile --overwrite
python3.3 mvf_translate.py --mvf Jalt_ortho_dna --out Jalt_ortho_codon
qsub mvf_paml.sh
python CombinedPAML.py Clade2_out Geneoutput_Clade2 GeneFunction.txt > Clade2_final.txt
```

## Introgression Analysis
##### Run ABBA using MVF
```
python ConcatSeq.py <inDIR> Jalt_concat_dna.fa
python3.3 fasta2mvf.py --fasta Jalt_concat_dna.fa --out jalt_concat_dna --overwrite
sh trios.sh
```
