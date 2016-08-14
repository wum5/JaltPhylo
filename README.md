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
qsub $SF/trim.sh
```
##### Remove the first 15-bp of reads due to non-random hexamer primers
```
qsub $SF/clip5end.sh
```
##### Check read quality using FastaQC
```
qsub $SF/FastaQC.sh
```
##### Build transcript assembly
```
qsub $SF/trinity.sh
```
##### Predict CDS using Transdecoder
```
qsub $SF/transdecoder.sh
```
##### Rename CDS files
```
for file in *_dir; do cp $file/longest_orfs.cds outDIR/$file'.cds'; done
python $SF/fix_names_from_transdecoder.py <DIR> <DIR>
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
```
##### Tree inference using RAxML
```
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
for file in Solyc*; do cp $file/MSA.PRANK.Without_low_SP_Col.With_Names <outDIR>; done
python find_unprocessed_files.py <processedDIR> <originalDIR> <unprocessedDIR>
qsub prank_plus.sh
```

##### Post-alignment treatment_1 (before constructing phylogeny)
```
python OrfBoundary.py $OF/post-guid/init $OF/post-guid/2nd_BoundaryFixed
qsub mask_bySW.sh
```

##### Post-alignment treatment_2 (before PAML analysis)
```
python $SF/seqformat_converter.py $OF/post-guid/3rd_maskedSW $OF/post-guid/4th_preSWAMP .fa
sh edit_phy2.sh
python $SF/codemlScript.py $OF/post-guid/4th_preSWAMP $OF/codeml_build $OF/treeFile
python $SF/directory_subpackage.py $OF/post-guid/4th_preSWAMP 1000 fileEnding
qsub paml.sh
python software/SWAMP-master/SWAMP.py -i $OF/post-guid/4th_preSWAMP -b $OF/branchcodes_withCap.txt -t 5 -w 15 -m 50
for file in Solyc*; do cp $file/*masked.phy $OF/post-guid/5th_postSWAMP/; done
for file in $OF/post-guid/5th_postSWAMP/fastaFile/Solyc*; do cp $file $OF/post-guid/6th_Gblocks; done
for file in Solyc*; do sed -i 's/N/-/g' $file; done
python $SF/DeleteSites.py $OF/post-guid/5th_postSWAMP/fastaFile/ $OF/post-guid/6th_Final/
python $SF/StopCodon.py $OF/post-guid/6th_Final/ .fa
```

##### Concatenate alignments
```
python $SF/seqformat_converter.py $OF/post-guid/5th_postSWAMP/phylipFile/ $OF/post-guid/5th_postSWAMP/fastaFile/ .phy
python $SF/ConcatSeq.py $OF/post-guid/5th_postSWAMP/fastaFile/ $OF/post-guid/5th_postSWAMP/concat_withCap.fa
python $SF/seqformat_converter.py $OF/post-guid/5th_postSWAMP $OF/post-guid/5th_postSWAMP .fa
```

## Adaptive Evolution Analysis
##### Run PAML using MVF
```
python3.3 $SW/mvftools-dev-master/fasta2mvf.py --fasta $OF/post-guid/6th_Final/* --out $OF/MVF_PAML/withCap/Jalt_ortho_dna --contigbyfile --overwrite
python3.3 $SW/mvftools-dev-master/mvf_translate.py --mvf $OF/MVF_PAML/withCap/Jalt_ortho_dna --out $OF/MVF_PAML/withCap/Jalt_ortho_codon
qsub $SF/mvf_paml.sh
python $SF/CombinedPAML.py $OF/MVF_PAML/withCap/Clade2_out $OF/MVF_PAML/withCap/Geneoutput_Clade2 $OF/GeneFunction.txt > $OF/MVF_PAML/Clade2_final.txt
```

##### Run PAML on Untested Candidate Genes
```
sh $SF/prank/prank_plus.sh
```

## Introgression Analysis
##### Run ABBA using MVF
```
python $SF/ConcatSeq.py $OF/pranked_align/6th_Final $OF/D-stat/Jalt_concat_dna.fa
python3.3 $SW/mvftools-dev-master/fasta2mvf.py --fasta $OF/D-stat/Jalt_concat_dna.fa --out $OF/D-stat/jalt_concat_dna --overwrite
sh $OF/D-stat/trios.sh
```
