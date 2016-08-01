# Jaltomata Phylogenomics

## Description
* Raw scripts/Pipeline for the "Jaltomato Phylogenomics" Project.
* Some scripts were written by Ya Yang for her study (https://bitbucket.org/yangya/phylogenomic_dataset_construction)
* Still in updating!

## Authors: 
* Meng Wu
* https://github.com/wum5/JaltPhylo

## Folder Set
##### Python Script Folder
```
SF=/N/dc2/projects/jaltomt/scripts
```
##### Software Folder
```
SW=/N/dc2/projects/jaltomt/softwares
```
##### Orthologs Folder
```
OF=/N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog
```
##### Homologs Folder
```
HF=/N/dc2/projects/jaltomt/de_novo/e5_80/homologs
```

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
python $SF/fix_names_from_transdecoder.py trinity_OutFolder CDS_Folder
```
##### Reduce redundancy
```
qsub $SF/cd-hit-est.sh
```

## Homolog Inference
##### Make all-by-all blast
```
makeblastdb -in all.fa -parse_seqids -dbtype nucl -out all.fa
qsub $SF/blastn.sh
```
##### Cut ends that are fast-evolving, or using sequences from genome annotation
```
python $SF/cut_seq_ends.py all.fa all.rawblast
```
##### Inferring putative homolog groups using similarity
```
python $SF/blast_to_mcl.py all.rawblast <hit_fraction_cutoff>[0.4]
mcl all.rawblast.hit-frac0.4.minusLogEvalue --abc -te 5 -tf 'gq(10)' -I 2.5 -o hit-frac0.4_I2.5_e10
python $SF/write_fasta_files_from_mcl.py <fasta files without ends cut> <mcl_outfile> <minimal_taxa> <outDIR>
```
##### Delete empty files
```
find . -size 0 -delete
```
##### make initial alignments
```
qsub $SF/mafft.sh
qsub $SF/phyutility.sh
qsub $SF/fasttree.sh
```
##### Cut long internal branch
```
python cut_long_branches_iter.py inDIR outDIR 0.3 0.1
```
##### refine the final clusters
```
qsub $SF/mafft.sh
qsub $SF/phyutility.sh
```
##### Tree inference using RAxML
```
qsub $SF/raxml.sh
```
##### Cut long internal branches
```
python cut_long_internal_branches.py inDIR internal_branch_length_cutoff[0.2] minimal_taxa outDIR
```
##### Trim spurious tips
```
python trim_tips.py <treDIR> <outDIR> <relative_cutoff>[0.05] <absolute_cutoff1>[0.2] <absolute_cutoff1>[0.1]
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
python write_ortholog_fasta_files.py <fasta file with all seqs> <ortholog tree DIR> outDIR MIN_TAXA
```
##### Rename the sequence files based on Tomato Gene Model rather than Cluster ID
```
python $SF/SeqRename.py $OF/initial_ortholog_align $OF/without_Capsella $OF/Cluster2Gene.txt
```
##### Add Capsella-Tomato 1-to-1 orthologous sequence into sequence files
```
python $SF/CapsellaOrtholog.py $OF/without_Capsella $OF/Tomato_Capsella.txt $OF/Capsicum.annuum.L_Zunla-1_v2.0_CDS.fa ```
$OF/with_Capsella
```

## Alignment Construction and Quality Check
##### Run Guidance to make sequence alignments
```
python $SF/directory_subpackage.py $OF/with_Capsella/ 40 .fa
python build.py
python qsub_cmd.py
sh qsub.sh
```

##### Re-run Guidance on unprocessed sequences
```
cd $OF/guidAlign2
for file in Solyc*; do cp $file/MSA.PRANK.Without_low_SP_Col.With_Names ../post-guid/init/$file; done
python $SF/find_unprocessed_files.py $OF/post-guid/init $OF/with_Capsella $OF/unprocessed
qsub prank_plus.sh
```

##### Post-alignment treatment_1 (before constructing phylogeny)
```
python $SF/OrfBoundary.py $OF/post-guid/init $OF/post-guid/2nd_BoundaryFixed
qsub $SF/mask_bySW.sh
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
