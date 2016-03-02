# Jaltomata Phylogenomics

## Description
* Raw scripts/Pipeline for the "Jaltomato Phylogenomics" Project.
* Some scripts were written by Ya Yang for her study (https://bitbucket.org/yangya/phylogenomic_dataset_construction)
* Still in updating!

## AUTHORS: 
* Meng Wu
* https://github.com/wum5/JaltPhylo

## Folder Set
##### Python Script Folder
SF=/N/dc2/projects/jaltomt/script
##### Orthologs Folder
OF=/N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog
##### Homologs Folder
HF=/N/dc2/projects/jaltomt/de_novo/e5_80/homologs

## Ortholog Inference
##### Paralogy pruning to infer orthologs
python prune_paralogs_MI.py <homologDIR> <tree_file_ending> <relative_long_tip_cutoff> <absolute_long_tip_cutoff> <minimal_taxa> <outDIR>
##### Write sequence files from ortholog trees
python write_ortholog_fasta_files.py <fasta file with all seqs> <ortholog tree DIR> outDIR MIN_TAXA
##### Rename the sequence files based on Tomato Gene Model rather than Cluster ID
python $SF/SeqRename.py $OF/initial_ortholog_align $OF/without_Capsella $OF/Cluster2Gene.txt
##### Add Capsella-Tomato 1-to-1 orthologous sequence into sequence files
python $SF/CapsellaOrtholog.py $OF/without_Capsella $OF/Tomato_Capsella.txt $OF/Capsicum.annuum.L_Zunla-1_v2.0_CDS.fa $OF/with_Capsella

## Alignment Construction and Quality Check
##### Run Guidance to make sequence alignments
* python $SF/directory_subpackage.py $OF/with_Capsella/ 40 .fa
* python build.py
* python qsub_cmd.py
* sh qsub.sh

##### re-run Guidance on unprocessed sequences
* cd $OF/guidAlign2
* for file in Solyc*; do cp $file/MSA.PRANK.Without_low_SP_Col.With_Names ../post-guid/init/$file; done
* python $SF/find_unprocessed_files.py $OF/post-guid/init $OF/with_Capsella $OF/unprocessed
* sh prank_plus.sh

##### post-alignment treatment
* python $SF/OrfBoundary.py $OF/post-guid/init $OF/2nd_BoundaryFixed
* 
