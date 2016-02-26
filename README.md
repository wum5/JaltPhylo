# JaltPhylo

## Description
* Raw scripts for the "Jaltomato Phylogenomics" project.
* Still in updating!

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
##### Write alignment files from ortholog trees
python write_ortholog_fasta_files.py <fasta file with all seqs> <ortholog tree DIR> outDIR MIN_TAXA
##### Rename the alignment files based on Tomato Gene Model rather than Cluster ID
python $SF/SeqRename.py $OF/initial_ortholog_align $OF/without_Capsella $OF/Cluster2Gene.txt
##### Add Capsella 1-to-1 orthologous sequence into alignment
