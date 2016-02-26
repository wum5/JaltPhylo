# JaltPhylo

Raw scripts for the "Jaltomato Phylogenomics" project.

Still in updating!

## Folder Set
##### Python Script Folder
SF=/N/dc2/projects/jaltomt/script
##### Orthologs Folder
OF=/N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog

## Ortholog Inference
##### Rename the alignment files based on Tomato Gene Model rather than Cluster ID
python $SF/SeqRename.py $OF/initial_ortholog_align $OF/without_Capsella $OF/Cluster2Gene.txt

##### Add Capsella 1-to-1 orthologous sequence into alignment
