# specI_utilities

### `genomes_to_taxonomy.py`
Script that given a file like `genome_ID\tNCBI_species_ID` print the full taxonomy. Example input:
```
155920.SAMN02441076	1562720
1559982.SAMD00020875	1559982
```
Another input is necessary: the NCBI tax dump (like: https://zenodo.org/record/3357977#.YFoeUuYo-fU).
The input should be changed in the script and the result is printed to the stdout:
```
155920.SAMN02441076	2 Bacteria	1117 Cyanobacteria	NA Cyanobacteria class incertae sedis	1118 Chroococcales	1890452 Cyanobacteriaceae	102234 Cyanobacterium	1562720 Cyanobacterium sp. IPPAS B-1200
1559982.SAMD00020875	2 Bacteria	201174 Actinobacteria	1760 Actinobacteria	85011 Streptomycetales	2062 Streptomycetaceae	1883 Streptomyces	1559982 Streptomyces sp. NBRC 109436
```
Note that we fill in the missing levels (example: `NA Cyanobacteria class incertae sedis`)
There are two possible issues:
1. `NOT_SPECIES_ID` if the tax id is not at species level
2. `NO_NCBI_ID` if the speciesID is not in the NCBI dump

### `evaluate_clustering.py`
Script that given a folder with different clusterings creates an evaluation table with precision, recall and f1_score. As input arguments it requires the path to the clustering files and the name of the evaluation table. Moreover, the script requires the installation of pandas and numpy. 
```
python evaluate_clustering.py path/to/clustering/files output_filename
