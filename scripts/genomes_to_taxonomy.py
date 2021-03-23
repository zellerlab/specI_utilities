import sys

# for mOTUs > 2.5.0 and progenomes 2: https://zenodo.org/record/3357977#.YFoeUuYo-fU
ncbi_dump = "~/NCBI_taxdump_2019_01_08/"
INPUT_FILE = "~/genomes_info"
# expected:
# genome_id\tNCBI_species_id
# example
# 155920.SAMN02441076	2371
# 1559982.SAMD00020875	1559982



################################################################################
# we save the parent and the level information
# FIRST: load the map from tax_num_id to tax level -----------------------------
#        and the parent
# 2 -> superkingdom
# 2 -> 1
o = open(ncbi_dump+"nodes.dmp","r")
tax_id_2_tax_level = dict()
parent = dict()
for i in o:
    vals = i.rstrip().split("\t|\t")
    tax_id_2_tax_level[vals[0]] = vals[2]
    parent[vals[0]] = vals[1]
o.close()

# SECOND: scientific names for the NCBI tax IDs --------------------------------
# 2 -> Bacteria
# note that the ID "2" have different names (all saved in names), but here we
# select one
o = open(ncbi_dump+"names.dmp","r")
sci_names = dict()
available_sci_names = dict()
for i in o:
    vals = i.rstrip().split("\t|\t")
    if vals[3] == "scientific name\t|":
        sci_names[vals[0]] = vals[1]
        available_sci_names[vals[1]] = vals[0]

# useful function --------------------------------------------------------------
def find_parent(current_clade,level):
    # from the current clade, we want to find a proper annotation for the
    # level passed as input
    u = parent[current_clade]
    while True:
        if u == "1":
            return "NA"
        if level == 0 and tax_id_2_tax_level[u] == "superkingdom":
            return u
        if level == 1 and tax_id_2_tax_level[u] == "phylum":
            return u
        if level == 2 and tax_id_2_tax_level[u] == "class":
            return u
        if level == 3 and tax_id_2_tax_level[u] == "order":
            return u
        if level == 4 and tax_id_2_tax_level[u] == "family":
            return u
        if level == 5 and tax_id_2_tax_level[u] == "genus":
            return u
        if level == 6 and tax_id_2_tax_level[u] == "species":
            return u
        u = parent[u]

# we open the file and fill in the values --------------------------------------
o = open(INPUT_FILE,"r")
genomes_tax = dict()
count_no_ncbi = 0
count_not_species = 0
count_correct = 0
for line in o:
    vals = line.rstrip().split()
    genomes_tax[vals[0]] = ["NA","NA","NA","NA","NA","NA","NA"]
    if not vals[1] in tax_id_2_tax_level:
        sys.stderr.write("Error: "+vals[1]+" not in the NCBI dump.\n")
        genomes_tax[vals[0]] = ["NO_NCBI_ID"]
        count_no_ncbi = count_no_ncbi + 1
    else:
        if tax_id_2_tax_level[vals[1]] != "species":
            sys.stderr.write("Error: "+vals[1]+" not a species tax id.\n")
            genomes_tax[vals[0]] = ["NOT_SPECIES_ID"]
            count_not_species = count_not_species + 1
        else:
            count_correct = count_correct + 1
            genomes_tax[vals[0]][6] = vals[1]
            genomes_tax[vals[0]][5] = find_parent(vals[1],5)
            genomes_tax[vals[0]][4] = find_parent(vals[1],4)
            genomes_tax[vals[0]][3] = find_parent(vals[1],3)
            genomes_tax[vals[0]][2] = find_parent(vals[1],2)
            genomes_tax[vals[0]][1] = find_parent(vals[1],1)
            genomes_tax[vals[0]][0] = find_parent(vals[1],0)
o.close()

sys.stderr.write("Number of correct genomes: " +str(count_correct)+"\n")
sys.stderr.write("Number of genomes with ID that is not species level: " +str(count_not_species)+"\n")
sys.stderr.write("Number genomes with wrong ID (not in NCBI): " +str(count_no_ncbi)+"\n")



# we add the names -------------------------------------------------------------
to_substitute = ["","phylum incertae sedis","class incertae sedis","order incertae sedis","fam. incertae sedis","gen. incertae sedis","species incertae sedis"]
# now we fill in from top to bottom, changing the NAs
for g in genomes_tax:
    if len(genomes_tax[g]) != 1:
        last_annotated = sci_names[genomes_tax[g][0]]
        for t in range(7):
            if genomes_tax[g][t] != "NA":
                last_annotated = sci_names[genomes_tax[g][t]]
                genomes_tax[g][t] = genomes_tax[g][t] + " " + last_annotated
            else:
                genomes_tax[g][t] = "NA "+last_annotated+" "+to_substitute[t]
    # we print the result
    print(g+"\t"+ "\t".join(genomes_tax[g]))
