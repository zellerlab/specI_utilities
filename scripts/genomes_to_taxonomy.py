import sys
from ete3 import NCBITaxa
ncbi = NCBITaxa()

# for mOTUs > 2.5.0 and progenomes 2: https://zenodo.org/record/3357977#.YFoeUuYo-fU
ncbi_dump = "/Users/milanese/Desktop/NCBI_taxdump_2019_01_08/"
INPUT_FILE = "/Users/milanese/Desktop/temp_euk/all_genomes"
# expected:
# genome_id\tNCBI_species_id
# example
# 155920.SAMN02441076	2371
# 1559982.SAMD00020875	1559982

add_kingdom = True
# in the result it will add the kingdom to the superkingdom
# useful for euk

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
        if u == "NA":
            return "NA"
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
        # add kingdom
        if level == "kingdom" and tax_id_2_tax_level[u] == "kingdom":
            return u
        u = parent[u]

# we open the file and fill in the values --------------------------------------
o = open(INPUT_FILE,"r")
genomes_tax = dict()
genome_kingdom = dict()
count_no_ncbi = 0
count_not_species = 0
count_correct = 0
for line in o:
    vals = line.rstrip().split()
    genomes_tax[vals[0]] = ["NA","NA","NA","NA","NA","NA","NA"]
    genome_kingdom[vals[0]] = "NA"
    if not vals[1] in tax_id_2_tax_level:
        sys.stderr.write("Error: "+vals[1]+" not in the NCBI dump.\n")
        genomes_tax[vals[0]] = ["NO_NCBI_ID", vals[1]]
        count_no_ncbi = count_no_ncbi + 1
    else:
        if tax_id_2_tax_level[vals[1]] != "species":
            vals[1] = find_parent(vals[1],6)
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
            # kingdom
            genome_kingdom[vals[0]] = find_parent(vals[1], "kingdom")
o.close()

sys.stderr.write("Number of correct genomes: " +str(count_correct)+"\n")
sys.stderr.write("Number of genomes with ID that is not species level: " +str(count_not_species)+"\n")
sys.stderr.write("Number genomes with wrong ID (not in NCBI): " +str(count_no_ncbi)+"\n\n")




# we add the names -------------------------------------------------------------
to_substitute = ["","phylum incertae sedis","class incertae sedis","order incertae sedis","fam. incertae sedis","gen. incertae sedis","species incertae sedis"]
# now we fill in from top to bottom, changing the NAs
for g in genomes_tax:
    if len(genomes_tax[g]) > 2:
        last_annotated = sci_names[genomes_tax[g][0]]
        for t in range(7):
            if genomes_tax[g][t] != "NA":
                last_annotated = sci_names[genomes_tax[g][t]]
                genomes_tax[g][t] = genomes_tax[g][t] + " " + last_annotated
            else:
                genomes_tax[g][t] = "NA "+last_annotated+" "+to_substitute[t]



# we search the missing one in ETE ---------------------------------------------
def find_name(id):
    return list(ncbi.get_taxid_translator([str(id)]).values())[0]
def find_parent_ete(id,level):
    if id == "NA":
        return "NA"
    all_lineages = ncbi.get_lineage(str(id))
    for i in all_lineages:
        if list(ncbi.get_rank([str(i)]).values())[0] == level:
            return i
    return "NA"

sys.stderr.write("We try to identify the missing species in ETE\n")
for g in genomes_tax:
    if genomes_tax[g][0] == "NO_NCBI_ID":
        ID = genomes_tax[g][1]
        if not list(ncbi.get_rank([ID]).values())[0] == "species":
            ID = find_parent_ete(ID,"species")
        if list(ncbi.get_rank([ID]).values())[0] == "species":
            genome_kingdom[g] = find_parent_ete(ID,"kingdom")
            genomes_tax[g] = ["NA","NA","NA","NA","NA","NA","NA"]
            genomes_tax[g][0] = find_parent_ete(ID,"superkingdom")
            genomes_tax[g][1] = find_parent_ete(ID,"phylum")
            genomes_tax[g][2] = find_parent_ete(ID,"class")
            genomes_tax[g][3] = find_parent_ete(ID,"order")
            genomes_tax[g][4] = find_parent_ete(ID,"family")
            genomes_tax[g][5] = find_parent_ete(ID,"genus")
            genomes_tax[g][6] = ID
            # complete
            last_annotated = find_name(genomes_tax[g][0])
            for t in range(7):
                if genomes_tax[g][t] != "NA":
                    last_annotated = find_name(genomes_tax[g][t])
                    genomes_tax[g][t] = str(genomes_tax[g][t]) + " " + last_annotated
                else:
                    genomes_tax[g][t] = "NA "+last_annotated+" "+to_substitute[t]
        else:
            sys.stderr.write("Found not species ID: "+str(ID)+"\n")
            genomes_tax[g] = ["NOT_SPECIES_ID"]



# we print the result ----------------------------------------------------------
for g in genomes_tax:
    if add_kingdom:
        genomes_tax[g][0] = genomes_tax[g][0]+"("+str(genome_kingdom[g])+")"
    print(g+"\t"+ "\t".join(genomes_tax[g]))
