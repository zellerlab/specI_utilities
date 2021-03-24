import pandas as pd
import numpy as np
from collections import defaultdict
import os
import re
import sys 

input_path_to_files = sys.argv[1]
output_name_of_file = sys.argv[2]
# example: python evaluate_clustering.py /g/scb2/bork/mende/specI_euks/vsearch//combined3//clustering_comb3/ evaluation_clusterings.csv

def main():
    path_to_files = input_path_to_files
    all_files = os.listdir(path_to_files)
    all_files.sort(reverse=True)
    scores_all_cluster = pd.DataFrame(columns = ['filename','treshold', 'total_cluster','TP', 'FP', 'FN','precision', 'recall', 'f1_score'])
    for file in all_files:
        
        ####### Create two different dictionaries; all_cluster lists for each cluster number all IDs;  all_ids lists for each ID all cluster numbers, which contain the specific ID 
        all_cluster = defaultdict(list)
        all_ids = defaultdict(list)
        with open(path_to_files+file,'r') as f:
            for n_cluster,line in enumerate(f):
                for genome in line.split(';'):
                    all_cluster[n_cluster].append(genome.split('|')[2].strip())
                    all_ids[genome.split('|')[2].strip()].append(n_cluster)
                    
        ####### Identify False Positives: False Positives are identified as clusters containing different IDs; different species are merged into one cluster
        equal_ids = []
        different_ids = []
        for cl,list_ids in all_cluster.items():
            if len(set(list_ids))== 1:
                equal_ids.append(cl)
            else:
                different_ids.append(cl)

        ####### Identify False Negatives: False Negatives are identified as clusters containing IDs, which appear in several clusters; the same species is split into several clusters
        equal_clust = []
        different_clust = []
        clust_ids_oneclust = []
        clust_ids_diffclust = []
        for id_n,list_clust in all_ids.items():
            if len(set(list_clust))== 1:
                equal_clust.append(list_clust[0])
            else:
                for i in set(list_clust):
                    different_clust.append(i)
        for cl in range(len(all_cluster)):
            if cl in different_clust:
                clust_ids_diffclust.append(cl)
            else:
                clust_ids_oneclust.append(cl)

        TP = len(np.setdiff1d(clust_ids_oneclust,different_ids)) # TP: all clusters, which do not belong to FP or FN
        FP = len(different_ids)
        FN = len(np.setdiff1d(clust_ids_diffclust,different_ids)) # remove clusters, which also contain merged species and are listed already in FP 

        precision = TP / (TP + FP)
        recall = TP / (TP + FN)
        f1_score = 2 * precision * recall / (precision + recall)

        m = re.search('0.[0-9]+', file) #extract treshold from filename
        treshold = str(1.0 - float(m.group(0))) 
        
        scores_all_cluster = scores_all_cluster.append({'filename':file, 'treshold':treshold, 'total_cluster': len(all_cluster), 'TP':TP, 'FP':FP, 'FN':FN, 'precision':precision, 'recall':recall, 'f1_score':f1_score}, ignore_index=True)
        scores_all_cluster.to_csv(output_name_of_file, sep='\t', index=False)

if __name__ == '__main__':
  main()