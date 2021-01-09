#plot confusion matrix for DRUID output to assess its accuracy
import seaborn
import matplotlib.pyplot as plt
import numpy as np
import math
import argparse
from collections import defaultdict

# acronym names used by DRUID for some very close relatives
# see github page for details
# I will use 0 to denote unrelated pairs
relationship = {'FS':1, 'PC':1, 'P':1, 'C':1, 'HS':2, 'GP':2, 'GC':2, 'AU':2, 'NN':2, 'UN':-1, 'DC':2, 'MZ':0, 'AV':2}
N_COLUMNS = 11
true_label_file = '/fs/cbsubscb09/storage/yilei/simulate/SAMAFS/safs.ped.labels'
infer_label_file = './safs.DRUID'
removed_inds = ['802213', 'A09222', 'A19021', 'A38060', 'A45074', 'A28102',
'840914', 'A09222', 'A27186', 'A28171', '808134']

def read_true_labels(file):
    rels = defaultdict(lambda : {})
    with open(file) as f:
        f.readline()
        line = f.readline()
        while line:
            _, _, id1, id2, degree, _, _ = line.strip().split('\t')
            if degree == 'Inf':
                degree = -1
            elif degree.find('.') != -1:
                line = f.readline()
                continue
            else:
                degree = int(degree)
            id1, id2 = min(id1, id2), max(id1, id2)
            if degree > 10:
                degree = -1
            rels[id1][id2] = degree
            line = f.readline()
    # take care of twins
    rels['A03081']['A03082'] = 0
    rels['A03442']['A03443'] = 0
    rels['813139']['813140'] = 0
    rels['A00336']['A07012'] = 0
    rels['816230']['816231'] = 0
    rels['A23133']['A23134'] = 0
    return rels

def read_druid_labels(file):
    rels = defaultdict(lambda : {})
    with open(file) as f:
        f.readline() #ignore header line
        line = f.readline()
        while line:
            id1, id2, degree = line.strip().split('\t')
            id1, id2 = min(id1, id2), max(id1, id2)
            rels[id1][id2] = relationship[degree] if str.isalpha(degree) else int(degree)
            line = f.readline()
    return rels

def main():
    matrix = np.zeros((N_COLUMNS, N_COLUMNS))
    true_labels = read_true_labels(true_label_file)
    inferred_labels = read_druid_labels(infer_label_file)
    for id1 in true_labels.keys():
        for id2 in true_labels[id1].keys():
            if id1 in removed_inds or id2 in removed_inds:
                continue
            dr_ref = true_labels[id1][id2]
            dr_ref = 10 if dr_ref >= 10 else dr_ref
            if id1 not in inferred_labels or id2 not in inferred_labels[id1]:
                dr_inferred = 10
            else:
                dr_inferred = inferred_labels[id1][id2]
            matrix[dr_ref, dr_inferred] += 1

    #sum by row
    num_pairs_per_category = np.sum(matrix, axis=1)
    #divide by row sum, so each row should now sum up to 1
    matrix = matrix/num_pairs_per_category[:,np.newaxis]
    labels = ['MZ','1','2','3','4','5','6','7','8', '9', '>=10']
    mask_matrix = (matrix < 0.005)
    fig, ax = plt.subplots(figsize=(16,15))
    ax = seaborn.heatmap(matrix, annot=True, xticklabels=labels, yticklabels=labels, mask=mask_matrix,
                          fmt='.2f',cmap='Oranges')
    plt.title(f'Inference for SAMAFS using DRUID', fontsize=20, fontweight='bold')
    plt.xlabel('Inferred Degree of Relatedness', fontsize=15, fontweight='bold')
    plt.ylabel('Reference Degree of Relatedness', fontsize=15, fontweight='bold')
    plt.savefig(f'safs.druid.png', dpi=300)

if __name__ == '__main__':
    main()
