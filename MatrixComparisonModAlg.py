from Bio import SeqIO
from tqdm import tqdm
import numpy as np, matplotlib.pyplot as plt

filename = 'SourceFiles\generated_sequence.fa'

seq = SeqIO.parse(filename, 'fasta')
for rec in seq:
    example = rec
    break

CPG_SIZE = 300
NON_CPG = np.array([[0.300, 0.205, 0.285, 0.210], [0.322, 0.298, 0.078, 0.302], [0.248, 0.246, 0.298, 0.208], [0.177, 0.239, 0.292, 0.292]])
CPG = np.array([[0.180, 0.274, 0.426, 0.120], [0.171, 0.368, 0.274, 0.188], [0.161, 0.339, 0.375, 0.125], [0.079, 0.355, 0.384, 0.182]])
#целевые матрицы (МЦ)

nucList = ['A', 'C', 'G', 'T']

dict_nuc_cpg = {'A': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'C': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'G': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'T': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'N': {'N': 0}}

dict_nuc_noncpg =  {'A': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'C': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'G': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'T': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'N': {'N': 0}}

comparison = np.full((4,4), 0.0)
comparison_ncpg = np.full((4,4), 0.0)

x = np.linspace(0, len(example.seq), len(example.seq)) 
yCpG = []
yNon = []

def MatrixInit (def_matrix, iter_num:int):
    for j in range(CPG_SIZE): #заполнение матрицы CpG первыми 300 нукл.
        def_matrix[example.seq[(j*iter_num)]][example.seq[(j*iter_num)+1]] += 1
        
    for k in range(4):
        sum_of_row = (def_matrix[nucList[k]]['A'] + 
                    def_matrix[nucList[k]]['C'] +
                    def_matrix[nucList[k]]['G'] + 
                    def_matrix[nucList[k]]['T'])
        for l in range(4):
            comparison[k][l] = (CPG[k][l] - (def_matrix[nucList[k]][nucList[l]] / sum_of_row))**2
    return(-np.sum(comparison))

yCpG.append(MatrixInit(dict_nuc_cpg, 1))
yNon.append(MatrixInit(dict_nuc_noncpg, 1))

for i in tqdm(range(CPG_SIZE, len(example.seq)-1), leave=True): #scaner move (CpG comp)
    dict_nuc_cpg[example.seq[i-CPG_SIZE]][example.seq[i-CPG_SIZE+1]] -= 1
    dict_nuc_cpg[example.seq[i]][example.seq[i+1]] += 1
    dict_nuc_noncpg[example.seq[i-CPG_SIZE]][example.seq[i-CPG_SIZE+1]] -= 1
    dict_nuc_noncpg[example.seq[i]][example.seq[i+1]] += 1
    for k in range(4):
        sum_of_row = (dict_nuc_cpg[nucList[k]]['A'] + 
                dict_nuc_cpg[nucList[k]]['C'] +
                dict_nuc_cpg[nucList[k]]['G'] + 
                dict_nuc_cpg[nucList[k]]['T'])
        sum_of_row_ncgp = (dict_nuc_noncpg[nucList[k]]['A'] + 
                dict_nuc_noncpg[nucList[k]]['C'] +
                dict_nuc_noncpg[nucList[k]]['G'] + 
                dict_nuc_noncpg[nucList[k]]['T'])
        for l in range(4):
            comparison[k][l] = (CPG[k][l] - (dict_nuc_cpg[nucList[k]][nucList[l]] / sum_of_row))**2
            comparison_ncpg[k][l] = (NON_CPG[k][l] - (dict_nuc_noncpg[nucList[k]][nucList[l]] / sum_of_row))**2
    yCpG.append(-np.sum(comparison))
    yNon.append(-np.sum(comparison_ncpg))

dyn_row = []
dyn_row_coords = []
islands_list = []

for i in range(len(yCpG)-1):
    while yCpG[i] > yNon[i]:
        dyn_row.append(yCpG[i])
        i+=1

fig, ax = plt.subplots()
ax.plot(yCpG)
ax.plot(yNon)
plt.show()