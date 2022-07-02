import math
from Bio import SeqIO
from tqdm import tqdm
import numpy as np, matplotlib.pyplot as plt
#test commit

filename = 'SourceFiles\hg38_chr1_and_chr2.fa'
f = open('res.txt', 'w')
f.write(f"Loaded file - {filename}\n")

seq = SeqIO.parse(filename, 'fasta')
for rec in seq:
    example = rec
    break

CPG_SIZE = 200
NON_CPG = np.array([[0.32716698441593656, 0.17248167002809267, 0.24443678155548018, 0.25591456400049056], 
                    [0.3516689128695902, 0.2579304149667452, 0.046029645994422504, 0.34437102616924214], 
                    [0.28956746321110677, 0.20887634532686009, 0.2581403123573799, 0.2434158791046533], 
                    [0.21703089912397244, 0.2051678488685245, 0.24930659876773192, 0.3284946532397712]])

CPG = np.array([[0.18831755840546174, 0.27329978209686234, 0.42616928982466556, 0.11221336967301036], 
                [0.15740548738943166, 0.36455117988944646, 0.2875502852707851, 0.19049304745033682], 
                [0.1619506626548219, 0.35192909824399826, 0.36484986132906516, 0.1212703777721147], 
                [0.08708915961511532, 0.3666073202771276, 0.35580160968746455, 0.19050191042029246]])

nucList = ['A', 'C', 'G', 'T']

dict_nuc_cpg = {'A': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'C': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'G': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'T': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'N': {'N': 0, 'A': 0,   'C': 0,    'G': 0,     'T': 0}}

dict_nuc_noncpg =  {'A': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'C': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'G': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'T': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'N': {'N': 0, 'A': 0,   'C': 0,    'G': 0,     'T': 0}}

comparison = np.full((4,4), 0.0)
comparison_ncpg = np.full((4,4), 0.0)

x = np.linspace(0, len(example.seq), len(example.seq)) 
yCpG = []
yNon = []
yAdd = []
islands_list = []
counter = 0

def MatrixInit (def_matrix, iter_num:int):
    for j in range(CPG_SIZE): #заполнение матрицы CpG первыми 300 нукл.
        def_matrix[example.seq[(j*iter_num)]][example.seq[(j*iter_num)+1]] += 1
        
    for k in range(4):
        sum_of_row = (def_matrix[nucList[k]]['A'] + 
                    def_matrix[nucList[k]]['C'] +
                    def_matrix[nucList[k]]['G'] + 
                    def_matrix[nucList[k]]['T'])
        for l in range(4):
            if sum_of_row != 0:
                comparison[k][l] = (CPG[k][l] - (def_matrix[nucList[k]][nucList[l]] / sum_of_row))**2
    return(-np.sum(comparison))

yCpG.append(MatrixInit(dict_nuc_cpg, 1))
yNon.append(MatrixInit(dict_nuc_noncpg, 1))
for i in range(CPG_SIZE):
    yAdd.append(0.01)

def CatchIsland (coord, data, param):
    coord_add = coord
    value = data[coord]
    
    if param == 1:
        value = data[coord]
        while value < 0:
            value = data[coord_add]
            coord_add -= 1
    
    if param == 2:
        while value < 0:
            value = data[coord_add]
            if coord_add+1 == len(data):
                break
            else:
                coord_add += 1
    return coord_add

for i in tqdm(range(CPG_SIZE, len(example.seq)-1), leave=True): #scaner move (CpG comp)
# for i in tqdm(range(CPG_SIZE, 1000000), leave=True):
    dict_nuc_cpg[example.seq[i-CPG_SIZE].upper()][example.seq[i-CPG_SIZE+1].upper()] -= 1
    dict_nuc_cpg[example.seq[i].upper()][example.seq[i+1].upper()] += 1
    dict_nuc_noncpg[example.seq[i-CPG_SIZE].upper()][example.seq[i-CPG_SIZE+1].upper()] -= 1
    dict_nuc_noncpg[example.seq[i].upper()][example.seq[i+1].upper()] += 1
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
            if sum_of_row != 0:
                comparison[k][l] = (CPG[k][l] - (dict_nuc_cpg[nucList[k].upper()][nucList[l].upper()] / sum_of_row))**2
            if sum_of_row_ncgp != 0:
                comparison_ncpg[k][l] = (NON_CPG[k][l] - (dict_nuc_noncpg[nucList[k].upper()][nucList[l].upper()] / sum_of_row_ncgp))**2
    yCpG.append(-np.sum(comparison))
    yNon.append(-np.sum(comparison_ncpg))
    yAdd.append(1/CPG_SIZE*math.log(np.sum(comparison)/np.sum(comparison_ncpg)))

last_call = 2*CPG_SIZE
last_write = ""

for i in range(1, len(yAdd)):
    if yAdd[i-1] < -0.01 and i - last_call >= CPG_SIZE:
        now_write = f"{int(CatchIsland(i-1, yAdd, 1) - CPG_SIZE)} \t {int(CatchIsland(i-1, yAdd, 2) - CPG_SIZE)} \n"
        if last_write != now_write:
            f.write(f"{int(CatchIsland(i-1, yAdd, 1) - CPG_SIZE)} \t {int(CatchIsland(i-1, yAdd, 2) - CPG_SIZE)} \n")
            last_write = f"{int(CatchIsland(i-1, yAdd, 1) - CPG_SIZE)} \t {int(CatchIsland(i-1, yAdd, 2) - CPG_SIZE)} \n"
            last_call = i

f.close()
fig, ax = plt.subplots()
ax.plot(yAdd)
ax.plot(np.full(len(yAdd), -0.01), color = 'gray')
plt.show()