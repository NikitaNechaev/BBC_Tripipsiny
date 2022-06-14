import numpy as np, matplotlib.pyplot as plt
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd

#filename = 'BioBootCamp 2022, 2-ой этап\generated_sequence.fa'
filename = 'BioBootCamp 2022, 2-ой этап\hg38_chr1_and_chr2.fa'
#filename = 'myFasta.fa'

seq = SeqIO.parse(filename, 'fasta')
for rec in seq:
    example = rec
    break

exel_data = pd.read_excel('BioBootCamp 2022, 2-ой этап\CleanCoords.xlsx')
islands_coords = pd.DataFrame(exel_data).to_numpy()

# образцовые целевые матрицы (generated_sequence) :
NON_CPG = np.array([[0.300, 0.205, 0.285, 0.210], [0.322, 0.298, 0.078, 0.302], [0.248, 0.246, 0.298, 0.208], [0.177, 0.239, 0.292, 0.292]])
CPG = np.array([[0.180, 0.274, 0.426, 0.120], [0.171, 0.368, 0.274, 0.188], [0.161, 0.339, 0.375, 0.125], [0.079, 0.355, 0.384, 0.182]])

MIN_ISLAND_SIZE = 200
COUNT_OF_NUC_IN_CPG = 1942557

nucList = ['A', 'C', 'G', 'T']

dict_nuc = {'A': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'C': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'G': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'T': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
            'N': {'N': 0}}

dict_nuc_noncpg =  {'A': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
                    'C': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
                    'G': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
                    'T': {'A': 0,   'C': 0,    'G': 0,     'T': 0},
                    'N': {'N': 0}}

for i in tqdm(range(len(islands_coords)-1)):  # считаем количество переходов побуквенно (кол-во AA, CG, ...)
    for j in range(islands_coords[i][0], islands_coords[i][1]):
        dict_nuc[example.seq[j].upper()][example.seq[j+1].upper()] += 1
    for k in range(islands_coords[i][1], islands_coords[i+1][0]):
        dict_nuc_noncpg[example.seq[j].upper()][example.seq[j+1].upper()] += 1

     # уникальные случаи (от 0 до первого числа, 2 последних)   
# for j1 in range(islands_coords[-1][0], islands_coords[-1][1]):
#     dict_nuc[example.seq[j1].upper()][example.seq[j1+1].upper()] += 1
# for j2 in range(islands_coords[0][0]):
#     dict_nuc_noncpg[example.seq[j2].upper()][example.seq[j2+1].upper()] += 1
# for j3 in range(islands_coords[-1][1], len(example.seq)):
#     dict_nuc_noncpg[example.seq[j3].upper()][example.seq[j3+1].upper()] += 1

for i in range(4):      # считаем по формуле - значение перехода из X в Y опр. как отношение числа перех. из X в Y к сумме перех. из X
    sum_of_row = (dict_nuc[nucList[i]]['A'] + 
                dict_nuc[nucList[i]]['C'] +
                dict_nuc[nucList[i]]['G'] + 
                dict_nuc[nucList[i]]['T'])

    sum_of_row_nc = (dict_nuc_noncpg[nucList[i]]['A'] + 
                dict_nuc_noncpg[nucList[i]]['C'] +
                dict_nuc_noncpg[nucList[i]]['G'] + 
                dict_nuc_noncpg[nucList[i]]['T'])
    for j in range(4):
        dict_nuc[nucList[i]][nucList[j]] /= sum_of_row
        dict_nuc_noncpg[nucList[i]][nucList[j]] /= sum_of_row_nc

print(dict_nuc)
print(dict_nuc_noncpg)