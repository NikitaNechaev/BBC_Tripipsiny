# Импорт библиотек:
# math - вычисление логарифма, 
# Biopython - чтение, обработка данных FASTA файлов, 
# tqdm - вывод прогрессбара в консоль (удобное отслеживание процесса выполнения скрипта),
# numpy - работа с матрицами, массивами, некоторые вычисления, 
# matplotlib - работа с графиками (визуализация выходных данных)
import math
from Bio import SeqIO
from tqdm import tqdm
import numpy as np, matplotlib.pyplot as plt

# Инициализация FASTA файла, файла с выводом координат обнаруженных островков
filename = 'SourceFiles\hg38_chr1_and_chr2 copy.fa' # используется копия исходного файла только с фрагментом последовательности chr2 (chr1 удален)
f = open('res.txt', 'w')
f.write(f"Loaded file - {filename}\n")
seq = SeqIO.parse(filename, 'fasta')
for rec in seq:
    example = rec
    break

# Инициализация переменных и массивов
# Создаем константы, определяющие размер CpG островка, матрицы переходов нуклеотидов, характерные для CpG островков и для не CpG участков
triggeer_point = 0 # Переменная, задающее значение выходных данных, после которого фрагмент можно считать CGI
CPG_SIZE = 200
NON_CPG = np.array([[0.32716698441593656, 0.17248167002809267, 0.24443678155548018, 0.25591456400049056], 
                    [0.3516689128695902, 0.2579304149667452, 0.046029645994422504, 0.34437102616924214], 
                    [0.28956746321110677, 0.20887634532686009, 0.2581403123573799, 0.2434158791046533], 
                    [0.21703089912397244, 0.2051678488685245, 0.24930659876773192, 0.3284946532397712]])

CPG = np.array([[0.18831755840546174, 0.27329978209686234, 0.42616928982466556, 0.11221336967301036], 
                [0.15740548738943166, 0.36455117988944646, 0.2875502852707851, 0.19049304745033682], 
                [0.1619506626548219, 0.35192909824399826, 0.36484986132906516, 0.1212703777721147], 
                [0.08708915961511532, 0.3666073202771276, 0.35580160968746455, 0.19050191042029246]])

# Словари dict_nuc_cpg и dict_nuc_noncpg позволят записывать количество подсчитанных комбинаций динуклеотидов
# (e.g. последовательность ...ATTC... будет записана как dict_nuc_cpg['A']['T'] += 1, dict_nuc_cpg['T']['T'] += 1, dict_nuc_cpg['T']['C'] += 1)
# Массив nucList позволит обратиться к определенной клетке словаря не по названию ячейки(char), а по номеру(int) 
# (e.g. dict_nuc_cpg['A']['C'] == dict_nuc_cpg[nucList[0]][nucList[1]], при этом вариант dict_nuc_cpg[0][1] не является рабочим)
nucList = ['A', 'C', 'G', 'T']
dict_nuc_cpg = {'A': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'C': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'G': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'T': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'N': {'N': 0, 'A': 0,   'C': 0,    'G': 0,     'T': 0}}
dict_nuc_ncpg =  {'A': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'C': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'G': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'T': {'A': 0,   'C': 0,    'G': 0,     'T': 0, 'N':0},
            'N': {'N': 0, 'A': 0,   'C': 0,    'G': 0,     'T': 0}}

# В матрицы comparison и comparison_ncpg представляет из себя разность константной и временной матриц
# Временная матрица обновляется каждую итерацию и представляет из себя матрицу вероятностей перехода от нуклеотида Xn к Xn+1, подсчитанную в диапазоне {Xn, Xn+1, Xn+2, ..., Xn+CPG_SIZE}
comparison = np.full((4,4), 0.0)
comparison_ncpg = np.full((4,4), 0.0)

# Создаем одномерные массивы:
# x - массив, представляющий из себя диапазон {0, 1, 2, ... , len(<массив-исходная последовательность нуклеотидов>)}
# yCpG, yNon - массивы, заполняющиеся данными о схожести участка длиной <CPG_SIZE> нуклеотидов на CpG островок или на не CGI участок
# islands_list - массив, заполнется координатами найденного CpG островка, данные в нем подлежат оброботке
# (e.g. если CpG островок по координатам (1000, 1300), то массив islands_list представляет из себя [1000, 1001, 1002, ... ,1300])
# res - двумерный массив, представляющий из себя таблицу из двух столбцов, содержащих координаты островков
# (e.g. [[1000, 1300], [2000, 2300], [5500, 5800]] == res, если CpG островки в последовательности расположены по координатам (1000, 1300), (2000, 2300), (5500, 5800))
x = np.linspace(0, len(example.seq), len(example.seq)) 
yCpG = []
yNon = []
yAdd = []
islands_list = []
counter = 0

# Функция MatrixInit заполняет входную матрицу данными о первых <CPG_SIZE> нуклеотидах 
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
yNon.append(MatrixInit(dict_nuc_ncpg, 1))
for i in range(CPG_SIZE): # заполнение первых <CPG_SIZE> ячеек массива выходных данных для корректного отображения графика
    yAdd.append(0)

def CatchIsland (coord, data, param): # функция, определяющая начало и конец островка после его фиксации (островок фиксируется в середине)
    coord_add = coord
    value = data[coord]
    
    if param == 1:
        value = data[coord]
        while value > 0:
            value = data[coord_add]
            coord_add -= 1
    
    if param == 2:
        while value > 0:
            value = data[coord_add]
            if coord_add+1 == len(data):
                break
            else:
                coord_add += 1

    return coord_add - CPG_SIZE

n2i = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4} # Словарь, обратный по функционалу nucList

def NucToIdx(nuc):
    return n2i[nuc]

prob_cpg = 1
prob_ncpg = 1

for i in range(CPG_SIZE): # заполнение переменных произведения (prob_cpg & prob_ncpg) первыми <CPG_SIZE> нуклеотидами
    if(example.seq[i].upper() != 'N' and example[i + 1].upper() != 'N'):
        prob_cpg_cur = CPG[NucToIdx(example.seq[i].upper())][NucToIdx(example.seq[i + 1].upper())] * 10
        prob_ncpg_cur = NON_CPG[NucToIdx(example.seq[i].upper())][NucToIdx(example.seq[i + 1].upper())] * 10
        prob_cpg *= prob_cpg_cur
        prob_ncpg *= prob_ncpg_cur  


# Основной цикл программы (tqdm позволяет отобразить наглядный progressbar в консоли, отображающий прогресс и скорость выполнения программы)
# Основной идеей скрипта является движение "сканера", размером <CPG_SIZE> нуклеотидов по исходной последовательности, сравнивая участок "сканера" с CGI участком или с nonCGI участком

#for i in tqdm(range(CPG_SIZE, len(example.seq)-1), leave=True): #scaner move (CpG comp) # строка позволяет переключать входные данные на весь файл
for i in tqdm(range(10000000, 10100000), leave=True): 
    if(example.seq[i - CPG_SIZE].upper() != 'N' and example[i - CPG_SIZE + 1].upper() != 'N'):
        prob_cpg_cur = CPG[NucToIdx(example.seq[i - CPG_SIZE].upper())][NucToIdx(example.seq[i - CPG_SIZE + 1].upper())] * 10
        prob_ncpg_cur = NON_CPG[NucToIdx(example.seq[i - CPG_SIZE].upper())][NucToIdx(example.seq[i - CPG_SIZE + 1].upper())] * 10
        prob_cpg /= prob_cpg_cur
        prob_ncpg /= prob_ncpg_cur 
    if(example.seq[i].upper() != 'N' and example[i + 1].upper() != 'N'):
        prob_cpg_cur = CPG[NucToIdx(example.seq[i].upper())][NucToIdx(example.seq[i + 1].upper())] * 10
        prob_ncpg_cur = NON_CPG[NucToIdx(example.seq[i].upper())][NucToIdx(example.seq[i + 1].upper())] * 10
        prob_cpg *= prob_cpg_cur
        prob_ncpg *= prob_ncpg_cur  

    yAdd.append(math.log(prob_cpg/prob_ncpg))

last_call = 2*CPG_SIZE
last_write = ""

# нормализация выходных данных (приведение их к общему для любых данных масштабу)
top_peak = max(yAdd)
for i in range(len(yAdd)):
    yAdd[i] /= top_peak/10
triggeer_point = max(yAdd) - np.mean(yAdd)*1.5
print(triggeer_point)

# запись в файл координат зафиксированных островков
for i in range(1, len(yAdd)):
    if yAdd[i-1] > triggeer_point and i - last_call >= CPG_SIZE:
        now_write = f"{int(CatchIsland(i-1, yAdd, 1))} \t {int(CatchIsland(i-1, yAdd, 2))} \n"
        if last_write != now_write:
            f.write(f"{int(CatchIsland(i-1, yAdd, 1))} \t {int(CatchIsland(i-1, yAdd, 2))} \n")
            last_write = f"{int(CatchIsland(i-1, yAdd, 1))} \t {int(CatchIsland(i-1, yAdd, 2))} \n"
            last_call = i

# закрытие текстового файла, вывод графика
f.close()
fig, ax = plt.subplots()
ax.plot(yAdd)
ax.plot(np.full(len(yAdd), triggeer_point), color = 'gray')
plt.show()