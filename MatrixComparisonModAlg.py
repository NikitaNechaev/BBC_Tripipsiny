# Импорт необходимых библиотек:
from Bio import SeqIO
from tqdm import tqdm
import numpy as np, matplotlib.pyplot as plt

# Инициализация файлов, открытие/создание выходного файла
filename = 'SourceFiles\generated_sequence.fa'
f = open('res.txt', 'w')

# Записываем данные из исходного fasta файла в переменную example
seq = SeqIO.parse(filename, 'fasta')
for rec in seq:
    example = rec
    break

# Создаем константы, определяющие размер CpG островка, матрицы переходов нуклеотидов, характерные для CpG островков и для не CpG участков
CPG_SIZE = 300
NON_CPG = np.array([[0.300, 0.205, 0.285, 0.210], [0.322, 0.298, 0.078, 0.302], [0.248, 0.246, 0.298, 0.208], [0.177, 0.239, 0.292, 0.292]])
CPG = np.array([[0.180, 0.274, 0.426, 0.120], [0.171, 0.368, 0.274, 0.188], [0.161, 0.339, 0.375, 0.125], [0.079, 0.355, 0.384, 0.182]])

# Словари dict_nuc_cpg и dict_nuc_noncpg позволят записывать количество подсчитанных комбинаций динуклеотидов
# (e.g. последовательность ...ATTC... будет записана как dict_nuc_cpg['A']['T'] += 1, dict_nuc_cpg['T']['T'] += 1, dict_nuc_cpg['T']['C'] += 1)
# Массив nucList позволит обратиться к определенной клетке словаря не по названию ячейки(char), а по номеру(int) 
# (e.g. dict_nuc_cpg['A']['C'] == dict_nuc_cpg[nucList[0]][nucList[1]], при этом вариант dict_nuc_cpg[0][1] не является рабочим)
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
islands_list = []
res = [[]]

# Функция MatrixInit заполняет входную матрицу данными о первых <CPG_SIZE> нуклеотидах 
def MatrixInit (def_matrix, iter_num:int):
    for j in range(CPG_SIZE):
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


# Основной цикл программы (tqdm позволяет отобразить наглядный progressbar в консоли, отображающий прогресс и скорость выполнения программы)
# Основной идеей скрипта является движение "сканера", размером <CPG_SIZE> нуклеотидов по исходной последовательности, сравнивая участок "сканера" с CGI участком или с nonCGI участком
for i in tqdm(range(CPG_SIZE, len(example.seq)-1), leave=True):
    dict_nuc_cpg[example.seq[i-CPG_SIZE]][example.seq[i-CPG_SIZE+1]] -= 1 # вычитаем элемент за сканером, данные о переходе, который сканер уже не затрагивает
    dict_nuc_cpg[example.seq[i]][example.seq[i+1]] += 1 # добавляем данные о новом переходе перед сканером. В совокупности с верхней строкой получаем "движение сканера" по исходной последовательности со скоростью 1нукл./1итер.
    dict_nuc_noncpg[example.seq[i-CPG_SIZE]][example.seq[i-CPG_SIZE+1]] -= 1
    dict_nuc_noncpg[example.seq[i]][example.seq[i+1]] += 1 # аналогичные действия, обеспечивающие "движение сканера", фиксирующего nonCGI участки
    for k in range(4):
        sum_of_row = (dict_nuc_cpg[nucList[k]]['A'] + # переменная sum_of_row (как и sum_of_row_ncpg) является суммой количества подсчитанных переходов из нуклеотида X
                dict_nuc_cpg[nucList[k]]['C'] +
                dict_nuc_cpg[nucList[k]]['G'] + 
                dict_nuc_cpg[nucList[k]]['T'])
        sum_of_row_ncgp = (dict_nuc_noncpg[nucList[k]]['A'] + 
                dict_nuc_noncpg[nucList[k]]['C'] +
                dict_nuc_noncpg[nucList[k]]['G'] + 
                dict_nuc_noncpg[nucList[k]]['T'])
        for l in range(4):
            comparison[k][l] = (CPG[k][l] - (dict_nuc_cpg[nucList[k]][nucList[l]] / sum_of_row))**2 # добавляем в матрицу comparison информацию о различии константной матрицы и полученной матрицы
                # показатель разности матриц - сумма всех элементов матрицы, каждая клетка которой - квадрат разности константной вероятности конкретного перехода и полученной на участке сканера вероятности перехода
            comparison_ncpg[k][l] = (NON_CPG[k][l] - (dict_nuc_noncpg[nucList[k]][nucList[l]] / sum_of_row_ncgp))**2
    yCpG.append(np.sum(comparison)) # Добавляем полученные данные в массив с данными для графика
    yNon.append(np.sum(comparison_ncpg))
    if yCpG[i-1-CPG_SIZE] > yNon[i-1-CPG_SIZE]: # Производится запись данных о полученных координатах в выходной файл
        islands_list.append(i-1-CPG_SIZE)
    elif yCpG[i-1-CPG_SIZE] > yNon[i-1-CPG_SIZE] and yCpG[i-CPG_SIZE] < yNon[i-CPG_SIZE]:
        f.write(islands_list[int(CPG_SIZE/2)], '\t', islands_list[:(int(CPG_SIZE/2))], '\n')

fig, ax = plt.subplots() # Создание и вывод графика
ax.plot(yCpG)
ax.plot(yNon)
plt.show()
f.close()