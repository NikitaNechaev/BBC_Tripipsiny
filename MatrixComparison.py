from Bio import SeqIO
import numpy as np, matplotlib.pyplot as plt
from tqdm import tqdm

np.set_printoptions(threshold=15000)

filename = 'SourceFiles\generated_sequence.fa'
#filename = 'BioBootCamp 2022, 2-ой этап\hg38_chr1_and_chr2 copy.fa'
#filename = 'myFasta.fa'
seq = SeqIO.parse(filename, 'fasta')
for rec in seq:
    example = rec
    break

print("\n\n\nАнализируется файл - ", filename, "\n размер - ", len(example.seq), '\n', example, "\n\n")

def TransitionRegister (firstNuc, secondNuc, array):
    if firstNuc == 'A':
        if secondNuc == 'A':
            array[0][0] += 1
        elif secondNuc == 'C':
            array[0][1] += 1
        elif secondNuc == 'G':
            array[0][2] += 1
        elif secondNuc == 'T':
            array[0][3] += 1

    if firstNuc == 'C':
        if secondNuc == 'A':
            array[1][0] += 1
        elif secondNuc == 'C':
            array[1][1] += 1
        elif secondNuc == 'G':
            array[1][2] += 1
        elif secondNuc == 'T':
            array[1][3] += 1

    if firstNuc == 'G':
        if secondNuc == 'A':
            array[2][0] += 1
        elif secondNuc == 'C':
            array[2][1] += 1
        elif secondNuc == 'G':
            array[2][2] += 1
        elif secondNuc == 'T':
            array[2][3] += 1

    if firstNuc == 'T':
        if secondNuc == 'A':
            array[3][0] += 1
        elif secondNuc == 'C':
            array[3][1] += 1
        elif secondNuc == 'G':
            array[3][2] += 1
        elif secondNuc == 'T':
            array[3][3] += 1

CPG_SIZE = 300

NON_CPG = np.array([[0.300, 0.205, 0.285, 0.210], [0.322, 0.298, 0.078, 0.302], [0.248, 0.246, 0.298, 0.208], [0.177, 0.239, 0.292, 0.292]])
CPG = np.array([[0.180, 0.274, 0.426, 0.120], [0.171, 0.368, 0.274, 0.188], [0.161, 0.339, 0.375, 0.125], [0.079, 0.355, 0.384, 0.182]])
#целевые матрицы (МЦ)

dynAr = np.full((4,4),0.0) # динамическая матрица (МД)
compAr = np.full((4,4), 0.0) # матрица сравнений (МС)
# NON_CPG_mult = NON_CPG*300
# CPG_mult = CPG*300

x = []
y = []
yNon = []

f = open('testOutput.txt', 'w')
f.write("MatrixComparison IS WRITING...\n\n\n")

##########################

for p in tqdm(range(len(example.seq)-CPG_SIZE-1)):
    if p == 0:
        for k in range(p, CPG_SIZE+p):      # проверка на островок ПЕРВЫЙ РАЗ
            TransitionRegister(example.seq[k].upper(), example.seq[k+1].upper(), dynAr)

            for i in range(len(dynAr)):
                for j in range(len(dynAr[i])):
                    if k == 0:
                        compAr[i][j] = (CPG[i][j] - dynAr[i][j])**2
                    else:
                        compAr[i][j] = (CPG[i][j] - dynAr[i][j]/k)**2
            #f.write(f"{k} iteration. \n x = {np.sum(compAr)} \n array:\n{dynAr}\n\n")
            #print(p/(len(example.seq)-CPG_SIZE)*100, "%")

            i = 0
            j = 0

            x.append(k)
            y.append(np.sum(compAr))
    
    subArFirst = np.full((4,4), 0)  # инит временных матриц крайних переходов
    subArLast = np.full((4,4), 0)

    TransitionRegister(example.seq[p].upper(), example.seq[p+1].upper(), subArFirst)    # расчет крайних от сканера переходов
    TransitionRegister(example.seq[p+CPG_SIZE].upper(), example.seq[p+CPG_SIZE+1].upper(), subArLast)

    dynAr -= subArFirst # модификация МД с учетом крайних переходов
    dynAr += subArLast

    for i in range(len(dynAr)):     # проверка на островок НЕ в первый раз
        for j in range(len(dynAr[i])):
            compAr[i][j] = (CPG[i][j] - dynAr[i][j]/k)**2

    i = 0
    j = 0

    x.append(p+k)
    y.append(np.sum(compAr))

################################

for p in tqdm(range(len(example.seq)-CPG_SIZE-1)):
    if p == 0:
        for k in range(p, CPG_SIZE+p):      # проверка на островок ПЕРВЫЙ РАЗ
            TransitionRegister(example.seq[k].upper(), example.seq[k+1].upper(), dynAr)

            for i in range(len(dynAr)):
                for j in range(len(dynAr[i])):
                    if k == 0:
                        compAr[i][j] = (NON_CPG[i][j] - dynAr[i][j])**2
                    else:
                        compAr[i][j] = (NON_CPG[i][j] - dynAr[i][j]/k)**2
            #f.write(f"{k} iteration. \n x = {np.sum(compAr)} \n array:\n{dynAr}\n\n")
            #print(p/(len(example.seq)-CPG_SIZE)*100, "%")

            i = 0
            j = 0

            #x.append(k)
            yNon.append(np.sum(compAr))
    
    subArFirst = np.full((4,4), 0)  # инит временных матриц крайних переходов
    subArLast = np.full((4,4), 0)

    TransitionRegister(example.seq[p].upper(), example.seq[p+1].upper(), subArFirst)    # расчет крайних от сканера переходов
    TransitionRegister(example.seq[p+CPG_SIZE].upper(), example.seq[p+CPG_SIZE+1].upper(), subArLast)

    dynAr -= subArFirst # модификация МД с учетом крайних переходов
    dynAr += subArLast

    for i in range(len(dynAr)):     # проверка на островок НЕ в первый раз
        for j in range(len(dynAr[i])):
            compAr[i][j] = (NON_CPG[i][j] - dynAr[i][j]/k)**2

    i = 0
    j = 0

    #x.append(p+k)
    yNon.append(np.sum(compAr))

##################################

print("Полученное максимальное приближение к CpG - ", min(y), ". По координатам - ", y.index(min(y)), ". CpG island - ", y.index(min(y)) - 300, ':', y.index(min(y)))

midAr = np.full(len(x)-200, np.mean(y[200:])) # среднее значение по всему графику

yDiv = []               # массив с комбинированным поиском
for i in range(len(y)):
    yDiv.append(yNon[i] - y[i])

#plt.plot(x[200:], y[200:], color = 'blue')
#plt.plot(x[200:], yNon[200:], color = 'red')
#plt.plot(x[200:],midAr, color = 'green')

del yDiv[0:300]


# peaks = ss.find_peaks(yDiv, height=1, threshold=1, distance=1)

# xN = np.linspace(300, 14999, 14699)

# peak_height = peaks[1]['peak_heights']
# peak_pos = xN[peaks[0]]

# print(peaks)

# fig = plt.figure()
# ax = fig.subplots()

# ax.plot(xN, yDiv)
# ax.scatter(peak_pos, peak_height, color = 'r')

# plt.plot(yDiv[300:])
# plt.plot(yS, yDiv[yS], "o")
# plt.scatter(yS, x[yS], color = 'red')


q = np.linspace(0, 0.1)

st = -6
st = 1
plt.plot(np.power(yDiv, st))#-np.mean(np.power(yDiv, st)))

#f.write(' '.join(str(e) for e in yDiv))
plt.show()