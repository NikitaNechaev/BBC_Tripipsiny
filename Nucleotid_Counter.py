from itertools import count
from re import A
import matplotlib.pyplot as plt, numpy as np
from Bio import SeqIO

filename = 'BioBootCamp 2022, 2-ой этап\generated_sequence.fa'
seq = SeqIO.parse(filename, 'fasta')
for rec in seq:
    example = rec
    break

def counter (start, end, goal):
    cnt = 0
    for i in range(end-start):
        if example[start+i] == goal:
            cnt += 1
    return cnt

#частота появления - example[i:i+10].count(X) / 10

xA = [0]
xG = [0]
xC = [0]
xT = [0]
y = [0]

for i in range(len(example.seq)):
    xA.append(example.seq[i:i+10].count('A'))
    #print(f"from {i} to {i+10} is {example.seq[i:i+10].count('A')} A")
    xG.append(example.seq[i:i+10].count('G'))
    xC.append(example.seq[i:i+10].count('C'))
    xT.append(example.seq[i:i+10].count('T'))
    y.append(i)

fig, ax = plt.subplots()

ax.plot(y, xA)
ax.plot(y, xG)
ax.plot(y, xC)
ax.plot(y, xT)

plt.show()