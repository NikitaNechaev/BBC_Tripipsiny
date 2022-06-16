from itertools import count
from re import A
import matplotlib.pyplot as plt, numpy as np
from Bio import SeqIO

filename = 'SourceFiles\generated_sequence.fa'
seq = SeqIO.parse(filename, 'fasta')
for rec in seq:
    example = rec
    break

#частота появления - example[i:i+10].count(X) / 10

xA = [0]
xG = [0]
xC = [0]
xT = [0]
y = [0]

for i in range(len(example.seq)-20):
    xA.append(example.seq[i:i+20].count('CG'))
    y.append(i)

fig, ax = plt.subplots()

ax.plot(y, xA)

plt.show()