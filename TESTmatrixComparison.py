#сравнение матриц

import random
import numpy as np, matplotlib.pyplot as plt

a = np.array([[64,64,64], [64,64,64], [64,64,64]]) # целевая матрица
b = np.array([[0,0,0],[0,0,0],[0,0,0]])   # динамическая матрица
comp = np.full((3,3), 0)    # матрица сравнений

y = []
x = []

f = open('testOutput.txt', 'w')
f.write("FUNCTEST IS WRITING...\n\n\n")

for k in range(1000):   #количество итераций == кол-во нуклеотидов
    b += random.randint(-10, 10)    #обновление значений массива

    for i in range(len(b)):         #заполнение матрицы сравнений(МС) (по осям i и j)
        for j in range(len(b[i])):
            comp[i][j] = (a[i][j] - b[i][j])**2
    f.write(f"{k} iteration. \n x = {np.sum(comp)} \n array:\n{b}\n\n")
    i = 0       #сброс значений координат после одного заполнения МС
    j = 0
    
    y.append(np.sum(comp)) #перевод МС в одну переменную и добавление в массив
    x.append(k)

#x.pop(0) #костыльно убираем нули со времен создания массива 
#y.pop(0)

print(min(y), y.index(min(y))) # вывод координат максимального приближения

plt.plot(x, y) # график
plt.show()
f.close()