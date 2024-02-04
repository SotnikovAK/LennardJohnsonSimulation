import numpy as np
import matplotlib.pyplot as plt
import math

with open('./DATA/AbsV.txt', 'r') as file:
    A = [float(line) for line in file]

X = {}
accuracy = 0.05
for value in A:
    i = int(value / accuracy)
    if i in X:
        X[i] += 1
    else:
        X[i] = 1

plt.bar(X.keys(), X.values(), width=1, align='center')
plt.xlabel('Промежутки')
plt.ylabel('Количествоо вхождений в промежуток')
plt.show()