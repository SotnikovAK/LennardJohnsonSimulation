import numpy as np
import matplotlib.pyplot as plt

with open('./DATA/Energy.txt', 'r') as file:
    X, Y, Z = [], [], []
    E = []
    for line in file:
        A = [0] * 3
        A[0], A[1], A[2] = line.split()
        A = list(map(float, A))
        X.append(A[0])
        Y.append(A[1])
        Z.append(A[2])
        E.append(A[2] - A[1])


plt.plot(X,Y, label='T')
plt.plot(X,Z, label='U')
plt.plot(X,E, label='Полная энергия')

plt.ylabel('Энергия')
plt.xlabel('t,шаг')



plt.show()
