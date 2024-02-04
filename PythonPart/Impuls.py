import numpy as np
import matplotlib.pyplot as plt

with open('./DATA/P.txt', 'r') as file:
    X,Y = [],[]
    for line in file:
        A = list(map(float, line.split()))
        X.append(A[0])
        Y.append(np.power(A[1], 2) + np.power(A[2], 2) + np.power(A[3], 2))

plt.plot(X,Y)
plt.ylabel('Импульс, P')
plt.xlabel('t,шаг')

plt.show()