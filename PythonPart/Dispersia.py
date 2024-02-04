import numpy as np
import matplotlib.pyplot as plt

with open('./DATA/Diffusion.txt', 'r') as file:
    X,Y = [],[]
    for line in file:
        A = list(map(float, line.split()))
        X.append(A[0])
        Y.append(A[1])

plt.plot(X,Y)

plt.show()