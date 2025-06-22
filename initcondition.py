import numpy as np
import matplotlib.pyplot as plt

A = np.linspace((1,2), (3,4), 21)
B = np.linspace(0,1,21)
plt.plot(B,A)
plt.show()