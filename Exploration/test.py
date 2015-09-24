import matplotlib.pyplot as plt
import numpy as np
x = np.arange(0,10.0, 0.1)
fig, ax = plt.subplots(2,2)
i = 0
j = 0
ax[i,j].plot(x, x)
ax[0,1].plot(x, x**2)
ax[1,0].plot(x, np.sqrt(x))
ax[1,1].plot(x, 1/x)
plt.show()