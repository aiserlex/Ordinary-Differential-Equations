import numpy as np
import matplotlib.pyplot as plt

X = (0, 5)
Y = (0, 5)

# right side of the first equation in the system
def f(x, y):
    return x - x*y

# right side of the second equation in the system
def g(x, y):
    return -y + x*y

x, y = np.meshgrid(np.linspace(X[0], X[1], 25),
                   np.linspace(Y[0], Y[1], 25))

fig, ax = plt.subplots()
ax.set_title('Vector Field')
ax.set_xlim(X)
ax.set_ylim(Y)

# normalizing factor
hyp = np.hypot(f(x, y), g(x, y))

# plotting the vector field
ax.quiver(x, y, f(x, y)/hyp, g(x, y)/hyp, angles='xy')

plt.show()
