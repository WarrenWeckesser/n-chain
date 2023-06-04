import numpy as np
from scipy.linalg import eigh
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib import animation
from numpngw import AnimatedPNGWriter


def nmatrix(n):
    k = np.arange(n, 0, -1)
    return np.minimum.outer(k, k)


def nlevp(x, D, B):
    n = len(x) - 1
    omegasq = x[-1]**4
    f = np.zeros_like(x)
    f[:n] = omegasq * B @ np.sin(x[:n]) - D @ np.tan(x[:n])
    return f


# Number of links (i.e. the 'n' of 'n-chain')
n = 25
# Which branch of relative equilibria to follow
k = 6
num_iters = 16000
delta = 0.001
skip = 125

B = nmatrix(n)
D = np.diag(np.arange(n, 0, -1))

# This code assumes that the eigenvalues will be returned in
# increasing order.

evals, evecs = eigh(D, B)
print(f'{evals=}')

omega2 = evals[k]
# The k-th eigenvector gives the initial search direction.
evec = evecs[:, k]
direc = np.concatenate((evec, [0.0]))
direc = np.sign(direc[0]) * direc
direc[-1] = 1e-12  # FIXME: derive a better value.

phi = np.zeros(n + 1)
phi[-1] = np.sqrt(np.sqrt(omega2))

result = [phi]

for i in range(num_iters):
    direc = direc/np.linalg.norm(direc)
    phiguess = phi + delta*direc
    phinext = fsolve(nlevp, phiguess, args=(D, B), xtol=1e-11, factor=0.1)
    phinext = np.sign(phinext[0]) * phinext
    phinext[-1] = abs(phinext[-1])
    result.append(phinext)
    direc = phinext - phi
    phi = phinext

phi = np.array(result)

cosphi = np.cos(phi[:, :n])
sinphi = np.sin(phi[:, :n])

x = np.cumsum(sinphi, axis=1)
x = np.pad(x, ((0, 0), (1, 0)))
y = -np.cumsum(cosphi, axis=1)
y = np.pad(y, ((0, 0), (1, 0)))


print("Generating the animated PNG file:")


def update_line(num, x, y, line):
    """
    Animation "call back" function for each frame.
    """
    if num >= len(x):
        num = len(x) - num - 1
    line.set_data(x[num, :], y[num, :])
    return line,


fig = plt.figure()
ax = fig.gca()
ax.set_aspect('equal')
ax.set_title(f"Whirling mode {k}\nof the {n}-chain")
ax.set_xlim((1.15*x.min(), 1.15*x.max()))
ax.axis('off')

xx = x[::skip]
yy = y[::skip]

# Plot the initial condition. lineplot is reused in the animation.
lineplot, = ax.plot(xx[0, :], yy[0, :], 'g.-', alpha=0.8)
plt.tight_layout()

print('*** Instantiating FuncAnimation')
ani = animation.FuncAnimation(fig, update_line, frames=2*len(xx),
                              init_func=lambda: None,
                              fargs=(xx, yy, lineplot))
print('*** Creating AnimatedPNGWriter')
writer = AnimatedPNGWriter(fps=25)
filename = f'ani{n:03}m{k:03}.png'
print(f'*** Writing {filename}')
ani.save(filename, dpi=75, writer=writer)
print('*** Done')
