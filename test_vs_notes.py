from advection import *
import matplotlib.pyplot as plt

# Initial setup
# Setting all the initial conditions

nx = 1000  # number of points
u = 1.0  # velocity

CFL_list = [0.1, 0.5, 0.9]
x_range = [0., 1.]
init_range = [1. / 3., 2. / 3.]

# Creating the initial tophat distribution
thx, tha = tophat(nx, x_range=x_range, init_range=init_range)

########################

# Upwind advection calculation and plotting

# Calculation start
dx = 1. / 64.
T = 1.0
plt.plot(thx, tha, 'k--', label="Initial")
for CFL in CFL_list:
    x, a = solve_advection(dx=dx, x_range=x_range,
                           u=u, CFL=CFL, period=1.0,
                           init_range=init_range,
                           method='upwind')
# Calculation end

# Plotting start
    plt.plot(x, a,
             label=r"CFL = {}".format(round(CFL, 2)))

plt.legend(frameon=False)
plt.xlim(0, 1)
plt.xlabel(r"$x$", fontsize=20)
plt.ylabel(r"$a$", fontsize=20)
plt.title(r"Upwind$:\ T={},\ \Delta x={}$".format(T, round(dx, 2)),
          fontsize=22)
plt.tight_layout()
plt.show()
# Plotting end

########################

# FTCS advection calculation and plotting

# Calculation start
dx = 1. / 64.
T = 0.1
plt.plot(thx, tha, 'k--', label="Initial")
for CFL in CFL_list:
    x, a = solve_advection(dx=dx, x_range=x_range,
                           u=u, CFL=CFL, period=0.1,
                           init_range=init_range,
                           method='ftcs')
# Calculation end

# Plotting start
    plt.plot(x, a,
             label=r"CFL = {}".format(round(CFL, 2)))

plt.legend(frameon=False)
plt.xlim(0, 1)
plt.xlabel(r"$x$", fontsize=20)
plt.ylabel(r"$a$", fontsize=20)
plt.title(r"FTCS$:\ T={},\ \Delta x={}$".format(T, round(dx, 2)),
          fontsize=22)
plt.tight_layout()
plt.show()
# Plotting end

########################

# Implicit in-time advection calculation and plotting

# Calculation start
dx = 1. / 64.
T = 1.0
CFL_list = [0.5, 1.0, 10.0]
plt.plot(thx, tha, 'k--', label="Initial")
for CFL in CFL_list:
    x, a = solve_advection(dx=dx, x_range=x_range,
                           u=u, CFL=CFL, period=1.0,
                           init_range=init_range,
                           method='iit')
# Calculation end

# Plotting start
    plt.plot(x, a,
             label=r"CFL = {}".format(round(CFL, 2)))

plt.legend(frameon=False)
plt.xlim(0, 1)
plt.xlabel(r"$x$", fontsize=20)
plt.ylabel(r"$a$", fontsize=20)
plt.title(r"IIT$:\ T={},\ \Delta x={}$".format(T, round(dx, 2)),
          fontsize=22)
plt.tight_layout()
plt.show()
