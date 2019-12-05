from advection import *
import matplotlib.pyplot as plt

# Initial setup

# Setting all the initial conditions

nx = 1000
u = 1.0

CFL_list = [0.1, 1. / 3., 0.7]
dx_list = [0.05, 1. / 30., 0.01]
T_list = [1.0, 0.1]
x_range = [0., 1.]
init_range = [0.4, 0.6]
# Creating the initial tophat distribution
thx, tha = tophat(nx, x_range=x_range, init_range=init_range)

plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

########################

# Upwind advection calculation and plotting

# Calculation start
for dx in dx_list:
    for T in T_list:
        plt.figure(figsize=(8, 8))
        plt.clf()
        plt.plot(thx, tha, 'k--', label="Initial")
        for CFL in CFL_list:
            x, a = solve_advection(dx=dx, x_range=x_range,
                                   u=u, CFL=CFL, period=T,
                                   init_range=init_range,
                                   method='upwind')
# Calculation end

# Plotting start
            plt.plot(x, a,
                     label=r"CFL = {}".format(round(CFL, 2)))

        plt.legend(frameon=False)
        plt.xlim(0, 1)
        plt.xlabel(r"$x$", fontsize=28)
        plt.ylabel(r"$a$", fontsize=28)
        plt.title(r"Upwind$:\ T={},\ \Delta x={}$".format(T, round(dx, 2)),
                  fontsize=28)
        plt.tight_layout()
        fig_string = "upwind_T{}_dx{}".format(T, round(dx, 2)).replace(".", "")
        plt.savefig(fig_string + "  .pdf", format='pdf')
# Plotting end

########################

# FTCS advection calculation and plotting

# Calculation start
for dx in dx_list:
    for T in T_list:
        plt.figure(figsize=(8, 8))
        plt.clf()
        plt.plot(thx, tha, 'k--', label="Initial")
        for CFL in CFL_list:
            x, a = solve_advection(dx=dx, x_range=x_range,
                                   u=u, CFL=CFL, period=T,
                                   init_range=init_range,
                                   method='ftcs')
# Calculation end

# Plotting start
            plt.plot(x, a,
                     label=r"CFL = {}".format(round(CFL, 2)))

        plt.legend(frameon=False)
        plt.xlim(0, 1)
        plt.xlabel(r"$x$", fontsize=28)
        plt.ylabel(r"$a$", fontsize=28)
        plt.title(r"FTCS$:\ T={},\ \Delta x={}$".format(T, round(dx, 2)),
                  fontsize=28)
        plt.tight_layout()
        fig_string = "FTCS_T{}_dx{}".format(T, round(dx, 2)).replace(".", "")
        plt.savefig(fig_string + ".pdf", format='pdf')
# Plotting end

########################

# Implicit in-time advection calculation and plotting

# Calculation start
for dx in dx_list:
    for T in T_list:
        plt.figure(figsize=(8, 8))
        plt.clf()
        plt.plot(thx, tha, 'k--', label="Initial")
        for CFL in CFL_list:
            x, a = solve_advection(dx=dx, x_range=x_range,
                                   u=u, CFL=CFL, period=T,
                                   init_range=init_range,
                                   method='iit')
# Calculation end

# Plotting start
            plt.plot(x, a,
                     label=r"CFL = {}".format(round(CFL, 2)))

        plt.legend(frameon=False)
        plt.xlim(0, 1)
        plt.xlabel(r"$x$", fontsize=28)
        plt.ylabel(r"$a$", fontsize=28)
        plt.title(r"IIT$:\ T={},\ \Delta x={}$".format(T, round(dx, 2)),
                  fontsize=28)
        plt.tight_layout()
        fig_string = "IIT_T{}_dx{}".format(T, round(dx, 2)).replace(".", "")
        plt.savefig(fig_string + ".pdf", format='pdf')
# Plotting end
