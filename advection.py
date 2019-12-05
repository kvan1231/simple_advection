"""
An implementation of linear 1-d advection

The equation for advection is given by:
    a_t = -u a_x

Need to solve the equation using

    a) upwind
        (a_i^n+1 - a_i^n)/dt = -u(a_i^n - a_i-1^n)/dx

    b) FTSC
        (a_i^n+1 - a_i^n)/dt = -u(a_i+1^n - a_i-1^n)/2dx

    c) IIT (implicit in time)
        (a_i^n+1 - a_i^n)/dt = -u(a_i^n+1 - a_i-1^n+1)/dx

Considering a domain [0, 1],
              u = 1,
              periodic boundary conditions
              initial condition of top-hat
                  0 if        x <= 0.4
              a = 1 if 0.4 <= x <= 0.6
                  0 if 0.6 <= x

"""

# import numpy as np
import matrix_calc as mc
from copy import deepcopy


def solve_advection(dx=0.05, x_range=[0., 1.], u=1.0,
                    CFL=0.9, period=1.0, init_range=[0.4, 0.6],
                    method='upwind'):
    """
    Code for solving advection with an initial tophat distribution.
    The code currently only allows for upwind, FTCS and implicit-in-time
    methods in solving advection. If one of these methods is not chosen
    the code does not work. Currently uses an initial tophat distribution
    but other initial distributions could be included.

    INPUTS
    dx (float) - difference between two consective grid points in the
                 x direction
    x_range (list) - the domain of calculation
    u (float) - the velocity of advection
    CFL (float) - the CFL number to be used in the advection calculation
          CFL = dt * u / dx
    period (float) - the period of simulation, period = 1/u
    init_range (list) - the range where the initial distribution is set

    OUTPUTS
    x (list) - the x array representing the boundaries
    a (list) - the amplitudes of our advection solution
    """

    # Initialize the required variables/properties
    x_min = min(x_range)
    x_max = max(x_range)
    nx = int((x_max - x_min) / dx) + 1

    # Ensure that the method input doesn't depend on capitalization
    method = method.lower()

    # Create the grids
    x = [x_min]
    for i in xrange(1, int(nx)):
        x.append(x[i - 1] + dx)

    a = [0.] * len(x)

    # Define the initial variables
    dt = CFL * dx / u
    t = 0.0
    t_max = period  # only go for one period

    # Create the tophat initial conditions
    for index in xrange(len(x)):
        if x[index] >= min(init_range):
            if x[index] <= max(init_range):
                a[index] = 1.0

    # Create the temporary array to store calculations
    a_temp = [0.] * len(a)

    # If we are doing the implicit in time calculation need to create
    # a matrix as well.
    if method == "iit":
        a_matrix = [[0 for row in xrange(nx)] for col in xrange(nx)]

    # Loop for advection calculation
    while t < t_max:

        # upwind and FTCS can be calculated in a similar way
        if method == "upwind" or method == "ftcs":

            # Loop through the cells to calculate
            for i in xrange(nx):
                if method == "upwind":
                    # Upwind method as described by the equation above
                    a_temp[i] = -CFL * (a[i] - a[i - 1])

                elif method == "ftcs":
                    # FTCS method as described by the equation above
                    # Need special condition for the final point to
                    # represent periodic boundary conditions
                    if i == nx - 1:
                        a_temp[i] = -0.5 * CFL * (a[0] - a[i - 1])
                    else:
                        a_temp[i] = -0.5 * CFL * (a[i + 1] - a[i - 1])

            # Update the main array with values
            for i in xrange(len(a)):
                a[i] = a[i] + a_temp[i]

            # Update the system
            a[0] = a[-1]
            t += dt

        # Implict in time calculation is done differently
        elif method == "iit":

            # Populate the matrix used in implict in time calc
            for i in xrange(0, nx):
                a_matrix[i][i] = 1.0 + CFL
                a_matrix[i][i - 1] = -CFL

            # Define the right hand side of the matrix equation
            RHS = a[:]

            # numpy calculation for checking
            # a_temp = np.linsalg.solve(a_matrix, RHS)

            # create a copy of our a_matrix, if not done we encounter errors
            a_matrix_temp = deepcopy(a_matrix)

            # numpy calculation inverting then doing matrix multiplication
            # a_inv = np.linalg.inv(a_matrix)
            # a_temp = np.matmul(a_inv, RHS)

            # Calculating the solution using custom code

            # Invert the matrix
            a_inv = mc.invert(a_matrix_temp)

            # Matrix multiplication to find solution
            a_temp = mc.matrixmult(a_inv, RHS)

            # Apply calculated solution
            a[:] = a_temp[:]

            # Update the system
            a[0] = a[-1]
            t += dt

        else:
            print("Only supports upwind, FTCS and implicit-in-time (IIT)")

    return x, a


def tophat(nx=1000, x_range=[0., 1.],
           init_range=[0.4, 0.6]):
    """
    Simple code for creating a high resolution tophat distribution.
    Creates an array containing the various x gridpoints of interest.
    Populates an array which contains the relevant a values in a tophat
    formation.

    INPUTS
    nx (float) - number of gridpoints, higher means higher resolution
    x_range (list) - the domain where our initial distribution exists
    init_range (list) - the range where out tophat is non-zero

    OUTPUTS
    x (list) - the x values of our tophat distribution in the domain with
               nx values
    a (list) - the a values of our tophat distribution in the domain, should
               be zero everywhere except between our init_range values.

    """
    # Defining initial variables
    x_min = min(x_range)
    x_max = max(x_range)

    init_min = min(init_range)
    init_max = max(init_range)

    dx = (x_max - x_min) / nx
    x = [0.]

    # Loop to create the required x points
    for i in xrange(1, nx):
        x.append(x[i - 1] + dx)

    # Create an empty a array to hold the values
    a = [0.] * nx

    # Create the tophat distribution
    for index in xrange(nx):
        if x[index] >= init_min:
            if x[index] <= init_max:
                a[index] = 1.0

    return x, a
