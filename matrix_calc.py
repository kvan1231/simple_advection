from copy import deepcopy


def invert(matrix):
    """
    Use gauss-jordan elimination to invert a matrix.
    We will append an identity matrix to our input matrix,
    using basic row operations to turn the input matrix into
    row-echelon form, the initial identity matrix will be
    the inverse of our input matrix. This is not a general
    matrix inverter, this only works in our case where our
    matrix is in a nice form. We would need to include functions
    to rearrange the matrix into a nice form for this code to
    be more general.

    INPUTS
    matrix (2d list) - the intial matrix to be inverted

    OUTPUTS
    inv_matrix (2d list) - the inverse of the inital input matrix
    """

    # Determine the dimensions of the initial matrix'
    matrix = deepcopy(matrix)
    rows = len(matrix)
    cols = len(matrix[0])

    # Create the inverse matrix
    inv_matrix = [[0 for x in xrange(rows)] for y in xrange(cols)]

    # Create the identity matrix
    identity = [[0 for x in xrange(rows)] for y in xrange(cols)]

    for i in xrange(0, rows):
        for j in xrange(0, cols):
            if i == j:
                identity[i][j] = 1

    # Append the identity matrix to the right of our input
    for i in xrange(0, rows):
        matrix[i] += identity[i]

    # Loop through doing row operations to turn matrix into
    # row-echelon form
    i = 0
    for j in xrange(0, cols):
        # Divide the row by the current element to make the current
        # element unity
        matrix[i] = [element / matrix[i][j] for element in matrix[i]]

        # Rescale all other rows to make their values 0 below current
        # element
        for row in xrange(0, rows):
            if row != i:
                scaled_row = [matrix[row][j] * element
                              for element in matrix[i]]
                matrix[row] = [matrix[row][column] - scaled_row[column]
                               for column in xrange(0, len(scaled_row))]

        # stop looping if we've reached the end of the matrix
        if i == rows or j == cols:
            break
        i += 1

    # Only want right hand matrix which is the inverse of our input
    for i in xrange(0, rows):
        inv_matrix[i] = matrix[i][cols:len(matrix[i])]

    return inv_matrix


def matrixmult(A, B_list):
    """
    Simple code for matrix multiplication, currently works for multiplying
    nxn matrix with an nx1 matrix.

    INPUTS
    A - the nxn matrix to be multiplied
    B_list - the nx1 matrix to be multiplied

    OUTPUTS
    flat_C - a flattened nxn matrix
    """
    A = deepcopy(A)
    B_list = deepcopy(B_list)
    B = []
    for element in B_list:
        B.append([element])

    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    cols_B = len(B[0])

    if cols_A != rows_B:
        print "Cannot multiply the two matrices. Incorrect dimensions."
        return

    # Create the result matrix
    # Dimensions would be rows_A x cols_B
    C = [[0 for row in range(cols_B)] for col in range(rows_A)]

    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                C[i][j] += A[i][k] * B[k][j]
    flat_C = [item for sublist in C for item in sublist]
    return flat_C
