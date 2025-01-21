# ------------------ SNIPPET 1 ------------------
import math
import numpy as np

#NEEDS REFACTORING
def romberg_integration(func, a, b, n):

    h = b - a
    R = np.zeros((n, n), dtype=float)

    # R[0,0] = 0.5 * h * (f(a) + f(b))
    R[0, 0] = 0.5 * h * (func(a) + func(b))

    for i in range(1, n):
        # Halve the step size
        h /= 2.0

        # Calculate the new points to evaluate
        sum_term = 0.0
        # 2^i intervals, but we only do the odd indices
        for k in range(1, 2 ** i, 2):
            sum_term += func(a + k * h)

        # First column of Romberg table: trapezoid approximations
        R[i, 0] = 0.5 * R[i - 1, 0] + h * sum_term

        # Higher-order corrections
        for j in range(1, i + 1):
            R[i, j] = R[i, j - 1] + (R[i, j - 1] - R[i - 1, j - 1]) / (4 ** j - 1)

    return R[n - 1, n - 1]


def simpsons_rule(f, a, b, n):

    if n % 2 != 0:
        raise ValueError("Number of subintervals (n) must be even for Simpson's Rule.")

    h = (b - a) / n

    # Precompute all function values to avoid repeated calls
    values = [f(a + i * h) for i in range(n + 1)]

    # Apply Simpson's pattern: 4 for odd indices, 2 for even indices (excluding endpoints)
    s = values[0] + values[-1]  # endpoints
    s += sum(4 * values[i] for i in range(1, n, 2))
    s += sum(2 * values[i] for i in range(2, n, 2))

    return s * (h / 3.0)


def trapezoidal_rule(f, a, b, n):

    h = (b - a) / n
    # Precompute function values
    values = [f(a + i * h) for i in range(n + 1)]

    # Trapezoidal sum
    # T = (f(a) + f(b)) / 2 + sum(f(midpoints))
    # Then multiply by step size h
    trape_sum = 0.5 * values[0] + 0.5 * values[-1] + sum(values[1:-1])

    return trape_sum * h


def main():
    # Example function to integrate
    def f(x):
        return math.e ** (x ** 2)

    # Ask user for integration limits
    print("Enter the integration limits [a, b]:")
    a = float(input(" a = "))
    b = float(input(" b = "))

    # Ask user for number of subintervals / iterations
    print("\nChoose the number of subintervals (n). For Simpson's Rule, n must be even.")
    n = int(input(" n = "))

    # Ask user which method to use
    print("\nChoose an integration method:")
    print("  1) Romberg Integration")
    print("  2) Simpson's Rule")
    print("  3) Trapezoidal Rule")
    choice = input("Your choice: ")

    if choice == '1':
        # For Romberg, 'n' here indicates the level of Romberg iterations
        result = romberg_integration(f, a, b, n)
        print(f"\nRomberg Integration result over [{a}, {b}] with {n} iterations: {result}")

    elif choice == '2':
        # Simpson's rule requires n to be even
        if n % 2 != 0:
            print("Error: For Simpson's Rule, n must be even. Please rerun and choose an even n.")
            return
        result = simpsons_rule(f, a, b, n)
        print(f"\nSimpson's Rule result over [{a}, {b}] with {n} subintervals: {result}")

    elif choice == '3':
        # Trapezoidal rule
        result = trapezoidal_rule(f, a, b, n)
        print(f"\nTrapezoidal Rule result over [{a}, {b}] with {n} subintervals: {result}")

    else:
        print("Invalid choice. Please run again and choose 1, 2, or 3.")


if __name__ == '__main__':
    main()


# ------------------ SNIPPET 2 ------------------
def checkIfSquare(mat):
    """
    this function checks if the matrixis square.
    :param mat: matrix - type list
    :return: boolean
    """
    rows = len(mat)
    for i in mat:
        if len(i) != rows:
            return False
    return True


def isDDM(m, n):
    """
     check the if given matrix is Diagonally Dominant Matrix or not.
    :param m: the matrix, type list.
    :param n: size of the matrix (nxn)
    :return: boolean
    """
    # for each row
    for i in range(0, n):

        # for each column, finding sum of each row.
        sum1 = 0
        for j in range(0, n):
            sum1 = sum1 + abs(m[i][j])

        # removing the diagonal element.
        sum1 = sum1 - abs(m[i][i])

        # checking if diagonal element is less than sum of non-diagonal element.
        if (abs(m[i][i]) < sum1):
            return False
    return True


def rowSum(row, n, x):
    """
    caculates the rowws sum
    :param row: a single row from the matrix
    :param n: the row's size
    :param x: the x vector with results
    :return: the sum
    """
    sum1 = 0
    for i in range(n):
        sum1 += row[i] * x[i]
    return sum1


def checkResult(result, last_result, n, epsilon):
    """
    checking if the result is accurate enough
    :param result: the most recent result
    :param last_result: the previous result
    :param n: the size of the result vector
    :return: boolean
    """
    for i in range(n):
        if abs(result[i] - last_result[i]) > epsilon:
            return False
    return True


def Jacobi(mat, b, epsilon=0.000001):  # mat needs to be a list, example: l1 = [[2,3],[4,5]]
    """
    caculating matrix to find vareables vector accourding to yaakobi's algorithem
    :param mat: the matrix
    :param b: the result vector
    :return: the variables vector
    """
    # input check
    n = len(mat)
    if not checkIfSquare(mat):
        return "matrix is not square"
    if len(b) != n:
        return "b is not in the right size"

    # check if Diagonally Dominant Matrix
    if not isDDM(mat, n):
        print("matrix is not Diagonally Dominant")

    # taking a guess: all zeros
    last_result = list()
    for i in range(n):
        last_result.append(0)

    result = last_result.copy()

    print("all results:\nx\t\ty\t  z")
    count = 0
    while True:
        for i in range(n):  # for each variable
            result[i] = b[i] - (rowSum(mat[i], n, last_result) - mat[i][i] * last_result[i])
            result[i] /= mat[i][i]

        print("i = " + str(count) + ": " + str(result))
        count += 1
        if checkResult(result, last_result, n, epsilon):
            return result
        last_result = result.copy()


# ------------------ SNIPPET 3 ------------------
from colors import bcolors


def lagrange_interpolation(x_data, y_data, x):
    """
    Lagrange Interpolation

    Parameters:
    x_data (list): List of x-values for data points.
    y_data (list): List of y-values for data points.
    x (float): The x-value where you want to evaluate the interpolated polynomial.

    Returns:
    float: The interpolated y-value at the given x.
    """
    n = len(x_data)
    result = 0.0

    for i in range(n):
        term = y_data[i]
        for j in range(n):
            if i != j:
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
        result += term

    return result


if __name__ == '__main__':
    x_data = [1, 2, 5]
    y_data = [1, 0, 2]
    x_interpolate = 3  # The x-value where you want to interpolate
    y_interpolate = lagrange_interpolation(x_data, y_data, x_interpolate)
    print(bcolors.OKBLUE, "\nInterpolated value at x =", x_interpolate, "is y =", y_interpolate, bcolors.ENDC)

# ------------------ SNIPPET 4 ------------------
from colors import bcolors


def linearInterpolation(table_points, point):
    p = []
    result = 0
    flag = 1
    for i in range(len(table_points)):
        p.append(table_points[i][0])
    for i in range(len(p) - 1):
        if i <= point <= i + 1:
            x1 = table_points[i][0]
            x2 = table_points[i + 1][0]
            y1 = table_points[i][1]
            y2 = table_points[i + 1][1]
            result = (((y1 - y2) / (x1 - x2)) * point) + ((y2 * x1) - (y1 * x2)) / (x1 - x2)
            print(bcolors.OKGREEN, "\nThe approximation (interpolation) of the point ", point, " is: ", bcolors.ENDC,
                  round(result, 4))
            flag = 0
    if flag:
        x1 = table_points[0][0]
        x2 = table_points[1][0]
        y1 = table_points[0][1]
        y2 = table_points[1][1]
        m = (y1 - y2) / (x1 - x2)
        result = y1 + m * (point - x1)
        print(bcolors.OKGREEN, "\nThe approximation (extrapolation) of the point ", point, " is: ", bcolors.ENDC,
              round(result, 4))


if __name__ == '__main__':
    table_points = [(0, 0), (1, 0.8415), (2, 0.9093), (3, 0.1411), (4, -0.7568), (5, -0.9589), (6, -0.2794)]
    x = 1.28
    print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------\n", bcolors.ENDC)
    print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
    print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x)
    linearInterpolation(table_points, x)
    print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n",
          bcolors.ENDC)

# ------------------ SNIPPET 5 ------------------
from colors import bcolors


def neville(x_data, y_data, x_interpolate):
    n = len(x_data)

    # Initialize the tableau
    tableau = [[0.0] * n for _ in range(n)]

    for i in range(n):
        tableau[i][0] = y_data[i]

    for j in range(1, n):
        for i in range(n - j):
            tableau[i][j] = ((x_interpolate - x_data[i + j]) * tableau[i][j - 1] -
                             (x_interpolate - x_data[i]) * tableau[i + 1][j - 1]) / (x_data[i] - x_data[i + j])

    return tableau[0][n - 1]


if __name__ == '__main__':
    # Example usage:
    x_data = [1, 2, 5, 7]
    y_data = [1, 0, 2, 3]
    x_interpolate = 3

    interpolated_value = neville(x_data, y_data, x_interpolate)
    print(bcolors.OKBLUE, f"\nInterpolated value at x = {x_interpolate} is y = {interpolated_value}", bcolors.ENDC)

# ------------------ SNIPPET 6 ------------------
from colors import bcolors
from matrix_utility import *


def GaussJordanElimination(matrix, vector):
    """
    Function for solving a linear equation using gauss's elimination method
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=A(-1)b
    """
    # Pivoting process
    matrix, vector = RowXchange(matrix, vector)
    # Inverse matrix calculation
    invert = InverseMatrix(matrix, vector)
    return MulMatrixVector(invert, vector)


def UMatrix(matrix, vector):
    """
    :param matrix: Matrix nxn
    :return:Disassembly into a  U matrix
    """
    # result matrix initialized as singularity matrix
    U = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # pivoting process
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            # Finding the M(ij) to reset the organs under the pivot
            elementary[j][i] = -(matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)
    # U matrix is a doubling of elementary matrices that we used to reset organs under the pivot
    U = MultiplyMatrix(U, matrix)
    return U


def LMatrix(matrix, vector):
    """
       :param matrix: Matrix nxn
       :return:Disassembly into a  L matrix
       """
    # Initialize the result matrix
    L = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # pivoting process
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            # Finding the M(ij) to reset the organs under the pivot
            elementary[j][i] = -(matrix[j][i]) / matrix[i][i]
            # L matrix is a doubling of inverse elementary matrices
            L[j][i] = (matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)

    return L


def SolveLU(matrix, vector):
    """
    Function for deconstructing a linear equation by ungrouping LU
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=U(-1)L(-1)b
    """
    matrixU = UMatrix(matrix)
    matrixL = LMatrix(matrix)
    return MultiplyMatrix(InverseMatrix(matrixU), MultiplyMatrix(InverseMatrix(matrixL), vector))


def solveMatrix(matrixA, vectorb):
    detA = Determinant(matrixA, 1)
    print(bcolors.YELLOW, "\nDET(A) = ", detA)

    if detA != 0:
        print("CondA = ", Cond(matrixA, InverseMatrix(matrixA, vectorb)), bcolors.ENDC)
        print(bcolors.OKBLUE, "\nnon-Singular Matrix - Perform GaussJordanElimination", bcolors.ENDC)
        result = GaussJordanElimination(matrixA, vectorb)
        print(np.array(result))
        return result
    else:
        print("Singular Matrix - Perform LU Decomposition\n")
        print("Matrix U: \n")
        print(np.array(UMatrix(matrixA, vectorb)))
        print("\nMatrix L: \n")
        print(np.array(LMatrix(matrixA, vectorb)))
        print("\nMatrix A=LU: \n")
        result = MultiplyMatrix(LMatrix(matrixA, vectorb), UMatrix(matrixA, vectorb))
        print(np.array(result))
        return result


def polynomialInterpolation(table_points, x):
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points]  # Makes the initial matrix

    b = [[point[1]] for point in table_points]

    print(bcolors.OKBLUE, "The matrix obtained from the points: ", bcolors.ENDC, '\n', np.array(matrix))
    print(bcolors.OKBLUE, "\nb vector: ", bcolors.ENDC, '\n', np.array(b))
    matrixSol = solveMatrix(matrix, b)

    result = sum([matrixSol[i][0] * (x ** i) for i in range(len(matrixSol))])
    print(bcolors.OKBLUE, "\nThe polynom:", bcolors.ENDC)
    print('P(X) = ' + '+'.join(['(' + str(matrixSol[i][0]) + ') * x^' + str(i) + ' ' for i in range(len(matrixSol))]))
    print(bcolors.OKGREEN, f"\nThe Result of P(X={x}) is:", bcolors.ENDC)
    print(result)
    return result


if __name__ == '__main__':
    table_points = [(0, 0), (1, 0.8415), (2, 0.9093), (3, 0.1411), (4, -0.7568), (5, -0.9589), (6, -0.2794)]
    x = 1.28
    print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------\n", bcolors.ENDC)
    print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
    print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x, '\n')
    polynomialInterpolation(table_points, x)
    print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n",
          bcolors.ENDC)
