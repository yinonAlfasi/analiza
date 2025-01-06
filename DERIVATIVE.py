import sympy
import math

def initializeSympyPolynomialData():
    """
    gets a polynomial as input, the polynomial syntax should be sympy polynomial syntax,
    returns the necessary information for the rest of the work with the polynomial.

    @return: tuple: (sympy polynomial, derivative, polynomial lambdify, derivative lambdify, polynomial degree)
    """
    x = sympy.symbols('x')
    # get polynomial from user
    polynomial = sympy.cos((2 * x ** 3) + (5 * x ** 2) - 6) / (2 * math.e ** (-2 * x))
    print("Polynomial", polynomial)
    #polynomial = getSympyPoly()
    derivative = sympy.diff(polynomial, x)
    print("derivative: ",derivative)
    f = sympy.utilities.lambdify(x, polynomial)
    print("f: ",f)
    print("f(5) = ", f(5))
    fTag = sympy.utilities.lambdify(x, derivative)
    print("ftag: ",fTag)
    print("f'(5) = ", fTag(5))
    polynomialDegree = int(sympy.degree(polynomial))
    return polynomial, derivative, f, fTag, polynomialDegree



initializeSympyPolynomialData()