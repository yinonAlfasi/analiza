import sympy
import math


def bisection_method(left_boundary, right_boundary, f, f_prime=None):
    print(f"The bisection method is running on the range: {left_boundary} to {right_boundary}")

    if f_prime is None:
        f_prime = lambda x: 1

    max_iterations = math.ceil(-1 * (math.log(10**-10 / (right_boundary - left_boundary)) / math.log(2)))
    epsilon = 10**-10
    iterations = 0

    while iterations < max_iterations and (right_boundary - left_boundary) > epsilon:
        m = (right_boundary + left_boundary) / 2
        U_left = f(left_boundary) / f_prime(left_boundary)
        U_right = f(right_boundary) / f_prime(right_boundary)
        U_m = f(m) / f_prime(m)

        if U_left * U_m < 0:
            right_boundary = m
        else:
            left_boundary = m

        print(f" in {iterations} iterations: left boundary = {left_boundary} to right boundary = {right_boundary} m = {m}")

        iterations += 1

    if iterations >= max_iterations:
        print("Didn't find roots within the maximum number of iterations.")
    else:
        print(f"Root found: {m} in {iterations} iterations.")

def initializeSympyPolynomialData():
     x = sympy.symbols('x')
     polynomial = (x ** 3) - (x) - (1)
     derivative = sympy.diff(polynomial, x)
     f = sympy.utilities.lambdify(x, polynomial)
     f_prime = sympy.utilities.lambdify(x, derivative)
     return f, f_prime
def find_roots_in_range(f, f_prime, left_boundary, right_boundary):
    step = 0.1
    current_left = left_boundary
    while current_left < right_boundary:
        current_right = current_left + step
        if f(current_left) * f(current_right) < 0:
            bisection_method(current_left, current_right, f, f_prime)
        current_left = current_right


f, f_prime = initializeSympyPolynomialData()
find_roots_in_range(f, f_prime, 1, 2)