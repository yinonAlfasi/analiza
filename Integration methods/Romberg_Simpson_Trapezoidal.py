import math
import numpy as np

#code was refactored on the original code that was provided in the course together with the samples of pseudocode.
#rewriting was made with the use of copilot and some hand made changes to make the code more readable and understandable.
#exxesive comments were removed and the code was made more compact.
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
