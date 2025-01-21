import sympy


# -----------------------------------------------------------
# 1. Bisection Method
# -----------------------------------------------------------
def bisection_method(polynomial, start_point, end_point, epsilon=1e-10):

    if polynomial(start_point) * polynomial(end_point) > 0:
        print(f"No sign change found between x={start_point} and x={end_point} for the bisection method.")
        return None, None

    left, right = start_point, end_point
    iterations = 0

    # We stop when half the interval is smaller than epsilon
    while (right - left) / 2 > epsilon:
        iterations += 1
        mid = (left + right) / 2.0
        f_left = polynomial(left)
        f_mid = polynomial(mid)

        # Print iteration details
        print(
            f"Iteration {iterations}: "
            f"left={left:.6f}, right={right:.6f}, mid={mid:.6f}, "
            f"f(left)={f_left:.6e}, f(mid)={f_mid:.6e}"
        )

        if f_left * f_mid <= 0:
            right = mid
        else:
            left = mid

    root = (left + right) / 2.0
    return root, iterations

# -----------------------------------------------------------
# 2. Newton-Raphson Method
# -----------------------------------------------------------
def newton_raphson_method(polynomial, derivative, start_point, end_point, epsilon=1e-10, max_iter=1000):

    x_current = (start_point + end_point) / 2.0
    iterations = 0

    for i in range(max_iter):
        iterations += 1
        f_val = polynomial(x_current)
        f_prime_val = derivative(x_current)

        # Print iteration details
        print(
            f"Iteration {iterations}: "
            f"x_current={x_current:.6f}, f(x_current)={f_val:.6e}, "
            f"f'(x_current)={f_prime_val:.6e}"
        )

        if abs(f_prime_val) < 1e-12:
            print(f"Newton-Raphson: derivative is 0 at x={x_current}. No convergence.")
            return None, None

        x_next = x_current - f_val / f_prime_val

        # Check boundary
        if x_next < start_point or x_next > end_point:
            print(f"Newton-Raphson: stepping outside of [{start_point}, {end_point}]. Possibly no convergence.")
            return None, None

        # Check convergence
        if abs(x_next - x_current) < epsilon:
            return x_next, iterations

        x_current = x_next

    print(f"Newton-Raphson: did not converge after {max_iter} iterations.")
    return None, None

# -----------------------------------------------------------
# 3. Secant Method
# -----------------------------------------------------------
def secant_method(polynomial, start_point, end_point, epsilon=1e-10, max_iter=1000):

    x0 = start_point
    x1 = end_point

    f0 = polynomial(x0)
    f1 = polynomial(x1)

    if abs(f0) < epsilon:
        return x0, 0
    if abs(f1) < epsilon:
        return x1, 0

    iterations = 0
    for i in range(max_iter):
        iterations += 1
        # Print iteration details before computing the next approximation
        print(
            f"Iteration {iterations}: "
            f"x0={x0:.6f}, x1={x1:.6f}, "
            f"f(x0)={f0:.6e}, f(x1)={f1:.6e}"
        )

        if abs(f1 - f0) < 1e-12:
            print("Secant Method failed: difference f1 - f0 is too small.")
            return None, None

        x2 = x1 - f1 * ((x1 - x0) / (f1 - f0))
        f2 = polynomial(x2)

        # Check convergence
        if abs(x2 - x1) < epsilon:
            return x2, iterations

        x0, x1 = x1, x2
        f0, f1 = f1, f2

    print(f"Secant Method did not converge after {max_iter} iterations.")
    return None, None


# -----------------------------------------------------------
# 4. Main Program
# -----------------------------------------------------------
def main():
    x = sympy.symbols('x')
    # Define the polynomial: x^4 + x^3 - 3x^2
    polynomial_expr = x**4 + x**3 - 3*x**2

    derivative_expr = sympy.diff(polynomial_expr, x)

    # Convert to lambda functions
    f = sympy.utilities.lambdify(x, polynomial_expr, 'math')
    f_prime = sympy.utilities.lambdify(x, derivative_expr, 'math')

    # Hardcoded range and step
    left_range = -3.0
    right_range = 2.0
    step = 0.1

    # Ask the user which method to use
    print("Choose a method to find the roots:")
    print("1 - Bisection Method")
    print("2 - Newton-Raphson Method")
    print("3 - Secant Method")
    choice = input("Enter the method number (1/2/3): ")

    def find_root_method(a, b):
        if choice == '1':
            print(f"\n----- Bisection Method on [{a}, {b}] -----\n")
            return bisection_method(f, a, b, epsilon=1e-10)
        elif choice == '2':
            print(f"\n----- Newton-Raphson Method on [{a}, {b}] -----\n")
            return newton_raphson_method(f, f_prime, a, b, epsilon=1e-10)
        elif choice == '3':
            print(f"\n----- Secant Method on [{a}, {b}] -----\n")
            return secant_method(f, a, b, epsilon=1e-10)
        else:
            print("Invalid method choice.")
            return None, None

    # Scan the range in subintervals of size 'step'
    current_left = left_range
    found_roots = []

    while current_left < right_range:
        current_right = current_left + step
        if current_right > right_range:
            current_right = right_range

        # 1. Check sign change in f(x)
        if f(current_left) * f(current_right) < 0:
            root, it = find_root_method(current_left, current_right)
            if root is not None:
                found_roots.append((root, it))

        # 2. Check sign change in f'(x) for possible multiple root
        if f_prime(current_left) * f_prime(current_right) < 0:
            print(f"\n----- Checking derivative's sign change in [{current_left}, {current_right}] -----\n")
            root_deriv, it_deriv = bisection_method(f_prime, current_left, current_right, epsilon=1e-10)
            if root_deriv is not None:
                # If f(root_deriv) is close to zero, it might be a multiple root
                if abs(f(root_deriv)) < 1e-3:
                    found_roots.append((root_deriv, it_deriv))

        current_left = current_right
        if current_left >= right_range:
            break

    # Filter out duplicates (roots that are numerically very close)
    unique_roots = []
    for r, i in found_roots:
        if not any(abs(r - ur[0]) < 1e-7 for ur in unique_roots):
            unique_roots.append((r, i))

    # Final report
    if not unique_roots:
        print("\nNo roots were found in the specified range.")
    else:
        print("\n----- Summary of found roots -----")
        for idx, (root_val, iterations) in enumerate(unique_roots, start=1):
            print(f"Root #{idx}: x = {root_val}, iterations = {iterations}")


if __name__ == "__main__":
    main()
