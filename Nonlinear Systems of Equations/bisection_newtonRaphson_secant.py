import sympy
#code was refactored on the original code that was provided in the course together with the samples of pseudocode.
#rewriting was made with the use of copilot and some hand made changes to make the code more readable and understandable.
# ============================================================
# 1. Bisection Method
# ============================================================
def bisection_method(f, a, b, tol=1e-10, max_iter=1000, verbose=True):

    f_a = f(a)
    f_b = f(b)
    if f_a * f_b > 0:
        if verbose:
            print(f"[Bisection] No sign change found on [{a}, {b}].")
        return None, 0, []

    iter_data = []
    left, right = a, b
    iterations = 0

    for _ in range(max_iter):
        iterations += 1
        mid = 0.5 * (left + right)
        f_left = f(left)
        f_mid = f(mid)

        # Store iteration info
        iter_data.append({
            'iteration': iterations,
            'left': left,
            'right': right,
            'mid': mid,
            'f_left': f_left,
            'f_mid': f_mid
        })

        if verbose:
            print(f"[Bisection] Iter={iterations} "
                  f"left={left:.6f}, right={right:.6f}, "
                  f"mid={mid:.6f}, f(mid)={f_mid:.6e}")

        # Check convergence
        if abs(right - left) < tol:
            break

        # Narrow down the bracket
        if f_left * f_mid <= 0:
            right = mid
        else:
            left = mid

    root = 0.5 * (left + right)
    return root, iterations, iter_data


# ============================================================
# 2. Newton-Raphson Method
# ============================================================
def newton_raphson_method(f, df, a, b, tol=1e-10, max_iter=1000, verbose=True):

    x_current = 0.5 * (a + b)
    iter_data = []
    iterations = 0

    for i in range(max_iter):
        iterations += 1
        f_val = f(x_current)
        df_val = df(x_current)

        iter_data.append({
            'iteration': iterations,
            'x_current': x_current,
            'f': f_val,
            'df': df_val
        })

        if verbose:
            print(f"[Newton] Iter={iterations} x={x_current:.6f}, "
                  f"f(x)={f_val:.6e}, f'(x)={df_val:.6e}")

        if abs(df_val) < 1e-12:
            if verbose:
                print("[Newton] Derivative too small; stopping.")
            return None, iterations, iter_data

        x_next = x_current - f_val / df_val

        # If stepping outside the interval, stop
        if not (a <= x_next <= b):
            if verbose:
                print(f"[Newton] Stepped outside [{a}, {b}]. Stopping.")
            return None, iterations, iter_data

        if abs(x_next - x_current) < tol:
            return x_next, iterations, iter_data

        x_current = x_next

    if verbose:
        print("[Newton] Max iterations reached without convergence.")
    return None, iterations, iter_data


# ============================================================
# 3. Secant Method
# ============================================================
def secant_method(f, a, b, tol=1e-10, max_iter=1000, verbose=True):

    x0, x1 = a, b
    iter_data = []
    iterations = 0

    for i in range(max_iter):
        iterations += 1
        f0 = f(x0)
        f1 = f(x1)

        # Prevent division by zero
        if abs(f1 - f0) < 1e-12:
            if verbose:
                print("[Secant] Denominator too small; stopping.")
            return None, iterations, iter_data

        x_next = x1 - f1 * (x1 - x0) / (f1 - f0)

        iter_data.append({
            'iteration': iterations,
            'x0': x0,
            'x1': x1,
            'f0': f0,
            'f1': f1,
            'x_next': x_next
        })

        if verbose:
            print(f"[Secant] Iter={iterations} x0={x0:.6f}, "
                  f"x1={x1:.6f}, x_next={x_next:.6f}")

        if abs(x_next - x1) < tol:
            return x_next, iterations, iter_data

        x0, x1 = x1, x_next

    if verbose:
        print("[Secant] Max iterations reached without convergence.")
    return None, iterations, iter_data


# ============================================================
# 4. Unified root_finder
# ============================================================
def root_finder(method, f, df=None, a=None, b=None,
                tol=1e-10, max_iter=1000, verbose=True):

    if method == 'bisection':
        return bisection_method(f, a, b, tol, max_iter, verbose)
    elif method == 'newton':
        if df is None:
            raise ValueError("Newton's method requires a derivative function df.")
        return newton_raphson_method(f, df, a, b, tol, max_iter, verbose)
    elif method == 'secant':
        return secant_method(f, a, b, tol, max_iter, verbose)
    else:
        raise ValueError(f"Unknown method: {method}")


# ============================================================
# 5. Helper to find intervals with sign changes
# ============================================================
def find_brackets(f, left, right, step):

    brackets = []
    x_current = left
    # Use a while loop to increment in steps
    while x_current < right:
        x_next = x_current + step
        if x_next > right:
            x_next = right

        f_current = f(x_current)
        f_next = f(x_next)

        if f_current * f_next < 0:
            brackets.append((x_current, x_next))

        x_current = x_next
    return brackets


# ============================================================
# 6. Main usage example
# ============================================================
def main(method_choice='bisection', step=0.1, verbose=True):

    # Define the polynomial: x^4 + x^3 - 3*x^2
    x_sym = sympy.Symbol('x', real=True)
    polynomial_expr = x_sym ** 4 + x_sym ** 3 - 3 * x_sym ** 2
    derivative_expr = polynomial_expr.diff(x_sym)

    # Convert to callable functions with sympy as the backend
    f = sympy.utilities.lambdify(x_sym, polynomial_expr, 'sympy')
    f_prime = sympy.utilities.lambdify(x_sym, derivative_expr, 'sympy')

    # Range to scan for sign changes
    left_range, right_range = -3.0, 2.0

    # 1) Find intervals where sign changes
    brackets = find_brackets(f, left_range, right_range, step)
    found_roots = []

    # 2) Apply the chosen method on each bracket
    for (a, b) in brackets:
        root, iters, data = root_finder(
            method=method_choice,
            f=f,
            df=f_prime,
            a=a,
            b=b,
            tol=1e-10,
            max_iter=1000,
            verbose=verbose
        )
        if root is not None:
            found_roots.append(root)

    # 3) Optionally check the derivative for multiple roots
    d_brackets = find_brackets(f_prime, left_range, right_range, step)
    for (da, db) in d_brackets:
        root_d, iters_d, data_d = bisection_method(
            f_prime, da, db, 1e-10, 1000, verbose=verbose
        )
        if root_d is not None and abs(f(root_d)) < 1e-3:
            # likely a multiple root
            found_roots.append(root_d)

    # 4) Filter out duplicates (roots that are very close)
    unique_roots = []
    for r in found_roots:
        if not any(abs(r - ur) < 1e-7 for ur in unique_roots):
            unique_roots.append(r)

    # 5) Final report
    if not unique_roots:
        print("No roots were found in the specified range.")
    else:
        print("\n----- Summary of Found Roots -----")
        for idx, r in enumerate(unique_roots, 1):
            # Evaluate f at the found root
            val_at_root = f(r)
            print(f"Root #{idx}: x = {r:.10f}, f(x) ~ {val_at_root:.6e}")


# ============================================================
# Script entry point
# ============================================================
if __name__ == "__main__":
    # Change 'bisection' to 'newton' or 'secant' to try different methods.
    # Set verbose=False if you don't want iteration-by-iteration printouts.
    main(method_choice='bisection', step=0.1, verbose=True)
