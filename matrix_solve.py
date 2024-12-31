import numpy as np

def is_diagonally_dominant(A):
    n = A.shape[0]
    for i in range(n):
        row_sum = np.sum(np.abs(A[i, :])) - abs(A[i, i])
        if abs(A[i, i]) <= row_sum:
            return False
    return True

def pivot_for_diagonal_dominance(A, b):
    n = A.shape[0]
    A = A.copy()
    b = b.copy()
    for i in range(n):
        max_row = np.argmax(np.abs(A[i:, i])) + i
        if max_row != i:
            A[[i, max_row]] = A[[max_row, i]]
            b[[i, max_row]] = b[[max_row, i]]
    return A, b

def jacobi_iteration(A, b, x_current):
    n = A.shape[0]
    x_new = np.copy(x_current)
    for i in range(n):
        s = np.sum(A[i, :] * x_current) - A[i, i] * x_current[i]
        x_new[i] = (b[i] - s) / A[i, i]
    return x_new

def gauss_seidel_iteration(A, b, x_current):
    n = A.shape[0]
    x_new = np.copy(x_current)
    for i in range(n):
        s = 0
        for j in range(n):
            if j != i:
                s += A[i, j] * x_new[j]
        x_new[i] = (b[i] - s) / A[i, i]
    return x_new

def jacobi_solver(A, b, epsilon=1e-3, max_iter=1000):
    print("=== Jacobi Method ===")
    if is_diagonally_dominant(A):
        print("Matrix is already diagonally dominant.")
        diag_dom = True
    else:
        print("Matrix is not diagonally dominant. Attempting to reorder...")
        A_reordered, b_reordered = pivot_for_diagonal_dominance(A, b)
        if is_diagonally_dominant(A_reordered):
            print("Successfully obtained a diagonally dominant matrix via row swaps.")
            A, b = A_reordered, b_reordered
            diag_dom = True
        else:
            print("Still not diagonally dominant after row swaps.")
            A, b = A_reordered, b_reordered
            diag_dom = False

    x = np.zeros_like(b).flatten()
    for k in range(1, max_iter + 1):
        x_new = jacobi_iteration(A, b, x)
        print(f"Iteration {k}: x = {x_new}")

        if np.linalg.norm(x_new - x, ord=np.inf) < epsilon:
            if diag_dom:
                print(f"Converged in {k} iterations (Jacobi, diagonally dominant).")
            else:
                print(f"Converged in {k} iterations (Jacobi) despite NOT being diagonally dominant.")
            return x_new
        x = x_new

    # If we get here, we did not converge within max_iter
    print("Jacobi method did NOT converge within the maximum number of iterations.")
    return None


def gauss_seidel_solver(A, b, epsilon=1e-3, max_iter=1000):
    print("=== Gauss-Seidel Method ===")
    if is_diagonally_dominant(A):
        print("Matrix is already diagonally dominant.")
        diag_dom = True
    else:
        print("Matrix is not diagonally dominant. Attempting to reorder...")
        A_reordered, b_reordered = pivot_for_diagonal_dominance(A, b)

        if is_diagonally_dominant(A_reordered):
            print("Successfully obtained a diagonally dominant matrix via row swaps.")
            A, b = A_reordered, b_reordered
            diag_dom = True
        else:
            print("Still not diagonally dominant after row swaps.")
            A, b = A_reordered, b_reordered
            diag_dom = False

    # Now run the Gauss-Seidel iteration
    x = np.zeros_like(b).flatten()
    for k in range(1, max_iter + 1):
        x_new = gauss_seidel_iteration(A, b, x)
        print(f"Iteration {k}: x = {x_new}")

        # Check convergence
        if np.linalg.norm(x_new - x, ord=np.inf) < epsilon:
            if diag_dom:
                print(f"Converged in {k} iterations (Gauss-Seidel, diagonally dominant).")
            else:
                print(f"Converged in {k} iterations (Gauss-Seidel) despite NOT being diagonally dominant.")
            return x_new
        x = x_new

    # If we get here, we did not converge within max_iter
    print("Gauss-Seidel method did NOT converge within the maximum number of iterations.")
    return None

if __name__ == "__main__":
    # Example matrix and vector
    A = np.array([
        [4, 2, 0],
        [2, 10, 4],
        [0, 4, 5]
    ], dtype=float)
    b = np.array([2, 6, 5], dtype=float)

    print("Choose the method you want to use:")
    print("1) Jacobi")
    print("2) Gauss-Seidel")

    choice = input("Enter your choice (1 or 2): ").strip()

    if choice == '1':
        x_solution = jacobi_solver(A, b, epsilon=1e-4, max_iter=100)
        if x_solution is not None:
            print("Solution (Jacobi):", x_solution)
    elif choice == '2':
        x_solution = gauss_seidel_solver(A, b, epsilon=1e-4, max_iter=100)
        if x_solution is not None:
            print("Solution (Gauss-Seidel):", x_solution)
    else:
        print("Invalid choice. Please enter 1 for Jacobi or 2 for Gauss-Seidel.")
