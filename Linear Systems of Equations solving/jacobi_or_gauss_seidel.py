import numpy as np

#code was refactored on the original code that was provided in the course together with the samples of pseudocode.
#rewriting was made with the use of copilot and some hand made changes to make the code more readable and understandable.
#exxesive comments were removed and the code was made more compact.
# ----------------------- Diagonal Dominance & Pivoting -----------------------
def is_diagonally_dominant(A):
    """
    Checks if matrix A is diagonally dominant.
    A is assumed to be a square numpy array.
    """
    n = A.shape[0]
    for i in range(n):
        row_sum = np.sum(np.abs(A[i, :])) - abs(A[i, i])
        if abs(A[i, i]) < row_sum:
            return False
    return True

def pivot_for_diagonal_dominance(A, b):
    """
    Attempt to reorder rows so that A has stronger diagonal dominance.
    Returns a potentially re-ordered (A, b).
    """
    A = A.copy()
    b = b.copy()
    n = A.shape[0]
    for i in range(n):
        # Find row with the largest absolute pivot in column i
        max_row = np.argmax(np.abs(A[i:, i])) + i
        if max_row != i:
            # swap rows in A
            A[[i, max_row]] = A[[max_row, i]]
            # swap corresponding entries in b
            b[[i, max_row]] = b[[max_row, i]]
    return A, b

# ----------------------- Iterative Methods -----------------------
def jacobi_iteration(A, b, x_current):
    """
    Single Jacobi iteration step:
    Returns x_new based on x_current.
    """
    n = A.shape[0]
    x_new = np.copy(x_current)
    for i in range(n):
        s = np.sum(A[i, :] * x_current) - A[i, i] * x_current[i]
        x_new[i] = (b[i] - s) / A[i, i]
    return x_new

def gauss_seidel_iteration(A, b, x_current):
    """
    Single Gauss-Seidel iteration step:
    Returns x_new based on x_current.
    """
    n = A.shape[0]
    x_new = np.copy(x_current)
    for i in range(n):
        s = 0
        for j in range(n):
            if j != i:
                s += A[i, j] * x_new[j]
        x_new[i] = (b[i] - s) / A[i, i]
    return x_new

def jacobi_solver(A, b, epsilon=1e-4, max_iter=100):
    """
    Solve A x = b via Jacobi Iteration, printing each iteration.
    Checks for diagonal dominance; if not present, attempts reordering.
    """
    print("=== Jacobi Method ===")

    # Check diagonal dominance; try pivoting if needed
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

    # Jacobi iterations
    x = np.zeros_like(b)
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

    print("Jacobi method did NOT converge within the maximum number of iterations.")
    return None

def gauss_seidel_solver(A, b, epsilon=1e-4, max_iter=100):
    """
    Solve A x = b via Gauss-Seidel Iteration, printing each iteration.
    Checks for diagonal dominance; if not present, attempts reordering.
    """
    print("=== Gauss-Seidel Method ===")

    # Check diagonal dominance; try pivoting if needed
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

    # Gauss-Seidel iterations
    x = np.zeros_like(b)
    for k in range(1, max_iter + 1):
        x_new = gauss_seidel_iteration(A, b, x)
        print(f"Iteration {k}: x = {x_new}")

        if np.linalg.norm(x_new - x, ord=np.inf) < epsilon:
            if diag_dom:
                print(f"Converged in {k} iterations (Gauss-Seidel, diagonally dominant).")
            else:
                print(f"Converged in {k} iterations (Gauss-Seidel) despite NOT being diagonally dominant.")
            return x_new
        x = x_new

    print("Gauss-Seidel method did NOT converge within the maximum number of iterations.")
    return None

# ----------------------- Main Program -----------------------
def main():
    # Example system A x = b
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

if __name__ == "__main__":
    main()
