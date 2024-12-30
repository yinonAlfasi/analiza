import numpy as np

def pivoting(A):
    n = A.shape[0]
    for i in range(n):
        max_row = np.argmax(np.abs(A[i:n, i])) + i
        if max_row != i:
            A[[i, max_row]] = A[[max_row, i]]
    return A

def check_diagonal_dominance(A):
    n = A.shape[0]
    for i in range(n):
        row_sum = np.sum(np.abs(A[i, :])) - np.abs(A[i, i])
        if np.abs(A[i, i]) <= row_sum:
            return False
    return True

def jacobi_iteration(A, b, X_r, epsilon):
    n = A.shape[0]
    X_r1 = np.copy(X_r)
    for i in range(n):
        sum_ = 0
        for j in range(n):
            if i != j:
                sum_ += A[i, j] * X_r[j]
        X_r1[i] = (b[i] - sum_) / A[i, i]
    return X_r1

def gauss_seidel_iteration(A, b, X_r, epsilon):
    n = A.shape[0]
    X_r1 = np.copy(X_r)
    for i in range(n):
        sum_ = 0
        for j in range(n):
            if i != j:
                sum_ += A[i, j] * X_r1[j]
        X_r1[i] = (b[i] - sum_) / A[i, i]
    return X_r1

def solve_system(A, b, epsilon=0.001, max_iter=1000, method='jacobi'):
    n = A.shape[0]
    X_r = np.zeros(n)
    X_r1 = np.ones(n)
    num_of_runs = 0

    if method == 1:
        iteration_func = jacobi_iteration
    elif method == 2:
        iteration_func = gauss_seidel_iteration
    else:
        raise ValueError("Method must be 1 (jacobi) or 2 (gauss_seidel)")

    while np.linalg.norm(X_r1 - X_r, ord=np.inf) > epsilon and num_of_runs < max_iter:
        num_of_runs += 1
        X_r = np.copy(X_r1)
        X_r1 = iteration_func(A, b, X_r, epsilon)
        print(f"Iteration {num_of_runs}: X_r = {X_r1}")

    if np.linalg.norm(X_r1 - X_r, ord=np.inf) <= epsilon:
        print(f"Converged solution after {num_of_runs} iterations: {X_r1}")
        return X_r1
    else:
        print(f"Did not converge after {num_of_runs} iterations.")
        return None

A = np.array([[1, 2, 4],
              [1,3,3],
              [5,2,1]], dtype=float)

b = np.array([15, 10, 10], dtype=float)

A = pivoting(A)

if check_diagonal_dominance(A):
    print("Matrix is diagonally dominant.")
else:
    print("Matrix is not diagonally dominant.")
    print("not diagonal")

method_choice = int(input("Choose method (1 for Jacobi, 2 for Gauss-Seidel): "))

solution = solve_system(A, b, method=method_choice)
if solution is not None:
    print(f"Solved using method {method_choice}.")