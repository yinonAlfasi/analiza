import numpy as np


def lagrange_interpolation(x_data, y_data, x):
    # Filter out NaN values
    valid_indices = [i for i, y in enumerate(y_data) if not np.isnan(y)]
    x_data = [x_data[i] for i in valid_indices]
    y_data = [y_data[i] for i in valid_indices]

    n = len(x_data)
    if n == 0:  # Handle the case where all y_data are NaN
        return np.nan

    result = 0.0
    for i in range(n):
        term = y_data[i]
        for j in range(n):
            if i != j:
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
        result += term
    return result


def neville(x_data, y_data, x):
    # Filter out NaN values
    valid_indices = [i for i, y in enumerate(y_data) if not np.isnan(y)]
    x_data = [x_data[i] for i in valid_indices]
    y_data = [y_data[i] for i in valid_indices]

    n = len(x_data)
    if n == 0:
        return np.nan

    tableau = [[0.0] * n for _ in range(n)]
    for i in range(n):
        tableau[i][0] = y_data[i]

    for j in range(1, n):
        for i in range(n - j):
            numerator = ((x - x_data[i + j]) * tableau[i][j - 1] -
                         (x - x_data[i]) * tableau[i + 1][j - 1])
            denominator = (x_data[i] - x_data[i + j])
            tableau[i][j] = numerator / denominator

    return tableau[0][n - 1]


def linear_interpolation(points, x):
    # Filter out points with NaN y-values
    points = [(xi, yi) for xi, yi in points if not np.isnan(yi)]
    points.sort(key=lambda p: p[0])

    if not points:  # Handle the case where all y values are NaN
        return np.nan

    for i in range(len(points) - 1):
        x1, y1 = points[i]
        x2, y2 = points[i + 1]
        if x1 <= x <= x2:
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)

    # Extrapolation
    if x < points[0][0]:
        x1, y1 = points[0]
        x2, y2 = points[1]
    else:
        x1, y1 = points[-2]
        x2, y2 = points[-1]
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def polynomial_interpolation(points, x):
    # Filter out points with NaN y-values
    points = [(xi, yi) for xi, yi in points if not np.isnan(yi)]
    n = len(points)
    if n == 0:
        return np.nan
    X = np.zeros((n, n))
    Y = np.zeros(n)

    for i, (xi, yi) in enumerate(points):
        for j in range(n):
            X[i, j] = xi ** j
        Y[i] = yi

    try:
        coeffs = np.linalg.solve(X, Y)
    except np.linalg.LinAlgError:
        print("Matrix is singular. Polynomial interpolation may not be possible.")
        return None
    return sum(coeffs[j] * x ** j for j in range(n))


def cubic_spline_interpolation(x_data, y_data, x):
    # Filter out NaN values
    valid_indices = [i for i, y in enumerate(y_data) if not np.isnan(y)]
    x_data = [x_data[i] for i in valid_indices]
    y_data = [y_data[i] for i in valid_indices]

    n = len(x_data)
    if n < 2:  # Need at least two points for cubic spline
        return np.nan

    h = [x_data[i + 1] - x_data[i] for i in range(n - 1)]

    # Compute alpha for the system
    alpha = [0] * (n - 1)
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (y_data[i + 1] - y_data[i]) - (3 / h[i - 1]) * (y_data[i] - y_data[i - 1])

    # Tridiagonal system arrays
    l = [1] * n
    mu = [0] * n
    z = [0] * n

    for i in range(1, n - 1):
        l[i] = 2 * (x_data[i + 1] - x_data[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    # Back-substitution
    c = [0] * n
    b = [0] * (n - 1)
    d = [0] * (n - 1)

    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = ((y_data[j + 1] - y_data[j]) / h[j]) - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    # Find interval for x
    i = 0
    if x < x_data[0]:
        i = 0
    elif x > x_data[-1]:
        i = n - 2
    else:
        for j in range(n - 1):
            if x_data[j] <= x <= x_data[j + 1]:
                i = j
                break

    dx = x - x_data[i]
    # Spline polynomial
    return y_data[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3


def main():
    x_data = [0, 1, 2, 3, 4]
    y_data = [0, 0.84, 1, 0.14, -0.76]

    points = list(zip(x_data, y_data))

    while True:
        try:
            x_test = float(input("Enter the x value to interpolate: "))
            if x_test in x_data:
                print("x value already exists in the data. Please enter a different value to interpolate.")
            else:
                break
        except ValueError:
            print("Invalid input for x value. Please enter a number.")

    while True:
        print("\nChoose an interpolation method:\n"
              "1) Lagrange\n"
              "2) Neville\n"
              "3) Linear\n"
              "4) Polynomial\n"
              "5) Cubic Spline\n"
              "0) Exit")
        choice = input("Enter choice [0-5]: ")

        if choice == '0':
            print("Exiting...")
            break
        elif choice == '1':
            result = lagrange_interpolation(x_data, y_data, x_test)
            print(f"[Lagrange]    Value at x={x_test} is {result}")
        elif choice == '2':
            result = neville(x_data, y_data, x_test)
            print(f"[Neville]     Value at x={x_test} is {result}")
        elif choice == '3':
            result = linear_interpolation(points, x_test)
            print(f"[Linear]      Value at x={x_test} is {result}")
        elif choice == '4':
            result = polynomial_interpolation(points, x_test)
            if result is not None:
                print(f"[Polynomial]  Value at x={x_test} is {result}")
        elif choice == '5':
            result = cubic_spline_interpolation(x_data, y_data, x_test)
            print(f"[Cubic Spline] Value at x={x_test} is {result}")
        else:
            print("Invalid choice. Please enter a number between 0 and 5.")


if __name__ == "__main__":
    main()