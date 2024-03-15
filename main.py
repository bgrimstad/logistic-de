"""
Created 13 March 2024
Bjarne Grimstad, bjarne.grimstad@gmail.no
"""

import casadi as ca
import numpy as np
import matplotlib.pyplot as plt


def logistic_curve(u):
    return 1 / (1 + np.exp(-u))


def solve_numerically(L, x0=1 / 2, n=100):
    """
    Solve the logistic differential equation numerically using the backward Euler method
    This leads to a system of nonlinear equations, which we can cast as a non-linear programming problem
    The resulting NLP problem is solved using Ipopt (via Casadi)

    :return: Solution
    """

    # We divide z into n intervals of length dz
    dz = L / n

    # Create lists to hold variables and equations/constraints
    x = list()  # Variables
    g = list()  # Constraints (system of equations)

    # Add variables and equations for all intervals
    for i in range(n + 1):  # Note that we add one to get a total length of L

        # Variable for interval i: x[i] = f(u[i])
        x.append(ca.SX.sym(f'x_{i}'))

        if i == 0:
            # Add left boundary condition as equality constraint: x[0] - x0 = 0 => x[0] = u0
            g.append(x[i] - x0)

        else:  # Cells 1,...,n
            # Discretized differential equations as equality constraint (using backward Euler scheme)
            de = x[i] - x[i-1] - dz * x[i] * (1 - x[i])
            g.append(de)

    # Create variable and constraint vectors (for Casadi)
    x_vec = ca.vertcat(*x)
    g_vec = ca.vertcat(*g)

    # Initial guess on solution
    x_guess = [x0] * len(x)

    # Constraint bounds
    # Equality constraints, g(x) = 0, are implemented as: 0 <= g(x) <= 0
    lbg = [0] * len(g)
    ubg = [0] * len(g)

    # Solve system of equations using Ipopt
    f = 0  # We set the objective function, f, to zero to solve a feasibility problem
    nlp = {'x': x_vec, 'f': f, 'g': g_vec}
    solver_config = {'ipopt.print_level': 0, 'print_time': 0}
    solver = ca.nlpsol('S', 'ipopt', nlp, solver_config)
    result = solver(x0=x_guess, lbg=lbg, ubg=ubg)
    if not solver.stats()['success']:
        print('Solver status:', solver.stats())
        raise Exception('Solve was not successful')

    x_opt = result['x']
    x_opt = x_opt.full().flatten().tolist()  # Convert solution to flat numpy array and then to a list

    return x_opt


if __name__ == '__main__':
    L = 4       # u is in [0, L]
    x0 = 1/2    # Boundary condition f(0) = x0
    n = 100     # Number of cells in discretization

    # Solve numerically
    x_numerical = solve_numerically(L, x0, n)

    # Solve analytically
    z = np.linspace(0, L, n + 1)
    x_analytical = logistic_curve(z)

    # Plot results
    plt.figure()
    plt.plot(z, x_numerical, label='Numerical')
    plt.plot(z, x_analytical, '--r', label='Analytical')
    plt.legend()
    plt.show()

