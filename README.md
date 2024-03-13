## Simulation of the logistic differential equation

In this project, you will find an example on how to solve the logistic differential equation numerically using the 
backward (implicit) Euler method. The method results in a nonlinear optimization problem, which is solved using Ipopt 
(via Casadi). 

### The logistic differential equation
The logistic differential equation is given by:
d/dz f(z) = f(z) (1 - f(z)),
where z is a real number.

This equation has an analytical solution of the form: 
f(z) = exp(z) / (exp(z) + C).

The boundary condition f(0) = 1/2 yields the integration constant C = 1. 
The solution can then be written as f(z) = 1 / (1 + exp(-z)), which is known as the logistic curve.

### Numerical integration
We will solve the equation for z in [0, L], given the boundary condition f(0) = 1/2.

Numerical integration is performed using the backward/implicit Euler method. The solution method has the following steps: 
1. Discretize [0, L] into n intervals of length dz = L/n. Denote z[i] = i * dz and x[i] = f(z[i]) for i=0,...,n. 
2. Create Casadi variables: x[i] for i=0,...,n. 
3. Create Casadi expressions for the boundary condition, g[0] = x[0] - 1/2, and discretized differential equation: g[i] = x[i] - x[i-1] - dz * x[i] * (1 - x[i]) for i=1,...,n.
4. Solve the following feasibility problem using Ipopt: minimize 0 subject to g[i] = 0 for i=0,...,n, where the optimization variables are x[0],...,x[n].
5. Return the solution x[0],...,x[n]
