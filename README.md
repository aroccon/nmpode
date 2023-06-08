# Welcome to nmpode!

Git repository of the course "Numerical methods for ODEs and PDEs".
University of Udine, 2023.

## Handbook with notes

Currently under development :-)

## List of exercises (with source code)

There are two folders: src_ode and src_pde. The first contains the source code of the exercises relative to the ODEs (P1 prefix) while the second contains the exercises relative to the PDEs (P2 prefix).

## Brief description of each code script/function.

### ODEs part (P1) - folder src_ode.
- P1_euler_integration.m: Integration of first order equation with Explicit Eulero.
- P1_derivative_check.m:  Error evaluation on the derivative using Eulero explicit method.
- P1_error.m: Convergence of Eulero.
- P1_error_midpoint.m: Convergence of midpoint method.
- P1_mass_spring.m: Mass-spring system with Eulero (divergent result).
- P1_mass_spring_call_to_eulersolver.m: Same as above but with call to function for time integration.
- P1_mass_spring_call_to_midpointsolver.m: Mass-spring system with midpoint method (function is used).
- P1_mass_spring_call_to_rksolver.m: Mass-spring system with Runge-Kutta method (function is used).
- 5 functions (called by the other subroutines): P1_derivs.m, P1_deriv1st.m, P1_rksolver.m, P1_midpointsolver.m, P1_eulersolver.m

### PDEs part (P2) - folder src_pde.

- P2_jacobi.m: Jacobi iterative method for Laplace equation in 2D (fixed aspect ratio=1 of the domain).
- P2_jacobi_extended.m: Same as above but support for different Lx, Ly and grid resolutions along the two directions.
- P2_gauss_seidel.m: Gauss-Siedel iterative method for Laplace equation in 2D (fixed aspect ratio=1 of the domain).
- P2_black_red_gauss_seidel.m: Red-Black  Gauss-Siedel iterative method for Laplace equation in 2D (fixed aspect ratio=1 of the domain).
- P2_poisson_gauss_seidel.m:  Gauss-Siedel iterative method for Poisson equation in 2D (fixed aspect ratio=1 of the domain).
- P2_poisson_drop_gauss_seidel.m: Compute the pressure field obtained from a static droplet (Poisson equation in 2D, derived from P2_poisson_gauss_seidel.m).
- P2_cavity.m: Solves the incompressible Navier-Stokes equations in a rectangular domain with prescribed velocities along the boundary. Equations are solved using finite difference and wit a Projection-correction method (Credit: Benjamin Seibold).
- P2_quasi2d.m: Pseudo-spectral method for the solution of Navier-Stokes equations in a 2D box, external forcing is used to obtain an artificial turbulent-like flow.

### Final 

Feel free to modify and use the codes as you wish. Please let me know if you find any bugs.
