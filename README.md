# PinTOpt: Parallel-in-time solvers for linear optimality systems

This implementation contains both the ParaOpt and ParaDiag algorithms to
perform parallel-in-time (PinT) optimal control. Detailed documentation
can be found in the M-files of the methods themselves, as well as in the
thesis which this code accompanies.

## ParaDiag: Solving the all-at-once system with circulant preconditioners
[The `paradiag` function](paradiag.m) uses the ParaDiag method for solving
the all-at-once system arising in a time discretisation of the optimality
system. See [this paper by Wu et al.](https://doi.org/10.1051/cocv/2020012)
for the base method, and the thesis for extensions and improvements.

## ParaOpt: A Parareal variant for boundary value problems
[The `paraopt` function](paraopt.m) uses the ParaOpt method that sub-divides
the time interval into smaller parts. The original terminal-cost method is
described in [this paper by Gander et al.](https://doi.org/10.1137/19M1292291),
and extensions are offered in the thesis. In particular, the function supports
subspace enhancement and ParaDiag-based preconditioners.
