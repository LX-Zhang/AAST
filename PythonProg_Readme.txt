Copyright @ Niloufar Abhari, Simon Fraser University

This python program is the implemention of the least squares method to estimate the branch lengths for a tree-child network
from the branch weights of the trees displayed in the tree-child network.  
The input file contains a linear equation system, where an equation is given in a row.
The output contains the solution if sucessess. 

.con is the equality constraints residuals.

.fun is the objective function value at the optimum (if found).

.message is the status of the solution.

.nit is the number of iterations needed to finish the calculation.

.slack is the values of the slack variables, or the differences between the values of the left and right sides of the constraints.

.status is an integer between 0 and 4 that shows the status of the solution, such as 0 for when the optimal solution has been found.

.success is a Boolean that shows whether the optimal solution has been found.

.x is a NumPy array holding the optimal values of the decision variables.