The project produced a program named AAST for inferring tree-child networks with the smallest hybridization number (HN) 
that display all the input trees given in Newick format.

Arb_Ntk_Approx_NW.c    the C-code for the AAST

20lfs_10trs_50sets.zip   a collection of 50 datasets which each contains 10 tree with 20 taxa.  
40lfs_10trs_50sets.zip   a collection of 50 datasets which each contains 10 tree with 40 taxa.
50lfs_50trs_50sets.zip   a collection of 50 datasets which each contains 10 tree with 40 taxa.

Equations_Generation.c  A c-code to compute a linear equation system for estimating branche weights for a tree-child network
                          one input file is the list of tree edges with weight. Another file contains the list of the edges 
                          of the tree-child networks.


PythonProg_Readme    It contains a breif explanation of the Python program
sample_ntk2_lfs6_r3   A sample network input file to Equations_Generation.c
sample_trees_wt_ntk2   A smaple tree edge weight input file to Equations_Generation.c
sample_EquationSet_ntk2  A sample output file from Equations_Generation.c.

LeastSquaresMethod   A Python program to solve an optimization problem using the least squares method.
Sample_Ntk2_Eqs    A sample input file to LeastSquaresMethod.pyn
Sample_Ntk2_Wt_Output  A smaple output file to LeastSquaresMethod.pyn
