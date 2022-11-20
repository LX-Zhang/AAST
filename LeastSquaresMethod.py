# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 18:45:42 2020

@author: niloofar
"""
from scipy.optimize import linprog
import numpy as np
from numpy import savetxt

def read_coeff(n , filename):
    #n is the number of equations
    with open(filename,'r') as fp:
        content = fp.read()
    #removing spaces
    c=content.replace(' ','')
    #putting every element of the long string as an element of an array
    s = np.array([char for char in c])
    #divide the long array into number_of_rows=72 chunks
    l = np.array_split(np.array(s),n)
    A=[]
    rhs=[]
    for item in l:
        #stick 4 last strings of each row to form rhs of each eq
        b = item[-5]+item[-4]+item[-3]+item[-2]
        print(item,'-----',b)
        #change str to float and add it to an array called rhs
        rhs.append(float(b))
        #add all coeffs except the last 5 since they are
        #rhs value (4 chars e.g. 0.90) and '\n' as one char
        temp=np.array(item[:-5])
        #change the type to float
        A.append(temp.astype(np.float))
    #save the processed data as a csv file
    #savetxt(filename+'_csv_output.csv', A, delimiter=',')
    return A,rhs
    

def boundary(number_of_var):
    #set the boundaries of all variables as tuples
    bnd=[]
    for i in range(number_of_var):
        bnd.append((0.01,float("inf")))
    return bnd
        
            
def objective_func(number_of_var):
    #set the vector of objective function coefficients to -1
    #because the linprog slover only solves minimization problem
    return -np.ones(number_of_var)

def save_file_astxtx(file0,filename):
    file2write = open(filename,'a')
    file2write.write(file0)
    file2write.close()
        
#NOTE: you need to add an empty line at the end of the rows in your input file
#NOTE: Specify the number of equations in the input file in the number_of_eqs variable    

#54 is the number of equations in  'equation_dec21.txt'
#to run for 'equation_dec21.txt' uncomment the following two lines
number_of_eqs = 54
lhs_ineq ,rhs_ineq = read_coeff(number_of_eqs,'equation_dec21.txt')

#88 is the number of equations in  'eq_ntk3.txt'
#to run for 'equation_dec21.txt' uncomment the following two lines
#number_of_eqs = 88
#lhs_ineq ,rhs_ineq = read_coeff(number_of_eqs,'eq_ntk3.txt')
print(np.array(lhs_ineq).shape)

var = len(lhs_ineq[0])

    
obj = objective_func(var)
bnd = boundary(var)

opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq, bounds=bnd,
              method="revised simplex")

print(opt)

#to save the output uncomment the following line
#save_file_astxtx(str(opt), 'equation_dec21_weights.txt')

