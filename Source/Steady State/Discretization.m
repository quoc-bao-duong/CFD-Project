function [n]=Discretization(delta_x,delta_y)
%{
INPUT
 A: Matrix Phi with size n x n
 n: Number of segment along 1 row (x-axis) or 1 column (y-axis)
OUTPUT
 delta_x,delta_y
%}
    n=1/delta_x;
end