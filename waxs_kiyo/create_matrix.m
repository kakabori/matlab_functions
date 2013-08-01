function [row_label2, result] = create_matrix(matrix, row_vector, row, row_label)
%UNTITLED Summary of this function goes here
%   
%   matrix is the matrix variable to which row_vector will be concatenated.
%   Two dimensional matrix is assumed.
%   
%   row_vector is the input row vector that will be concatenated to the last row
%   of the matrix.
%   
%   row is the value of row label associated with the input row vector.
%
%   row_label is a column vector for the input label for the row of the matrix. 

result = [matrix; row_vector];
row_label2 = [row_label; row];

end

