%%
% Nathan Boyd
% ID: 107311566
%
%---------------------------- ---------------------------- MATLAB 501A Project ---------------------------- ----------------------------
% The purpose of this project is to solve two different linear algebra problems using different solutions. 
%------------------------------------------------------------------------------------------------------------------------------------------------

clc;
clear all;
%% Problem 1: Given the set of linear equations and parameters outlined below 

A = [10 -1 2 0; 
        -1 11 -1 3;
        2 -1 10 -4;
        0 3 -1 8]

b = [6 ; 25 ; -11; 15]

init_approx = [0;0;0;0]
err_tolerance = 1.0e-05

% Solve the set of linear equaitons wtih an approximate solution after 4 iterations.
% What is the minimum number of iterations required to meet error tolerance?
% Prove that Jacobi Iteration will converge to the exact solution for any initial approximation.
% Solve using the following methods:

% a) Jacobi Iteration Method
clear x;
fprintf(1,'\n')
disp('---------------------------- Jacobi Iteration Method ----------------------------')
max_num_iter = 4;
disp(['Initial Approximation=[0;0;0;0]. Find solution after ', int2str(max_num_iter), ' Iterations'])
[x,~] = jacobi_iter(A,b, init_approx,err_tolerance, max_num_iter);
fprintf('Approximate Value: \n')
fprintf('%7.5f \n' , x)
fprintf('\n')

disp(['Initial Approximation=[0;0;0;0]. Find the minimum number of iterations'])
max_num_iter = 100;
[x,num_iter] = jacobi_iter(A,b, init_approx,err_tolerance, max_num_iter);
fprintf('Minimum iterations: %7.5f \n' , num_iter)
fprintf('Approximate Value: \n')
fprintf('%7.5f \n' , x)
clear x

% b) Gauss-Seidel Iteration Method
disp('---------------------------- Gauss-Seidel Iteration Method ----------------------------')
max_num_iter = 4;
disp(['Initial Approximation=[0;0;0;0]. Find solution after ', int2str(max_num_iter), ' Iterations'])
[x,~] = gauss_seidel_iter(A,b, init_approx,err_tolerance, max_num_iter);
fprintf('Approximate Value: \n')
fprintf('%7.5f \n' , x)
fprintf('\n')

disp(['Initial Approximation=[0;0;0;0]. Find the minimum number of iterations'])
max_num_iter = 100;
[x,num_iter] = gauss_seidel_iter(A,b, init_approx,err_tolerance, max_num_iter);
fprintf('Minimum iterations: %7.5f \n' , num_iter)
fprintf('Approximate Value: \n')
fprintf('%7.5f \n' , x)
clear x

%% Problem 2: Calculate the dominant eigenvalue and associated eigenvector for the matrix

A = [2, -1, 0;
       -1, 2, -1;
       0, -1, 2 ];

err_tolerance = 1.0e-06;
% b) calculate the 3 eigenvalues and corresponding eigenvectors of the matrix
disp('---------------------------- MATLAB "eig" Method ----------------------------')
disp('Calculate the 3 eigenvalues and eigenvectors of the matrix with the "eig" method in MATLAB')
[V,D] = eig(A)

% c) Use the power method to solve the matrix with the initial guess and tolerance
init_approx = [1;1;1];

% d) Answer (c) except with a new initial approximation
init_approx = [1;0;0];

% e) Why is the rate of convergece more rapid in (c) than (d)