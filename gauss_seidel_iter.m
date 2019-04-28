function  [x, num_iter] = gauss_seidel_iter(A, b, x_0, err_tolerance, max_num_iter)
%% This function solves a given system of linear equations using the Gaus-Seidel iteration method

% INPUTS: 
%       A - Diagonally dominant matrix
%       b - Answer matrix for linear set of equations
%       x_0 - Initial approximation for the solver
%       err_tolerance - The minimum error allowable for solution
%       max_num_iter - The maximum number of iterations allowable

% OUTPUTS
%       x - The solution for the linear set of equations
%       num_iter - The number of iterations required to reach solution within tolerance
num_iter = 0;
n = length(b);
x = zeros(n,1);
error_val = Inf;
prev_x = x_0;
% Create an initial approximate solution
% Gauss-Seidel Method: xnew_i = (1/a_ii)*( b_i - (sum_[j=1->i-1](a_ij*xnew_j)) - (sum_[j=i+1->n](a_ij*xold_j)) )
while error_val > err_tolerance
    prev_x=x;
    for i=1:n
        s=0;
        % solve forward substitution based on traingular forms
        for j=1:i-1
                s=s+A(i,j)*x(j);
        end
        for j=i+1:n
                s=s+A(i,j)*prev_x(j);
        end
        x(i)=(1/A(i,i))*(b(i)-s);
    end
    % check max condition
    if num_iter > max_num_iter
        disp(['ERROR : A solution within tolerance could not be found after ', int2str(max_num_iter),' iterations'])
        break;
    end
    num_iter=num_iter+1;
    error_val = norm(A*x-b);
end