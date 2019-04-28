function [eig_vals, eig_vec, n] = power_iter(A, x_0, err_tolerance, max_num_iter)
%% This function calculates the dominant eigenvalue and eigen vectors of a matrix via the power method

% INPUTS: 
%       A - Matrix
%       x_0 - Initial guess for eigen vector 
%       err_tolerance - Error tolerance
%       max_num_iter - Maximum number of iterations

% OUTPUTS:
%       eig_vals - Estimates of dominant eigenvalues
%       eig_vec - Estimates of corresponding normalized eigen vectors.
%       n - The nmber of iterations performed to meet error tolerance


% Power method to calculate approximate eigenvector and eigenvalues
for k=1:max_num_iter
    x=A^k*x_0;
    x=abs(x/min(abs(x))); 
    Ax=A*x; % Dominant eigen vector
    
    Axdotx=dot(Ax,x);
    xdotx=dot(x,x);
    eigval=Axdotx/xdotx; % Dominant eigen value
    eigvalx=eigval*x;
    n=norm(Ax-eigvalx);
    % Check if max iterations have been exceeded
    if n>=e && k==i 
      disp('ERROR!!! ITERATION LIMIT REACHED.')
      break
    end
    % Check if dominant eigen vector/value are within error tolerance
    if n<e
        disp('Iteration')
         fprintf('%1.0f\n',k)
         fprintf('\n')
        disp('Eigen Value')
         fprintf('%4.5f\n',eigval)
         fprintf('\n')
        disp('Eigen Vector')
         fprintf('%4.5f\n',x)
         fprintf('\n')
          break
    end

end