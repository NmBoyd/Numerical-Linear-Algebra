function  [x, num_iter] = jaboi_iter(A, b, x_0, err_tolerance, max_num_iter)
%% This function solves a given system of linear equations using the Jacobi iteration method

% INPUTS: 
%       A - Diagonally dominant matrix
%       b - Answer matrix for linear set of equations
%       init_approx - Initial approximation for the solver
%       err_tolerance - The minimum error allowable for solution
%       max_num_iter - The maximum number of iterations allowable

% OUTPUTS
%       x - The solution for the linear set of equations
%       num_iter - The number of iterations required to reach solution within tolerance
num_iter = 0;
n = length(b);
x = zeros(n,1);
prev_x = x_0;
error_val = Inf;
% Jacobi Method: x_i = (1/a_ii)*(b_i-sum_[i!=j](a_ij*x_ij))
% Create an initial approximate solution
while error_val>err_tolerance
    prev_x=x;
    for i=1:n
        % find sum of A approx
        s=0;
        for j=1:n      
            if j~=i
                s=s+A(i,j)*prev_x(j);
            end
        end
        % Solve for x by implementing new a approximate
        x(i)=(1/A(i,i))*(b(i)-s);
    end
    num_iter=num_iter+1;
    if num_iter > max_num_iter
        disp(['ERROR : A solution within tolerance could not be found after ', int2str(max_num_iter),' iterations'])
        break;
    end
    error_val = norm(A*x-b);
end

% x =  b+(eye(n)-A)*init_approx;
% curr_x = x;
% % start at 1 since there was already an initial guess made
% i=1;
% % loop until max number of itertions
% while i < max_num_iter 
%     
%     % Substitute the new approximate values into the a matrix solution
%     x =  b+(eye(n)-A)*curr_x;
%     
%     prev_x = curr_x;
%     curr_x = x;
%     
%     num_iter = i;
%     % Check convergence tolerance and break out if met
%     if norm(curr_x-prev_x,1) < err_tolerance 
%         disp(['Found an approximate answer within toelrance using Jacobi Iteration after ', int2str(num_iter), ' Iterations'])
%         break;
%     end
%     i = i+1;
% end    
% if norm(curr_x-prev_x,1) > err_tolerance
%     disp([' ERROR: A value was not found within the error tolerance after ',int2str(max_num_iter),' iterations'])
% end

