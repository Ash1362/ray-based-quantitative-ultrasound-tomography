function [x, res_norm, residual_norm, solve_time] = solveTransmissionLinearSystem(...
    system_matrix, b, x, ls_itetation_number, step_length, method)
%SOLVETRANSMISSIONLINEARSYSTEM solves a linearised system of equations
%associated with Ultrasound transmission imaging
%
%
% DESCRIPTION:
%     solveTransmissionLinearSystem solves a linearised system of equations
%     Ax = b associated with Ultrasound transmission imaging

%
% USAGE:
%
%
% INPUTS:
%     system_matrix        - the system matrix for ultrasound transmission imaging
%     b                    - the right-hand-side of the equation Ax = b
%     x                    - the intial guess
%     ls_itetation_number  - a fixed number of iterations
%     step_length          - the step size applied on the updates of the
%                            sound speed
% OUTPUTS:
%     x                    - solution
%     res_norm             - the L2 norm of res updates for the cg algorithm
%                            and the L2 norm of update directions of the sound speed
%                            for the steepest descent algorithm at each inner iteration
%     residual_norm        - the first component of res_norm
%     solve_time           - the run time
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian 


ts = tic;

% determine whether the approach is SART or not.
if isstruct(step_length)  &&  ~strcmp(method, 'sart')
    error('If the step length is given as a struct, the used approach must be set SART.')
end
    

% allocate a vector for norm of res, which is the residual of linearised 
% subproblem
res_norm = zeros(1, ls_itetation_number);

switch method
    
    case 'steepest_descent'
    
    
    for iter = 1:ls_itetation_number
        
        % update the direction
        update_direction = system_matrix' * (b - system_matrix * x);
        
        % get the norm of the update direction
        res_norm(iter) = norm(update_direction);
        
        % update the solution
        x = x + step_length * update_direction;
        disp(['The number of inner iteration is:' num2str(iter)]);
    end
    


    case 'conjugate_gradient'
        
        % update the residual
         res = system_matrix' * (b - system_matrix * x);
         
         % update the update direction
         update_direction = res;
         
        for iter = 1:ls_itetation_number
            
            
            % update the action of the approximate Hessian on the update direction
            z = system_matrix' * (system_matrix * update_direction);
            
            % update the step length
            alpha = (res' * res) / (update_direction' * z);
            
             % store the current residual
            res_previous = res;
            
            % update the residual
            res = res_previous - alpha * z;
            
            % compute the L2 norm of the residual
            res_norm(iter) = norm(res);
            
            
            % terminate the cg algorithm, if the norm of residual increases
            if iter > 1
                if res_norm(iter) > res_norm(iter-1)
                    break
                end
            end
            
            % update the solution
            x = x + step_length * alpha * update_direction;
            
            % display the L2 norm of the residual
            disp( ['The norm of cg residual is:' num2str(norm(res), '%1.3f')]);
            
             % update the step length for residual
            beta = (res' * res)/(res_previous' * res_previous);
            
            % update the update direction
            update_direction = res + beta * update_direction;
            
            % display the number of iteration
            disp(['The number of inner iteration is:' num2str(iter)]);
            
        end
        
    case 'sart'
        
        for iter = 1:ls_itetation_number
            
            
            % update the residual
            res = b - system_matrix * x;
            
            % getb the norm of the update direction
            res_norm(iter) = norm(res);
            
            % terminate the cg algorithm, if the norm of residual increases
            if iter > 1
                if res_norm(iter) > res_norm(iter-1)
                    break
                end
            end
            
            % update the direction
            update_direction = step_length.coeff_gridpoints' .* (system_matrix' *...
                (step_length.coeff_res .* res));
            
            % get the norm of the update direction
            res_norm(iter) = norm(update_direction);
            
            % update the solution
            x = x + step_length.step_length * update_direction;
            
            % display the number of iteration
            disp(['The number of inner iteration is:' num2str(iter)]);
        end
         
  
end

% get the residual norm
residual_norm = res_norm(1); 
  
% terminate the time counter
solve_time = toc(ts);

% calculate the first norm of the residual
disp(['The norm of the residual after solving the linear subproblem is:' num2str(residual_norm, '%1.3f')]);

end