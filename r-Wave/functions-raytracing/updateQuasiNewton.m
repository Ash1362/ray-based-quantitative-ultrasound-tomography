function [updated_polar_initial_direction, updated_polar_direction_residual,...
    updated_derivative_matrix, updated_norm_res, norm_search_direction] = updateQuasiNewton(solve_ray,...
    polar_direction_receiver, polar_initial_direction, polar_direction_residual,...
    derivative_matrix, bounds, zeta, kappa, tau, eta, sigma_max, sigma_min, para)
%UPDATEQUASINEWTON solves an iteration of ray linking inverse problem using the Newton's method
%
% DESCRIPTION:
% updateQuasiNewton solves an iteration of ray linking inverse problem using
% the Qausi-Newton method. This approach is used only for ray linking in 3D case.
% USAGE:
%
%
% INPUTS:
%      solve_ray                - the function handle for solving the ray and the
%                                 parameters associated with the end point of the ray
%      polar_direction_receiver - the polar direction of a geometrical
%                                 vector from emitter to receiver 
%      polar_initial_direction  - a 2 x 1 vector of the polar initial
%                                 direction of the last update of the ray 
%      polar_direction_residual - a 2 x 1 vector of the polar residual, the discrepancy
%                                 between the polar unit vector from emitter to 
%                                 the last point of the ray and a polar unit vector from emitter to the receiver
%      derivative_matrix        - a 2 x 2 derivative matrix, the derivative of
%                                 the polar direction of a geomterical vector
%                                 from emitter to the last point of the ray
%                                 with respect to the polar initial direction
%                                 of the ray
%      bounds                   - a 2 x 2 matrix of bounds for the update of the polar initial direction.
%                                 the first and second columns correspond to the lower and upper bounds, respectively.
%      zeta                    - a parameter controlling the distance from
%                                an associted bound after a projection to
%                                the interior of the feasible set
%      kappa                   - the maximum permissible reduction for the amplitude of search direction
%                                by a projection to the feasible set (safe-guard).
%      tau                     - a scalar controlling the update of the derivative matrix
%      eta                     - a parameter for controlling chages of tau
%      sigma_max               - a scalar bound for the maximum to minimim singular value 
%                                of the derivative matrix
%      sigma_min               - a scalar bound for the minimim singular value of the derivative matrix
%      para                    - the optional inputs for ray linking using
%                                Quasi-Newton approach
% OPTIONAL INPUTS:
%
% OUTPUTS:
%      updated_polar_initial_direction  - a 2 x 1 vector of the update of polar initial
%                                 direction of the last update of the ray
%      updated_polar_direction_endpoint - a 2 x 1 vector of the update of polar 
%                                 direction of a geometrical a unit vector from 
%                                 emitter to the last point of the updated ray  
%      updated_polar_direction_residual - a 2 x 1 vector of the update of polar residual, 
%                                 the discrepancy between the polar unit vector from emitter to 
%                                 the last point of the ray and a polar unit vector from emitter
%                                 to the receiver
%      updated_derivative_matrix - the 2 x 2 updated derivative matrix, the derivative of
%                                 the polar direction of a geomterical vector
%                                 from emitter to the last point of the ray
%                                 with respect to the polar initial direction
%                                 of the ray
%      updated_norm_res                 - the updated norm of residual 
%      norm_search_direction    - the size of the search direction. this is
%                                 used for termination of iterations, 
%                                 if it is very small, because it may
%                                 cause instability for the next iterates


% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 



% update the search direction
search_direction = - derivative_matrix \ polar_direction_residual;

% modify the search direction to project the updated initial direction
% to the interior of feasible set
low_indices = polar_initial_direction + search_direction < bounds(:,1);
up_indices = polar_initial_direction + search_direction > bounds(:,2);

if any(low_indices)
    % define a vector psi
    psi = zeta * (bounds(low_indices, 1) - polar_initial_direction(low_indices))...
        ./ search_direction(low_indices);
    % modify the search direction for a projection to the interior of the feasible set
    search_direction(low_indices) = sign(psi) .* max(abs(psi), kappa).*...
        search_direction(low_indices);
end


if any(up_indices)
    % define a vector psi
    psi = zeta*(bounds(up_indices, 2) - polar_initial_direction(up_indices))...
        ./search_direction(up_indices);
    % modify the search direction for a projection to the interior of the feasible set
    search_direction(up_indices) = sign(psi).* max(abs(psi), kappa) .*...
        search_direction(up_indices);
end


% calcualte the size of the search direction
norm_search_direction = norm(search_direction);

% update the polar initial direction of the ray using the
% bounded search direction
updated_polar_initial_direction = polar_initial_direction +  search_direction;

% calculate the ray using the updated polar initial direction,
% and update the end point of the ray
[updated_polar_direction_residual, ~, ~, ~]...
    = feval(solve_ray, updated_polar_initial_direction, polar_direction_receiver, false);


% update the norm of the residual
updated_norm_res = norm(updated_polar_direction_residual);


    
if strcmp(para.Method,'BFGS')
    
      % to include the error in inversion of the derivative matrix,
      % the previous update of the residual is  replaced by - derivative_matrix * search_direction 
       y = updated_polar_direction_residual + derivative_matrix * search_direction; 
end        
      


if ~para.smooth
    
      switch para.Method
            case 'Good-Broyden'
                updated_derivative_matrix = derivative_matrix + (updated_polar_direction_residual * search_direction')...
                    ./(search_direction' * search_direction);
            case 'BFGS'
                if y' * search_direction > 0
                updated_derivative_matrix = derivative_matrix ...
                    - (derivative_matrix *(search_direction * search_direction')* derivative_matrix)...
                    ./(search_direction' * derivative_matrix * search_direction)...
                    + (y * y')./(y' * search_direction);
                end
      end
else
    
     % choose a lower bound for the eigenvalues of the derivative matrix
        lambda = min(updated_norm_res, sigma_min);
        
     % start the update factor 
        if strcmp(para.Method, 'Good-Broyden')
               % increment factor for changes of eta
                eta_c = eta;
                eta = 0;
                % the sigm of eta
                eta_s = -1;
        end   
      
  for i=1:20
            switch para.Method
                case 'Good-Broyden'      
                    updated_derivative_matrix = derivative_matrix + ...
                        (tau + eta_s*eta) * (updated_polar_direction_residual * search_direction')...
                    ./(search_direction' * search_direction); 
                
                % change the sign 
                eta_s = - eta_s;
                % increase eta for positive eta_s
                if eta_s >0
                    eta = eta + eta_c;
                end
                
                
                case 'BFGS'
                     updated_derivative_matrix = derivative_matrix ...
                   + tau * ( -(derivative_matrix *(search_direction * search_direction')* derivative_matrix)...
                    ./(search_direction' * derivative_matrix * search_direction)...
                    + (y * y')./(y' * search_direction));  
                
                % reduce tau by a reduction factor eta
                tau = eta*tau ;
            end
            
            if all(isfinite(updated_derivative_matrix(:)))
            if  cond(updated_derivative_matrix) < sigma_max &&  all(svd(updated_derivative_matrix)>lambda )
                break;
            end
            else
                % reject the updated derivative matrix, if it has at least
                % one infinite element
                updated_derivative_matrix = derivative_matrix;
                break;
            end
                
  end
    

 % reject the BFGS-based updated derivative matrix, if it violates the
 % conditions for the singular values of the deriavtive matrix
 if strcmp(para.Method,'BFGS')  && ( (cond(updated_derivative_matrix)>sigma_max...
         ||  any(svd(updated_derivative_matrix(:))<lambda)) )
    updated_derivative_matrix = derivative_matrix;
 end
  
                 
end

% reject the updated derivative matrix, if it has at least one infinite elements 
if  any(~isfinite(derivative_matrix(:)))
            updated_derivative_matrix = derivative_matrix;
end










end