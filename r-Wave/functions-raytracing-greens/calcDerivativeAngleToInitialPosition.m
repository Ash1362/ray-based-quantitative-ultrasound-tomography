function [directions] = calcDerivativeAngleToInitialPosition(directions,...
    dim, grid_position, transducer_position, perturbation_method)
%CALCDERIVATIVEANGLETOINITIALPOSITION computes the derivative of the
%direction of the ray to the initial position
%
% DESCRIPTION:
%       calcDerivativeAngleToInitialPosition computes the derivative of
%       the direction of the ray in the polar coordinates, i.e., perturbation to
%       angle of the ray, with respect to the initial position of the ray.
%       For each forward (resp. adjoint) ray, this derivative is computed using
%       finite differences and using the the rays initialised from the adjacent
%       emitters (resp. receivers).

%
% USAGE:
%
%
% INPUTS:
%       directions        - a 1 x num_emitter (resp. 1 x num_receiver) cell
%                           each containing the the Cartesian direction of the
%                           forward (resp. adjoint) rays initialised from an emitter
%                           (resp. receiver)
%       dim                - the number of dimensions
%       grid position       - if given, the num_points x dim Cartesian position of the grid
%                            points
%       transducer_position - if given, a num_transducer x dim matrix containing the
%                           Cartesian position of transducers along the
%                           cartesian coordintes
%       perturbation_method - get the method for computing perturbation to
%                           angle of the ray, with respect to the initial
%                           position of the ray. This can be set 'adjacent',
%                           or 'chain_rule'.

%
% Optional INPUTS:
%
%
% OUTPUTS:
%       directions     - a 1 x num_emitter (resp. 1 x num_receiver) cell
%                           each containing the the Cartesian direction of the
%                           forward (resp. adjoint) rays initialised from an emitter
%                           (resp. receiver), and the derivative of the
%                           angle of the rays with respect to the initial
%                           position of the ray
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 23.08.2020
%       last update     - 14.01.2023
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian




if nargin > 2

    % get the number of transducers
    num_transducer = size(transducer_position, 1);

    % get the number of dimensions
    if size(directions{1}, 2)-dim > 0

        directions_bent = cell(1, num_transducer);

        for ind_transducer = 1: num_transducer

            % create a struct for rays' directions perturbations and remove them
            % from the struct for rays' directions
            directions_bent{ind_transducer} = directions{ind_transducer}(:,end);
            directions{ind_transducer}(:,end)= [];
        end

    end

    if length(directions) ~= num_transducer
        error('The length of the cell array for ray directions is not correct.')
    end


    for ind_transducer = 1: num_transducer

        % add a third and fourth colums for straight rays' directions
        directions{ind_transducer} =  [directions{ind_transducer}, ...
            grid_position - transducer_position(ind_transducer, :)];
    end

else

    % get the perturbation method
    perturbation_method = 'adjacent';

    % get the number of transducers
    num_transducer = size(directions, 1);

    if size(directions{1}, 2)> dim

    for ind_transducer = 1: length(directions)
    directions{ind_transducer}(:,end) = [];
    end

    end


end


switch perturbation_method

    case 'adjacent'

        % make the rays' direction for all initial positions a matrix array
        directions = cat(3, directions{:});

        % add the last (resp. the first) initial positions as the first (resp.
        % last) index of the last dimension of array.
        % This is required for the centred finite difference along the last dimension
        % of the array for rays' directions
        directions = cat(3, directions(:,:, end), directions, directions(:, :, 1));

        % compute the signed angles between the two dimensional vectors included in the
        % second dimension of arrays by applying finite differences

        if nargin > 2
            angles = calcDirectionalAngle(permute(directions(:, dim+1:2*dim, 3:end), [2,1,3]), ...
                permute(directions(:, dim+1:2*dim, 1:end-2), [2,1,3]));
        else
            angles = calcDirectionalAngle(permute(directions(:,:, 3:end), [2,1,3]), ...
                permute(directions(:,:, 1:end-2), [2,1,3]));
        end
        % compute the absolute values of the computed angle derivatives with
        % respect to initial poistions, and then rechange size of the array
        angles = permute(1/2 * angles, [2,1,3]);

        % add the computed angles with respect to changes to the initial positions
        % to the original array of rays' directions
        directions = cat(2, directions(:, 1:dim, 2:end-1), angles);

        % make the resulting array matrix a cell array with the original length,
        % i.e., with length of the number of initial positions
        directions = permute(squeeze(mat2cell(directions, size(directions, 1), size(directions, 2),...
            ones(1, size(directions, 3)))), [2, 1]);

    case 'chain_rule'

        %error('This option was for testing and was deprectaed.')

        % get the perturbation factor
        perturbation_factor = 2*pi/num_transducer;

        % get the number of grid points
        num_point = size(grid_position, 1);

        for ind_transducer = 1: num_transducer


            % show the index of the transducer
            disp(['Computing the perurbation to rays angles for transducer:'...
                num2str(ind_transducer)])

            % get the initial vector perturbation
            initial_vector_perturb = - perturbation_factor *...
                [-transducer_position(ind_transducer, 2),...
                transducer_position(ind_transducer, 1)].';





            switch dim
                case 2

                    % get the direction for the current transducer and grid point
                    % directions_current = directions{ind_transducer}(:, dim + 1 : 2 * dim);

                    % get the angle perturbation on the grid point
                    %angles = 1./(vecnorm(directions_current).^2) .*...
                     %   ([-directions_current(:, 2), directions_current(:, 1)] *...
                     %   initial_vector_perturb);

                    % get the direction for the current transducer and grid point
                    directions_current = directions{ind_transducer}(:, 1:2);
                    directions_current = 1./vecnorm(directions_current, 2, 2) .* directions_current;
                    directions_current2 = directions{ind_transducer}(:, dim + 1 : 2 * dim);
                    distances = vecnorm(directions_current2, 2,2);
                    
                    directions_perturbation = 1./distances .* (repmat(initial_vector_perturb.', [size(directions_current, 1), 1])...
                        - ((directions_current * initial_vector_perturb) .* directions_current));

                    directions_current2 = 1./distances .* directions_current2;


                    % get the angle perturbation on the grid point
                    angles = (sum([-directions_current2(:, 2), directions_current2(:, 1)] .* directions_perturbation, 2)).';
 

                case 3

                    % allocate a zero vector for the angles on the grid points
                    angles = zeros(1, num_point);

                    for ind_point = 1: num_point

                        % get the direction for the current transducer and grid point
                        directions_current = directions{ind_transducer}(...
                            ind_point, dim + 1 : 2 * dim);

                        % get the angle perturbation on the grid point
                        angles(ind_point) = 1/(norm(directions_current)^2) *...
                            [-directions_current(2), directions_current(1)] *...
                            initial_vector_perturb;
                    end

            end


            % add the angles as a third column to rays' directions for the curret
            % tranducer
            directions{ind_transducer} = [directions{ind_transducer}(:, 1:dim),...
                angles.'];

        end

end


if exist('directions_bent')
    include_bent = true;
    disp(['The refraction effects are included in computing ' ...
        'the derivative of ray angles with respect to' ...
        'the initial positions'])
else
    include_bent = false;
end


if include_bent

    for ind_transducer = 1: num_transducer

        directions{ind_transducer}(:, end) = abs(directions{ind_transducer}(:, end) +...
            directions_bent{ind_transducer}(:,end));

    end

else

    for ind_transducer = 1: num_transducer

        directions{ind_transducer}(:, end) = abs(directions{ind_transducer}(:, end));

    end

end





end