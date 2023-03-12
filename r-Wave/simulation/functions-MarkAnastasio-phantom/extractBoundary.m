function boundary_mask = extractBoundary(mask, BC)
%EXTRACTSURFOFLOGVOL Extract boundary mask of logical matrix or volume
%
% DESCRIPTION:
%       extractBoundary extracts the boundaray of a logical volume or matrix
%       by assuming that the background is false and any true voxels whos
%       4/6-neighbourhood is not all true are boundary
%
% USAGE:
%       boundary_mask = extractBoundary(mask)
%       boundary_mask = extractBoundary(mask,BC)
%
% INPUTS:
%       mask - a 2D/3D logical array
%
% OPTIONAL INPUTS:
%       BC - boundary conditions, 'natural', 'Neumann', 'NB' or '0': determines
%            how to deal with voxels on the boundary whose neighbours are
%            all "true". For 'natural', 'Neumann', 'NB' they are not included
%            in the boundary mask, for "0", they are.
%
% OUTPUTS:
%       boundary_mask - a 2D/3D logical volume
%
% ABOUT:
%       author          - Felix Lucka
%       date            - 15th March 2017
%       last update     - 14th Nov   2017
%
%
% See also

% check user defined value for BC, otherwise assign default value
if(nargin < 2)
    BC = 'natural';
end

% get size of volume and create output mask
dim           = ndims(mask);
[nx, ny, nz]  = size(mask);
boundary_mask = false(size(mask));

% loop over inner voxels
for i=2:nx-1
    for j=2:ny-1
        
        switch dim
            
            case 2
                if(mask(i, j))
                    
                    % its an element of the mask, check whether is lies on the
                    % boundary by checking if one of its 4 neighbours is not mask
                    if( ~mask(i-1, j) || ...
                        ~mask(i+1, j) || ...
                        ~mask(i,   j-1) || ...
                        ~mask(i,   j+1))
                        
                        % the voxel is part of the boundary
                        boundary_mask(i, j) = true;
                        
                    end
                    
                end
                
            case 3
                for k=2:nz-1
                    if(mask(i, j, k))
                        
                        % its an element of the mask, check whether is lies on the
                        % boundary by checking if one of its 6 neighbours is not mask
                        if( ~mask(i-1, j, k)     || ...
                            ~mask(i+1, j, k)     || ...
                            ~mask(i,   j-1, k)   || ...
                            ~mask(i,   j+1, k)   || ...
                            ~mask(i,   j,   k-1) || ...
                            ~mask(i,   j,   k+1))
                            
                            % the voxel is part of the boundary
                            boundary_mask(i, j, k) = true;
                            
                        end
                        
                    end
                end
                
            otherwise
                
                error('dimension of mask must be 2 or 3.')
                
        end
        
    end
end

% treat the boundaries
switch BC
    
    case {'Neumann', 'NB', 'natural'}
        
    case '0'
        
        switch dim
            case 2
                boundary_mask([1, end], :) = mask([1, end], :);
                boundary_mask(:, [1, end]) = mask(:, [1, end]);
            case 3
                boundary_mask([1, end], :, :) = mask([1, end], :, :);
                boundary_mask(:, [1, end], :) = mask(:, [1, end], :);
                boundary_mask(:, :, [1, end]) = mask(:, :, [1, end]);
                
        end
        
    otherwise
        
        error('unknown boundary condition.')
        
end

end
