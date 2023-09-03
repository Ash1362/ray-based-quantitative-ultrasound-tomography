function [indices, tri, bc, nan_grid_binary] = gridDataFast2DModified(x, y, xi, yi)
% Copyright (c) 2012, Chao Huang
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% Note that the original version of the function has been modified according
% to the permission given by the owner of the copyright. Modifier: Ashkan Javaherian

% enforce x, y, xi, yi to be column vectors
x = x(:);
y = y(:);
xi = xi(:);
yi = yi(:);

% check if the MATLAB version is new enough to have DelaunayTri
if exist('DelaunayTri') %#ok<EXIST>
    % triangulize the data
    tri = DelaunayTri(x, y);  
else
    % call the old version of gridDataFast using tsearch
    [tri, bc] = gridDataFast2D_tsearch(x, y, xi, yi);
    return
end

% catch trinagulation error
if isempty(tri)
    error('Data cannot be triangulated.');
end

% find the nearest triangle and the corresponding Barycentric coordinates
[t, bc] = pointLocation(tri, [xi yi]);

% check points are valid
% if any(isnan(t))
  %  error('Cartesian points must lie within the k-space grid defined by kgrid');
% end

nan_grid_binary = ~isnan(t);

% get the indices of the points on which the interpolation is done
indices = repmat((1:nnz(nan_grid_binary))', 1, 3);

% get the barycebtric coordinates for the vertices of the triangle (voxel)
bc = bc(nan_grid_binary, :);

% get the vertices of the triangle (voxel)
tri = tri(t(nan_grid_binary), :);
end