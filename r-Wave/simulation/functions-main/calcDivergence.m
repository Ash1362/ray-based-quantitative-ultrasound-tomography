function [spatial_variation] = calcDivergence(field_x, field_y, field_z,...
     mask, grid_spacing)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


switch ndims(field_x)
    
    case 2
               
%spatial_variation_x = cat(1, field_x(1, :), diff(field_x(1:end-1, :), 1, 1), - field_x(end, :));
%spatial_variation_y = cat(2, field_y(:, 1), diff(field_y(:, 1:end-1), 1, 2), - field_y(:, end));

spatial_variation_x = cat(1, field_x(1, :), diff(field_x, 1, 1) );
spatial_variation_y = cat(2, field_y(:, 1), diff(field_y, 1, 2) );

spatial_variation = - (spatial_variation_x(mask) + spatial_variation_y(mask))/grid_spacing;


    case 3
        
%spatial_variation_x = cat(1, field_x(1, :, :), diff(field_x(1:end-1, :, :), 1, 1), - field_x(end, :, :));
%spatial_variation_y = cat(2, field_y(:, 1, :), diff(field_y(:, 1:end-1, :), 1, 2), - field_y(:, end, :));
%spatial_variation_z = cat(3, field_z(:, :, 1), diff(field_z(:, :, 1:end-1), 1, 3), - field_z(:, :, end));

spatial_variation_x = cat(1, field_x(1, :, :), diff(field_x, 1, 1));
spatial_variation_y = cat(2, field_y(:, 1, :), diff(field_y, 1, 2));
spatial_variation_z = cat(3, field_z(:, :, 1), diff(field_z, 1, 3));

spatial_variation = - (spatial_variation_x(mask) + spatial_variation_y(mask) + spatial_variation_z(mask))/grid_spacing;

end



end