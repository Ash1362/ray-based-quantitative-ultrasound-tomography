function [gradient_x, gradient_y, gradient_z, diffusion_coefficient] = ...
    calcDirectionalGradient(field_direction, field_background,...
    smoothing_parameter, grid_spacing)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


switch ndims(field_background)
    case 2
        
        
        % compute the diffusion coefficient using the background field
        gradient_x = padarray(diff(field_background, 1, 1), [1, 0], 0, 'post')/ grid_spacing;
        gradient_y = padarray(diff(field_background, 1, 2), [0, 1], 0, 'post')/ grid_spacing;
        
        diffusion_coefficient = 1./ (1 + (gradient_x.^2 +...
            + gradient_y.^2)./(smoothing_parameter^2)  );
        
        % compute and overwite the directional gradients using the update direction field
        gradient_x = padarray(diff(field_direction, 1, 1), [1, 0], 0, 'post')/ grid_spacing;
        gradient_y = padarray(diff(field_direction, 1, 2), [0, 1], 0, 'post')/ grid_spacing;
        
        
        gradient_z = 0;
        
    case 3
        
        % compute the diffusion coefficient using the background field
        gradient_x = padarray(diff(field_background, 1, 1), [1, 0, 0] , 0, 'post')/ grid_spacing;
        gradient_y = padarray(diff(field_background, 1, 2), [0, 1, 0] , 0, 'post')/ grid_spacing;
        gradient_z = padarray(diff(field_background, 1, 3), [0, 0, 1] , 0, 'post')/ grid_spacing;
        
        
        diffusion_coefficient = 1./ ( 1 + (gradient_x.^2 +...
            + gradient_y.^2  + gradient_z.^2)./(smoothing_parameter^2) );
        
        % compute and overwite the directional gradients using the update direction field
        gradient_x = padarray(diff(field_direction, 1, 1), [1, 0, 0] , 0, 'post')/ grid_spacing;
        gradient_y = padarray(diff(field_direction, 1, 2), [0, 1, 0] , 0, 'post')/ grid_spacing;
        gradient_z = padarray(diff(field_direction, 1, 3), [0, 0, 1] , 0, 'post')/ grid_spacing;
        
        
end

end