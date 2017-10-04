function [ X, Y ] = computeEllipse( mx, my, a, b, angle )
% Returns the 360 X and Y coordinates of an ellipse
% Code inspired by wikipedia article on ellipses
% Inputs:
% mx     -- X coordinate of the midpoint
% my     -- Y coordinate of the midpoint
% a      -- Big semi-axis 
% b      -- Small semi-axis
% angle  -- Rotation of the ellipse

     alpha = linspace(0, 2*pi, 360);
     X = mx + (a * cos(alpha) * cos(angle) - b * sin(alpha) * sin(angle));
     Y = my + (a * cos(alpha) * sin(angle) + b * sin(alpha) * cos(angle));
    
end

