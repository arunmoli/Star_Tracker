function L = cartesian(alpha, delta) 
% Compute unit vector of star position 
L = [cos(delta) * cos(alpha) ; % usually on x axis */
     cos(delta) * sin(alpha) ; % usually on y axis */ 
     sin(delta) ];  % usually on z axis */ 
