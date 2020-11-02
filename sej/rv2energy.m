function energy = rv2energy(mu,r,v)
%%%RV2ENERGY Orbital energy from position and velocity.
%
% Given the position magnitude r, velocity magnitude v, and gravitational
% parameter mu, compute the orbital energy. The function accepts either 
% dimensional or normalized inputs and returns outputs accordingly; however,
% units have to be consistent between mu and r, v.
% 
% Author:
%    Davide Amato, CU Boulder, davide.amato@colorado.edu
%

energy = 0.5 * v.^2  -  mu ./ r;

end