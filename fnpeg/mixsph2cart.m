function xCart = mixsph2cart(xSph)
%%MIXSPH2CART Transform from "mixed spherical" coordinates to Cartesian
%%coordinates.
% 
% Transform the "mixed spherical" state (r, \theta, \phi, v, \gamma, \psi)
% to Cartesian position and velocity in the planet-fixed frame (rx, ry, rz,
% vx, vy, vz).
% No transport terms are included. The mixed spherical state is assumed to
% be already in the planet-fixed frame. See section 2.3.1 of Reference [1].
% 
% Angles are in radians.
% 
% The terminology "mixed spherical" stems from the fact that the position
% and velocity are both expressed in spherical coordinates, but these are
% referred to different reference frames.
%
% INPUTS:
%    * xSph, double Nx6 array, mixed spherical state
% 
% OUTPUTS:
%    * xCart, double Nx6 array, Cartesian state
%
% AUTHOR:
%    Davide Amato, CU Boulder, davide.amato@colorado.edu.
%
% REFERENCES:
%    [1] Amato, Davide. “Fall 2019 individual report.” CU Boulder, Boulder,
%    CO. September 2019.
% 
%% Unpack vector and save trig functions

r = xSph(:,1); theta = xSph(:,2); phi = xSph(:,3);
v = xSph(:,4); gamma = xSph(:,5); psi = xSph(:,6);

cth = cos(theta); sth = sin(theta);  
cph = cos(phi);   sph = sin(phi);

%% Cartesian position vector

rx = r .* cph .* cth;
ry = r .* cph .* sth;
rz = r .* sph;

%% Cartesian velocity (ENR frame)

vE = v .* cos(gamma) .* sin(psi);
vN = v .* cos(gamma) .* cos(psi);
vR = v .* sin(gamma);

vx = - vE .* sin(theta) - vN .* sin(phi) .* cos(theta) + vR .* cos(phi) .* cos(theta);
vy =   vE .* cos(theta) - vN .* sin(phi) .* sin(theta) + vR .* cos(phi) .* sin(theta);
vz =                      vN .* cos(phi)               + vR .* sin(phi)              ;

%% Pack output
xCart = [rx, ry, rz, vx, vy, vz];