function xSph = cart2mixsph(xCart)
%%CART2MIXSPH Transform from "mixed spherical" coordinates to Cartesian
%%coordinates.
% 
% Transform the Cartesian state (rx, ry, rz, vx, vy, vz) to "mixed
% spherical" coordinates (r, \theta, \phi, v, \gamma, \psi). No transport
% terms are included. The Cartesian state is assumed to be in the
% planet-fixed frame. See section 4.2 of Reference [1].
%
% The terminology "mixed spherical" stems from the fact that the position
% and velocity are both expressed in spherical coordinates, but these are
% referred to different reference frames.
% 
% All angles are in radians.
%
% INPUTS:
%    * xCart, double Nx6 array, mixed spherical state
% 
% OUTPUTS:
%    * xSph, double Nx6 array, Cartesian state
%
% AUTHOR:
%    Davide Amato, CU Boulder, davide.amato@colorado.edu.
%
% REFERENCES:
%    [1] Amato, Davide. “Fall 2019 individual report.” CU Boulder, Boulder,
%    CO. September 2019.
% 
%% Unpack vector

rx = xCart(:,1); ry = xCart(:,2); rz = xCart(:,3);
vx = xCart(:,4); vy = xCart(:,5); vz = xCart(:,6);

%% Radius, velocity, longitude, latitude

r = zeros(size(xCart,1),1); v = zeros(size(xCart,1),1); 
for i=1:size(xCart,1)
    r(i) = norm(xCart(i,1:3));
    v(i) = norm(xCart(i,4:6));
end
theta = atan2(ry, rx);
phi = asin(rz ./ r);

%% Velocity, heading, flight-path angle

vE = -vx .* sin(theta)             + vy .* cos(theta)                             ;
vN = -vx .* sin(phi) .* cos(theta) - vy .* sin(phi) .* sin(theta) + vz .* cos(phi);
vR =  vx .* cos(phi) .* cos(theta) + vy .* cos(phi) .* sin(theta) + vz .* sin(phi);

psi   = atan2(vE, vN);
gamma = asin(vR ./ v);

%% Pack outputs

xSph = [r, theta, phi, v, gamma, psi];

end

