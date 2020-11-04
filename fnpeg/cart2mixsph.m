%% cart2mixsph.m
%   Transform from "mixed spherical" coordinates to Cartesian coordinates.
% 
% Transform the Cartesian state (rx, ry, rz, vx, vy, vz) to 
% (r, theta, phi, v, gamma, psi).
% 
% The Cartesian state is assumed to be in the planet-fixed frame.
% 
function xSph = cart2mixsph(xCart)
%extract vars
rx = xCart(:,1);
ry = xCart(:,2);
rz = xCart(:,3);
vx = xCart(:,4);
vy = xCart(:,5);
vz = xCart(:,6);

% init
r = zeros(size(xCart,1),1);
v = zeros(size(xCart,1),1);
for i = 1:size(xCart,1)
    r(i) = norm(xCart(i,1:3));
    v(i) = norm(xCart(i,4:6));
end

theta = atan2(ry, rx);  %lon
phi = asin(rz ./ r);    %lat

%v
vE = -vx .* sin(theta)             + vy .* cos(theta);
vN = -vx .* sin(phi) .* cos(theta) - vy .* sin(phi) .* sin(theta) + vz .* cos(phi);
vR =  vx .* cos(phi) .* cos(theta) + vy .* cos(phi) .* sin(theta) + vz .* sin(phi);

psi   = atan2(vE, vN);  %heading
gamma = asin(vR ./ v);  %fpa

xSph = [r, theta, phi, v, gamma, psi];
end

