%% 
% Given initial lat/lon,
% the great circle angle "s", and the azimuth "psi" (measured clockwise
% from the North), return the latitude and longitude of the final point
% "phi", "lambda".
% Inputs:
%   long0: initial longitude, deg
%   lat0: initial latitude, deg
%   s: great circle angle, deg (aka range)
%   psi: azimuth, deg
% Outputs:
%   lambda: latitude of final point
%   phi:    longitude of final point
function [lambda, phi] = track_greatCirc(long0, lat0, s, az)
%% Store trig functions

sphi0 = sind(lat0); cphi0 = cosd(lat0);
ss = sind(s); cs = cosd(s);
spsi = sind(az); cpsi = cosd(az);

%% Latitude, longitude
phi = 90 - acosd(sphi0*cs + cphi0*ss*cpsi);
lambda = long0 + atan2d(ss*spsi, cphi0*cs - sphi0*ss*cpsi);

end