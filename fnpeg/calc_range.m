%% calc_range.m
%   calculate range from spherical coords
% Inputs:
%   lat0: initial latitude, deg
%   lon0: initial longitude, deg
%   lat: current latitude, deg
%   lon: current longitude, deg
% 
function range = calc_range(lat0, lon0, lat, lon)

range = distance(lat0, lon0, lat, lon);
end