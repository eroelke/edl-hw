%% eom_3dof.m
%  edl 3dof equations of motion
% 
% Inputs:
%   x: state vector
%   veh: vehicle struct
%   atm: atmospheric data
% 
function [xdot, rho] = eom_3dof(x, veh, p)

% extract params
h = x(1); %altitude
v = x(2); %inertical v magnitude
gamma = x(3);   %fpa, rad
% gamma = -gamma; %fpa, rad, (+) below horizon

m = veh.m;  %current mass, kg
cd = veh.cd;    %drag coeff
cl = veh.cl;    %lift coeff
aref = veh.aref;    %ref area, m2
beta = m / (aref * cd); %ballistic coefficient, kg/m2

rho = p.rho0 * exp(-h / p.H);

% eqs
rv2b = (rho * v^2) / (2 * beta);    %common eom param

rdot = -v * sin(gamma);
vdot = -rv2b + p.g*sin(gamma);
gammadot = (1/v) * (-(v^2 * cos(gamma))/(h + p.r) - ((cl/cd) * rv2b) + p.g*cos(gamma));

% concat dxdt array
xdot = [rdot;vdot;gammadot];
end