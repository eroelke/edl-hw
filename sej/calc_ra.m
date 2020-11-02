%% calc_apoapsis.m
%   janky method of calculating apoapsis and apoapsis error
%   for dma with scalars instead of vectors
% 
% Assumptions:
%   1. rp is calculated from the minimum altitude in the aerocapture
%   trajectory 
% 
% Inputs:
%   rmag:   vehicle position scalar, m
%   vmag:   vehicle velocity scalar, m/s
%   rp:     periapsis radius, m 
%   ha_tgt: apoapsis target altitude, m
%   p:      planet struct
% 
function [ra, dr_a] = calc_ra(rmag, vmag, rp, ha_tgt, p)
energy = (vmag^2/2 - p.mu/rmag);
a = -(p.mu/(2 * energy));
ra = 2*a - rp;
e = (ra - rp)/(ra + rp);
r_tgt = ha_tgt + p.r;  % apoapsis target radius, m
dr_a = ra - r_tgt;
end