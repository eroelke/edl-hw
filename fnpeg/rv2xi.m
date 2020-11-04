%% rv2xi.m
%   calc negative of energy for fnpeg guidance
% 
function energy = rv2xi(mu,r,v)
energy = 0.5 * v.^2  -  mu ./ r;
end