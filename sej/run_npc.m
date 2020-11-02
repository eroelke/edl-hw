%% run_npc.m
%   numerical predictor corrector code
% 
% Inputs:
%   x0:     current vehicle state
%   t0:     current time
%   veh:    vehicle struct
%   p:      planet struct
%   g:      guidance struct
% 
function g = run_npc(x0, t0, veh, p, g)

e_tj = g.dt;
[ra, dra] = run_predictor(x0, t0, veh, p, g, 0);
[ra2, dra2]  = run_predictor(x0, t0, veh, p, g, e_tj);

%update term
drdt = (dra2 - dra)/e_tj;
dtj = dra/drdt; %f/f'

%bound
if dtj >= 20
    dtj = 20;
elseif dtj <= -20
    dtj = -20;
end

g.dtj = dtj;
g.ra = ra / 1000;   %km
g.dra = dra / 1000; %km
g.tj = g.tj - dtj;

end

% predictor for dma
function [ra, dra] = run_predictor(x0, t0, veh, p, g, e_tj)

flag = 0;
tj = g.tj + e_tj;   %for newton method
x = x0;
ind = 1;
t = t0;
jflag = false;
while ~flag
    x = run_rk4(x, veh, p, g.dt);
    h(ind) = x(1) / 1000;    %km
    v(ind) = x(2) / 1000;    %km/s
    
    % check jettison
    if (t >= tj && ~jflag)
        jflag = true;
        veh.m = g.m2;
        veh.cd = g.cd2;
        veh.aref = g.aref2;
    end
        

    % check exit conditions
    if (h(ind) >= p.hatm)
        flag = 1;
        break;
    elseif (h(ind) <= p.hmin || v(ind) <= 0)
        flag = 2;
        break;
    elseif (t >= g.tmax)
        flag = 3;
        break;
    end
    
    ind = ind + 1;
    t = t + g.dt;
end

rmag = (h(end) * 1e3) + p.r;
vmag = v(end) * 1e3;
rp = (min(h) * 1e3) + p.r;  %periapsis radius, m
[ra, dra] = calc_ra(rmag, vmag, rp, g.h_tgt, p);
end
