%% run_rk4.m
%   4th runge kutta integrator for edl
% 
function x_new = run_rk4(x0, veh, p, dt)

k1 = eom_3dof(x0, veh, p);
k2 = eom_3dof(x0 + (dt * k1/2), veh, p);
k3 = eom_3dof(x0 + (dt * k2/2), veh, p);
k4 = eom_3dof(x0 + k3 * dt, veh, p);

x_new = x0 + ... 
    (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end
