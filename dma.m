%% ASEN 6015 HW
% Drag Modulated Aerocapture Code
%  simplified, planar equations of motion
clear;clc;
%% custom parameters
% planet data
[H, rP, mu, rho0, g, h0] = get_planet_data('Earth');

% entry state
vatm = 13;  %entry velocity, km/s
efpa = 5.35;  % entry flight path angle, deg (+ below horizon for 3dof eqs)

% vehicle
beta = 200; % choose aref from beta and mass
mass = 1000; %kg
cD = 1.05; % drag coefficient
cL = 0; % lift coefficient
Aref = mass / (cD * beta);

% guidance
g_rate = 0.1;   %guidance call rate, Hz
ha_tgt = 42164;    %GEO cause why not, km

% sim
rate = 10;   %simulation integration rate, Hz
t_max = 1000;    %maximum time, s
t_init = 0; %initial time, s
h_min = 25; %min altitude (for quicker testing)

%% set up sim
% vehicle struct
veh.m = mass;
veh.cd = cD;
veh.cl = cL;
veh.aref = Aref;
% atmosphere
p.H = H;
p.rho0 = rho0;
p.r = rP;
p.mu = mu;
p.g = g;
% state
gamma0 = efpa * pi/180;
r0 = (h0 * 1000);  %planet relative position
v0 = vatm * 1000; %planet relative velocity
x0 = [r0;v0;gamma0];    %inertial state vector

% sim
dt = 1/rate; %integration step size
tvec = t_init : dt : t_max;   %time vector
sg_rate = round(rate/g_rate);
% preallocate stuff
x = nan(length(x0), length(tvec));
h = nan(length(x), 1); fpa = h; t = h; rho = h; v = h; r = h;
% v = nan(1,length(tvec)); r = v;

%% run sim
x(:,1) = x0;
h(1) = h0;
v(1) = vatm;
fpa(1) = -efpa;
t(1) = tvec(1);
control_flag = false;
for i = 2:length(x)
    % integrate - 4th order runge kutta
       
    [k1, rho(i-1)] = eom_3dof(x(:,i-1), veh, p);
    k2 = eom_3dof(x(:,i-1) + (dt * k1/2), veh, p);
    k3 = eom_3dof(x(:,i-1) + (dt * k2/2), veh, p);
    k4 = eom_3dof(x(:,i-1) + k3 * dt, veh, p);

    x(:,i) = x(:,i-1) + ... 
        (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    % check appreciable atmosphere
    if (abs(k1(2)) > 0.5 && ~control_flag)
        control_flag = true;
    end

    % try to run guidance
%     if (control_flag && mod(i - 1, sg_rate) == 0)
%         %todo
%     end

    % isolate individual params for easier output
    h(i) = x(1,i) / 1000;
    v(i) = x(2,i) / 1000;
    fpa(i) = -x(3,i) * 180/pi;
    
    % check exit conditions
    if (h(i) <= h_min || ... 
    (control_flag && h(i) > h0) || ...
    v(i) <= 0)
        break;
    end
end

%% post process
if (~isnan(h(end)))
    idxend = length(tvec);
else
    idxend = find(isnan(h) == true, 1) - 1;
end
hF = h(idxend) * 1000;
rF = hF + p.r;
vF = v(idxend) * 1000;
fpaF = fpa(idxend) * pi/180;

%calc error
rp = (min(h) * 1000) + p.r; %periapsis is at min flight altitude
[r_a, dr_a] = calc_apoapsis(rF, vF, rp, ha_tgt, p);
h_a = r_a - (p.r / 1000);

fprintf(['Apoapsis: %4.2f km\n' ...
    'Target: %4.2f km\n' ...
    'Apoapsis Error: %4.2f km\n'],r_a, ha_tgt, dr_a);

plot(v,h)



% get sma and eccentricity since that s
function [oe, r_a, dr_a] = calc_ra(rv, mu, tgt)
r = rv(1:3);
v = rv(4:6);
    
rmag = norm(r);     %magintude of position vector (m)
vmag = norm(v);     %magnitue of velocity vector (m/s)

a = rmag/(2 - (rmag*vmag^2)/mu);        %semimajor axis (vis viva equation)
h = cross(r,v);         %angular momentum
% hmag = norm(h);
e = norm(cross(v,h)/mu - r/rmag); %eccentricity
oe = [a e];

r_a = (a * (1 + e));

end

function [H, rP, mu, rho0, g, h0] = get_planet_data(planet)

switch (planet)
    case 'Venus'
        H = 15900;
        rP = 6.0518e6;
        mu = 3.24858592079e14;
        g = 8.870;
        rho0 = 65;
        h0 = 150;
    case 'Earth'
        H = 8500;  %atm scale height, m
        rP = 6378e3;  %radius, m
        mu = 3.9860e14; %gravitational param, m^3/s^2
        rho0 = 1.225;  %surface density, kg/m3
        g = 9.81;   %grav accel, m/s2
        h0 = 125;
    case 'Mars'
        rP = 3.3962e6;
        mu = 4.2830e13;
        H = 11100;
        rho0 = 0.02;
        g = 3.710;
        h0 = 100;
    otherwise
        error('Bad Planet');
end
end

function [r_a, dr_a] = calc_apoapsis(rmag, vmag, rp, ha_tgt, p)
energy = (vmag^2/2 - p.mu/rmag);
a = -(p.mu/(2 * energy));
% rp = ;
r_a = 2*a - rp;
e = (r_a - rp)/(r_a + rp);
r_tgt = (ha_tgt * 1000) + p.r;  % apoapsis target radius, m
dr_a = r_a - r_tgt;
r_a = r_a / 1000;
dr_a = dr_a / 1000;
end




