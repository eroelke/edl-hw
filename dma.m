%% ASEN 6015 HW
% Drag Modulated Aerocapture Code
%  simplified, planar equations of motion
clear;clc;
%% custom parameters
% planet data
[H, rP, mu, rho0, g, hatm] = get_planet_data('Earth');

% entry state
vatm = 13;  %entry velocity, km/s
efpa = 5.39;  % entry flight path angle, deg (+ below horizon for 3dof eqs)

% vehicle
beta = 200; % choose aref from beta and mass
mass = 1000; %kg
cD = 1.05; % drag coefficient
cL = 0; % lift coefficient
Aref = mass / (cD * beta);

b21 = 4;    %beta ratio
cD_2 = cD;  %assume same cd
mass_2 = 200;
beta_2 = b21 * beta;
Aref_2 = mass_2 / (cD_2 * beta_2);

% guidance
g_rate = 0.5;   %guidance call rate, Hz
ha_tgt = 42164;    %GEO cause why not, km
tj0 = 100;  %s, initial guess

% sim
rate = 100;   %simulation integration rate, Hz
t_max = 1000;    %maximum time, s
t_init = 0; %initial time, s
hmin = 25; %min altitude for termination, km

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
p.hatm = hatm;
p.hmin = hmin;
% state
gamma0 = efpa * pi/180;
r0 = (p.hatm * 1000);  %planet relative position
v0 = vatm * 1000; %planet relative velocity
x0 = [r0;v0;gamma0];    %inertial state vector
% guidance
guid.tj = tj0;
guid.dt = 1/rate;
guid.tmax = 1000;   %max integration time
guid.h_tgt = ha_tgt * 1e3;
guid.dtj = 0;
guid.jflag = false;    %jettison flag
guid.m2 = mass_2;   % should be an array but we're just doing SEJ
guid.cd2 = cD_2;
guid.aref2 = Aref_2;

% sim
dt = 1/rate; %integration step size
tvec = t_init : dt : t_max;   %time vector
sg_ratio = round(rate/g_rate);
% preallocate stuff
x = nan(length(x0), length(tvec));
h = nan(length(x), 1); fpa = h; t = h; rho = h;
v = h; r = h;
% v = nan(1,length(tvec)); r = v;

%% run sim
x(:,1) = x0;
h(1) = p.hatm;
v(1) = vatm;
fpa(1) = -efpa;
t(1) = tvec(1);
control_flag = false;
for i = 2:length(x)
    % try to run guidance (if we haven't jettisoned)
    if (mod(i - 1, sg_ratio) == 0 && ~guid.jflag)
        guid = run_npc(x(:,i-1), tvec(i-1), veh, p, guid);
    end
    
    % integrate - 4th order runge kutta
    x(:,i) = run_rk4(x(:,i-1), veh, p, dt);
    
    % check jettison
    if (tvec(i) >= guid.tj && ~guid.jflag)
        guid.jflag = true;
        veh.m = guid.m2;
        veh.cd = guid.cd2;
        veh.aref = guid.aref2;
    end

    % isolate individual params for easier output
    h(i) = x(1,i) / 1000;
    v(i) = x(2,i) / 1000;
    fpa(i) = -x(3,i) * 180/pi;
    
    % check exit conditions
    if (h(i) <= p.hmin || h(i) >= p.hatm || v(i) <= 0)
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
rp = (min(h) * 1e3) + p.r; %periapsis is at min flight altitude
[r_a, dr_a] = calc_ra(rF, vF, rp, guid.h_tgt, p);
r_a = r_a / 1000;   %m -> km
dr_a = dr_a / 1000; %m -> km
h_a = r_a - (p.r / 1000);   %apoapsis altitude

fprintf(['Apoapsis: %4.2f km\n' ...
    'Target: %4.2f km\n' ...
    'Apoapsis Error: %4.2f km\n'],h_a, guid.h_tgt / 1000, dr_a);

plot(v,h)

% end




