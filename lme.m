%% ASEN 6015 HW
% Q2: Lift Modulation FNPEG guidance code
% 
clear; clc; clearvars;
addpath('util/');
addpath('fnpeg/');

% Planetary constants
planet = get_planet_data('Mars');

% Reference quantities
lRef = planet.r;                             % (m)
tRef = sqrt(lRef^3/planet.mu);               % (s)

% settings
s_tgt = 7.952;      % deg, range target
h_tgt = 11e3;       % target altitude, m
v_tgt = 1.2e3;      % target velocity, m/s

%% Initial conditions at EI and vehicle parameters

% Initial conditions in spherical coordinates
ICs_sph.ti   = 0;                % initial time, s
ICs_sph.h0   = 130E3;            % initial altitude, m
ICs_sph.lon0 = 90;               % initial longitude, deg
ICs_sph.lat0 = 45;               % initial latitude, deg
ICs_sph.v0   = 4000;             % initial velocity, m/s
ICs_sph.fpa0 = -15;              % initial fpa, deg
ICs_sph.psi0 = 70;               % initial heading, deg

% Vehicle parameters (corresponding to HIAD)
vehicle.m0    = 49000;                                  % kg (initial mass)
vehicle.mDry  = 20000;                                  % kg (dry mass)
vehicle.B0    = 155;                                    % kg/m^2 (ballistic parameter)
vehicle.A     = 200;                                    % m^2 (vehicle reference area - ASSUMED)
vehicle.LD    = 0.15;                                   % Lift-to-drag ratio
vehicle.CD    = vehicle.m0 / (vehicle.B0 * vehicle.A);  % Drag coefficient
vehicle.CL    = vehicle.CD * vehicle.LD;                % Lift coefficient
% vehicle.TMax  = 800E3;                                  % Maximum thrust magnitude (N)
vehicle.TMin  = 200E3;                                  % Minimum thrust magnitude (N)

%% Guidance parameters
% GENERAL GUIDANCE PARAMETERS
guid.phase = 'ballistic';    % Initial guidance phase
guid.bankTrim = 90;          % Trim bank angle (deg)

% FNPEG targets
guid.vT = v_tgt;                                       % FNPEG target velocity (m/s)
guid.rT = planet.r + h_tgt;                            % FNPEG target radius (m). 11 km altitude
guid.eT = planet.mu/guid.rT - ...               % FNPEG target negative energy (m^2/s^2)
    0.5*guid.vT^2;    
% get target lat, lon from initial lat/lon, range target, and initial
% azimuth
[guid.lonT, guid.latT] = track_greatCirc(ICs_sph.lon0,ICs_sph.lat0,s_tgt,ICs_sph.psi0);
% FNPEG settings
guid.dpsi_i = 2;        % FNPEG lateral deadband initial width (deg)
guid.dpsi_f = 1.5;      % FNPEG lateral deadband final width (deg)
% FNPEG lateral deadband constants (deg, deg * s/m)
[guid.c0, guid.c1] = deadbConst([ICs_sph.v0, guid.dpsi_i], [guid.vT, guid.dpsi_f]);
guid.bankProfile = 'linear';  % FNPEG bank angle profile ('constant' or 'linear' with energy)
guid.bank0 = 90 * pi/180;     % FNPEG bank angle initial guess (rad)
guid.sigmaF = 60 * pi/180;    % FNPEG target bank angle for linear profile (rad)
guid.GNTol  = 1e-4;    % FNPEG Gauss-Newton exit tolerance
guid.tolInt = 1e-8;    % FNPEG predictor integration tolerance (should be several orders smaller than GNTol)
guid.dt = 1;           % FNPEG time step

%% Integration settings for the dynamics
integ.tF = 600;     % Max final time for the integration (s)
integ.tol = 1e-9;	% integration tolerance
integ.perts = true;	% Perturbations switch for EOMs


%% Launch simulation loop and save to disk

tic;
fprintf('Running FNPEG Guidance...\n');
dat = run_lme(ICs_sph, planet, vehicle, guid, integ, lRef, tRef);

s = calc_range(ICs_sph.lat0, ICs_sph.lon0, dat.lat(dat.idxend), dat.lon(dat.idxend));
ds = s - s_tgt;

elapsed = toc;
fprintf('Range Error: %g deg\n', ds);
% delete(gcp('nocreate'))


%% Plots

