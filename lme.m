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
vehicle.Isp   = 360;                                    % s (specific impulse)
vehicle.TMax  = 800E3;                                  % Maximum thrust magnitude (N)
vehicle.TMin  = 200E3;                                  % Minimum thrust magnitude (N)

%% Guidance parameters
s_tgt = 7.952;              % deg, range target

% GENERAL GUIDANCE PARAMETERS
guid.phase = 'ballistic';    % Initial guidance phase
guid.bankTrim = 90;          % Trim bank angle (deg)
guid.aStart = 0.15 * 9.807;  % Threshold acceleration to start guidance (m/s^2). Set = Inf to avoid starting guidance.

% Density profile used by guidance
% guid.densType = 'exp';                    % Density profile type: 'exp', 'GRAM', or 'predicted'
% guid.atmo.h0 = 0;                               % Base altitude used in 'exp' profile (km)
guid.atmo.rho0 = planet.rho0;                   % Base density used in 'exp' profile (kg/m^3)
guid.atmo.scaleH = planet.H;               % Scale height used in 'exp' profile (km)

% FNPEG targets
guid.FNPEG.vT = 1214;                                       % FNPEG target velocity (m/s)
guid.FNPEG.rT = planet.r + 11e3;                            % FNPEG target radius (m)
guid.FNPEG.eT = planet.mu/guid.FNPEG.rT - ...               % FNPEG target negative energy (m^2/s^2)
    0.5*guid.FNPEG.vT^2;    
[guid.FNPEG.lonT, guid.FNPEG.latT] =...                     % FNPEG target longitude and latitude (deg)
    track_greatCirc(ICs_sph.lon0,ICs_sph.lat0,s_tgt,ICs_sph.psi0);
% FNPEG settings
guid.FNPEG.dpsi_i = 2;                                      % FNPEG lateral deadband initial width (deg)
guid.FNPEG.dpsi_f = 1.5;                                    % FNPEG lateral deadband final width (deg)
[guid.FNPEG.c0, guid.FNPEG.c1] = ...                        % FNPEG lateral deadband constants (deg, deg * s/m)
    deadbConst([ICs_sph.v0, guid.FNPEG.dpsi_i],...
    [guid.FNPEG.vT, guid.FNPEG.dpsi_f]);
guid.FNPEG.bankProfile = 'linear';                          % FNPEG bank angle profile ('constant' or 'linear' with energy)
guid.FNPEG.bank0 = 90 * pi/180;                             % FNPEG bank angle initial guess (rad)
guid.FNPEG.sigmaF = 60 * pi/180;                            % FNPEG target bank angle for linear profile (rad)
guid.FNPEG.GNTol  = 1e-4;                                   % FNPEG Gauss-Newton exit tolerance
guid.FNPEG.tolInt = 1e-8;                                   % FNPEG predictor integration tolerance (should be several orders smaller than GNTol)
guid.FNPEG.dt = 1;                                          % FNPEG time step

%% Navigation parameters
% Navigation biases and variances
nav.enabled = false;                                        % = true to introduce errors, = false for perfect navigation
nav.Rbias = 0;                                              % Navigation - position bias (m)
nav.Vbias = 0;                                              % Navigation - velocity bias (m)
nav.Rvar  = 0^2;                                            % Navigation - position variance (m^2)
nav.Vvar  = 0^2;                                            % Navigation - velocity variance (m^2/s^2)
nav.hLim  = 30;                                             % Navigation - threshold altitude for error nulling

%% Integration settings for the dynamics
integ.hF = 0;                                               % Minimum final altitude (m)
integ.tF = 600;                                             % Max final time for the integration (s)
integ.tol = 1e-9;                                           % EOM integration tolerance
integ.perts = true;                                         % Perturbations switch for EOMs

%% Load network for density predictions

TMEP = struct;
TMEP.alts = 5:2:80;


%% Launch simulation loop and save to disk

tic;
fprintf('Running FNPEG Guidance...\n');
[solution, controls, flag, events, navErrs, meas, atmo] =...
    run_lme(ICs_sph, planet, vehicle, guid,...
    nav, integ, lRef, tRef);

s = calc_range(ICs_sph.lat0, ICs_sph.lon0, solution.sph.lat(end), solution.sph.lon(end));
ds = s - s_tgt;

elapsed = toc;
fprintf(['Finished! Time elapsed: %g s.\n' ... 
    'Range: %g deg\n' ... 
    'Range Error: %g deg\n'],elapsed, s, ds);
delete(gcp('nocreate'))


%% Plots

