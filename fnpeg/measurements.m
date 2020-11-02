function [a,q,p,M] = measurements(t, x, lRef, tRef, sigma0, sigmaF, direction,...
    e0, eF, bankProfile, a_Thrust, TMax, planet, vehicle, guid)
%%%MEASUREMENTS Simulate perfect measurements of acceleration and air data
% 
% The acceleration measured by the vehicle is a_meas = f_aero +
% Thrust. Most of this function is taken from EOM_unified. The measured
% acceleration is expressed along the axes of the MCMF frame (equivalently,
% the axes of the rest frame in which the acceleration is expressed are 
% along the same direction as those of the MCMF frame).
%
% The other measurements consists of air data, in particular: freestream
% dynamic pressure q, freestream static pressure p, and Mach number M.
% The measurements are not affected by noise.
%
% Input units are equal to those in EOM_unified (t, x dimensionless).
% Outputs are dimensional.
% 
% AUTHOR:
%    Davide Amato, CU Boulder, davide.amato@colorado.edu.
%
%% Unpack inputs
r = x(1:3)'; v = x(4:6)'; m = x(7);

rMag = norm(x(1:3));
vMag = norm(x(4:6));

%% Bank angle magnitude and sign
if (strcmpi(bankProfile, 'constant') || strcmpi(guid.phase, 'ballistic'))
    sigma = sigma0;
    
elseif (strcmpi(bankProfile, 'linear'))
    e = 1/rMag - 0.5 * vMag^2;
    sigma = sigma0 + (sigmaF - sigma0) * (e - e0) / (eF - e0);

end

sigma = direction * sigma;
csig  = cos(sigma); ssig = sin(sigma);

%% Get density and winds

% Altitude and time (km, s)
h_km = (rMag * lRef - planet.r)*1E-3;
t_sec = t * tRef;

rho = atm_exponential(h_km, planet.rho0, 0, planet.H);
w = zeros(3,1);

% Velocity wrt flow
vInfty = v - w;

%% Wind unit vectors
f = zeros(3,1); vInfty_msMag = 0;
if vMag > 10*eps  % Neglect aerodynamics at very small velocity
    
    yW = vInfty ./ norm(vInfty);
    zW = cross(r,vInfty); zW = zW ./ norm(zW);
    xW = cross(yW, zW); % Check: norm(xW) = 1

    %% Aerodynamic accelerations
    % This section is in dimensional units.

    % Ballistic parameter (kg/m^2). Remember that 'm' is m(t)/m0.
    B = m * vehicle.B0;

    % Drag acceleration (m/s^2)
    vInfty_ms  = vInfty .* (lRef / tRef);
    vInfty_msMag = norm(vInfty_ms);
    D_ms2 = (-0.5 * (rho/B) * vInfty_msMag) .* vInfty_ms;

    % Lift acceleration (m/s^2)
    L_ms2 = vehicle.LD * norm(D_ms2) .* (csig .* xW - ssig .* zW);

    % Drag and lift acceleration and aerodynamic force (dimensionless)
    D = D_ms2 / (lRef / tRef^2);
    L = L_ms2 / (lRef / tRef^2);
    f = D + L;

end
%% Thrust

% Consider thrust saturation
TMag = m * norm(a_Thrust);
if TMag > TMax
    TMag = TMax;
    a_Thrust = TMag .* a_Thrust ./ norm(a_Thrust);
    
end

%% Generate measurements

% Acceleration (m/s^2)
a = (a_Thrust + f)'; a = a * (lRef/tRef^2);

p = 0; q = 0; M = 0;
end