function dxSphdt = EOM_FNPEG_long(t, xSph, sigma0, e0, sigmaF, eF, lRef, tRef, ...
    planet, vehicle, guid, dummy)
%%EOM_FNPEG_LONG Right-hand side of equations of motion for longitudinal
%%atmospheric entry dynamics for FNPEG.
%
% This function computes the value of the right-hand side of the equations
% of motion for the atmospheric entry of a vehicle, expressed as a function
% of time. *Only longitudinal (along-track) dynamics* are computed. In 
% particular, we only integrate equations for r, V, gamma, s 
% (Eqs. (1,4,5,17) in Ref. [1]).
% The planetary rotation rate and the offset between the heading and
% azimuth of the target site are neglected (Omega = 0, \Delta{\psi} = 0).
% The independent variable is time.
%
% All quantities used in the computation of the right-hand side are
% dimensionless.
% 
% AUTHOR:
%   Davide Amato, CU Boulder, davide.amato@colorado.edu
% 
%% Unpack, non-dimensionalize, create auxiliary variables
% State vector
r = xSph(1); v = xSph(4); gamma = xSph(5);

% Auxiliary variables
sgam = sin(gamma); cgam = cos(gamma);

% Reference acceleration (m/s^2)
g0_ms2 = lRef / tRef^2;

% Initialize RHS
dxSphdt = zeros(size(xSph));

%% Kinematic equations (r, s)

dxSphdt(1) = v * sgam;
dxSphdt(7) = - v * cgam / r;

%% Atmospheric density

% Height (km)
height = (r * lRef - planet.r)*1E-3;

% Density profile is exponential, with constants given by user.
% This case corresponds to "imperfect knowledge" of density.
dens = atm_exponential(height, guid.atmo.rho0,0, guid.atmo.scaleH);

%% Aerodynamic accelerations
% Dimensionalize velocity (m/s).
v_ms = v * lRef / tRef;

% Dimensional drag and lift acceleration magnitudes (m/s^2)
D_ms2 = 0.5 * dens * v_ms.^2 / vehicle.B0;
L_ms2 = vehicle.LD * D_ms2;

% Non-dimensionalize accelerations
L = L_ms2 / g0_ms2;
D = D_ms2 / g0_ms2;

%% Bank angle magnitude
if strcmpi(guid.FNPEG.bankProfile, 'constant')
    sigma = sigma0;
    
elseif strcmpi(guid.FNPEG.bankProfile, 'linear')
    e = 1/r - 0.5 * v^2;
    sigma = sigma0 + (sigmaF - sigma0) * (e - e0) / (eF - e0);

end
csig = cos(sigma);

%% Dynamic equations (V, gamma)

% dV/dt
dxSphdt(4) = -D - sgam/r^2;

% d(gamma)/dt
dxSphdt(5) = ( L * csig + ( v^2 - 1/r ) * ( cgam / r ) ) / v;

end