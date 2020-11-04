%%EOM_UNIFIED Equations of motion in Cartesian coordinates for entry
%%and powered descent with prescribed thrust vector and bank angle.
% Compute the right-hand side of the equations of motion in Cartesian
% coordinates in the powered descent [1]. Aerodynamic forces are included
% following the approach in [2].
% 
% Thrust is included according to either of the following 'mode's:
%    * mode = 'const', constant thrust is applied in a fixed direction in
%    the wind reference frame. This also includes the case of no thrust
%    (simply null the magnitude of the thrust vector).
%    * mode = 'guid', thrust direction is the same as the time-varying 
%    velocity costate [1].
% 
% INPUTS: %%FIXME
%    * t, double, dimensionless time.
%    * x, 7x1 double, dimensionless state and mass.
%    * mode, string, either 'const' or 'guid' according to the thrust mode.
%    * planet, struct, planetary data.
%    * vehicle, struct, vehicle parameters.
%    * omega, double, dimensionless planetary angular velocity.
%    * lRef, double, reference length (m). Used to compute dimensional
%      altitude for drag.
%    * tRef, double, reference time (s). Used to compute dimensional
%      aerodynamic accelerations.
%    * t1, double, first dimensionless switch time
%    * t2, double, second dimensionless switch time
%    * t0, double, reference time for the costate
%    * lambda0, 6x1, costate at t = t0
%    * thrustParams, struct, thrust parameters. This includes:
%          * t0, double, dimensionless time to which the guidance solution
%          is referred. Only used in 'guid' mode.
%          * lambda0, 6x1 double, dimensionless costate at time t0. Only
%          used in 'guid' mode.
%          * TMax, double, dimensionless maximum thrust magnitude. Used
%          in both 'guid' and 'const' modes.
%          * TMin, double, dimensionless minimum thrust magnitude. Used
%          in both 'guid' and 'const' modes.
%          * Isp, double, dimensionless specific impulse. Used in both
%          'guid' and 'const' modes.
%          * uWind, 3x1 double, thrust vector direction in wind frame. Used
%          only in 'const' mode.
%     * perts, logical, = true if aero and fictitious accelerations are
%     considered
% 
function dxdt = lme_eom(t, x, lRef, tRef, sigma0, sigmaF, direction, e0, eF, ...
    bankProfile, omega, planet, vehicle, guid, perts)

% extract state
r = x(1:3);
v = x(4:6);
m = x(7);

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

% Altitude and time (km, s)
h_km = (rMag * lRef - planet.r) * 1E-3;
t_sec = t * tRef;

rho = atm_exponential(h_km, planet.rho0, 0, planet.H);
w = zeros(3,1);

% Velocity wrt flow
vInfty = v - w;

%% Wind unit vectors
f = zeros(3,1);
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
    L_ms2 = vehicle.LD * norm(D_ms2) .* (cos(sigma) .* xW - sin(sigma) .* zW);

    % Drag and lift acceleration and aerodynamic force (dimensionless)
    D = D_ms2 / (lRef / tRef^2);
    L = L_ms2 / (lRef / tRef^2);
    f = D + L;

end

dxdt = zeros(6,1);
% Auxiliary quantities
r3          = rMag^3;
grav        = -r/r3;
coriolis    = - (2 * omega) .* [-v(2); v(1); 0];
centrifugal = - omega^2     .* [r(1); r(2); 0];

if ~perts
    % Null perturbing accelerations, linearize gravity
    grav = -r;
    f = 0;
    coriolis = 0;
    centrifugal = 0;
end

dxdt(1:3) = v;
dxdt(4:6) = grav + f + coriolis + centrifugal;  %F = ma
dxdt(7) = 0;    %delta mass = 0

end