%% eom_fnpeg.m
%   longitudal equations of motion for atmospheric re-entry
%       Eqs. (1,4,5,17) in 1 Lu, P., “Entry guidance: A unified Method,”
% 
function dxdt = eom_fnpeg(t, xSph, sigma0, e0, sigmaF, eF, lRef, tRef, ...
    planet, vehicle, guid, opts)
% State vector
r = xSph(1);
v = xSph(4);
gamma = xSph(5);

% Reference accel (m/s^2)
g0_ms2 = lRef / tRef^2;

% Initialize RHS
dxdt = zeros(size(xSph));

%% Kinematic equations (r, s)

dxdt(1) = v * sin(gamma);
dxdt(7) = - v * cos(gamma) / r;

% Atmospheric density
alt = (r * lRef - planet.r) * 1E-3;
dens = planet.rho0 * exp(-alt / planet.H);

% Dimensionalize velocity (m/s).
v_ms = v * lRef / tRef;

% Dimensional drag and lift acceleration magnitudes (m/s^2)
D_ms2 = 0.5 * dens * v_ms.^2 / vehicle.beta;
L_ms2 = vehicle.LD * D_ms2;

% Non-dimensionalize accelerations
L = L_ms2 / g0_ms2;
D = D_ms2 / g0_ms2;

% bank angle magnitude
if strcmpi(guid.bankProfile, 'constant')
    sigma = sigma0;
elseif strcmpi(guid.bankProfile, 'linear')
    e = 1/r - 0.5 * v^2;
    sigma = sigma0 + (sigmaF - sigma0) * (e - e0) / (eF - e0);
end
csig = cos(sigma);

% dV/dt
dxdt(4) = -D - sin(gamma)/r^2;

% d(gamma)/dt
dxdt(5) = (L*csig + (v^2 - 1/r)*(cos(gamma) / r)) / v;

end
