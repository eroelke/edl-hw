%% lme_eom.m
%   run dynamics for lifting entry
% 
function dxdt = lme_eom(t, x, lRef, tRef, sigma0, sigmaF, direction, e0, eF, ...
    bankProfile, planet, vehicle, guid)

% extract state
r = x(1:3);
v = x(4:6);
% m = x(7);

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

alt = (rMag * lRef - planet.r) / 1000;
rho = planet.rho0 * exp(-alt / planet.H);
w = zeros(3,1);
% Velocity wrt flow
vInfty = v - w;
yW = vInfty ./ norm(vInfty);    %unit y
zW = cross(r,vInfty);
zW = zW ./ norm(zW);    %unit z
xW = cross(yW, zW); % unit x

% Drag acceleration (m/s^2)
vInfty_ms  = vInfty .* (lRef / tRef);
vInfty_msMag = norm(vInfty_ms);
D_ms2 = (-0.5 * (rho/vehicle.beta) * vInfty_msMag) .* vInfty_ms;

% Lift acceleration (m/s^2)
L_ms2 = vehicle.LD * norm(D_ms2) .* (cos(sigma) .* xW - sin(sigma) .* zW);

% Drag and lift acceleration and aerodynamic force (dimensionless)
D = D_ms2 / (lRef / tRef^2);
L = L_ms2 / (lRef / tRef^2);
f = D + L;

dxdt = zeros(6,1);
grav = -r/rMag^3;

dxdt(1:3) = v;
dxdt(4:6) = grav + f;
dxdt(7) = 0;
end