%% run_lme.m
%   run lift modulation problem with FNPEG guidance
% 
function dat = run_lme(ICs, planet, veh, guid, integ, lRef, tRef)
% reference quantities
vRef = lRef/tRef;

% Initial range-to-go (deg)
sCur_deg = calc_range(ICs.lat0, ICs.lon0, guid.latT, guid.lonT);

% Pack initial state vector and normalize
x0.sph = [( ICs.h0 + planet.r )/lRef, ICs.lon0 * pi/180, ICs.lat0 * pi/180, ...
    ICs.v0 / (lRef/tRef), ICs.fpa0 * pi/180, ICs.psi0 * pi/180, sCur_deg * pi/180];
x0.cart = mixsph2cart(x0.sph(1:6));
xDyn.sph  = [x0.sph 1];   % [r, lon, lat, v, fpa, heading, range mass]
xDyn.cart = [x0.cart 1];  % [x, y, z, vx, vy, vz mass]

% Normalize final energy
eF   = guid.eT / (lRef/tRef)^2;

% target to radians
lonT_FNPEG = guid.lonT * pi/180;    %target latitude, rad
latT_FNPEG = guid.latT * pi/180;    %target longitude, rad

% Initial bank angle and direction guesses
bankGuess = guid.bank0;
prevDir = 1;

% Initial bank angle command
bank = guid.bankTrim * pi/180;
direction = 1;
azErr = nan;

% Set bank profile
bankProfile = guid.bankProfile;

% Normalize and initialize times and time steps
tF = integ.tF / tRef;   %normalized final time
dt = guid.dt/tRef;      %integration step size
t = 0;  %initial time
nMax = tF/dt;   %max integration steps

% allocate output data struct
ldat1 = nan(nMax, 1);
dat.t = ldat1;
dat.r = ldat1;
dat.lon = ldat1;
dat.lat = ldat1;
dat.v = ldat1;
dat.fpa = ldat1;
dat.heading = ldat1;
dat.s = ldat1;
dat.e = ldat1;
dat.deadband = ldat1;
dat.bank = ldat1;
dat.azErr = ldat1;

% Initialize exit criterion, step counter
iter = 1;
flag = false;

% Integration options
opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
    'Events', @(t,x) check_termination(t, x, lRef, tRef, 0, eF, planet, integ));

%% Run Simulation Loop
while (~flag)
    % Current energy
    eCur = -rv2xi(1, xDyn.sph(1), xDyn.sph(4));
    
    % try do guidance
    if (mod(iter-1, guid.g_rate) == 0)
        switch lower(guid.phase)
            case 'ballistic'
                % Trim bank angle
                bank = guid.bankTrim * pi/180;
                direction = 1;
                bankProfile = 'constant';
                deadband = guid.c1 * xDyn.sph(4) * vRef + guid.c0;
                if (t >= 0.05)
                    % Ensure correct bank profile
                    bankProfile = guid.bankProfile;

                    % Enable termination for final energy
                    opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
                    'Events', @(t,x) check_termination(t, x, lRef, tRef, 0, eF, planet, integ));

                    % Switch to entry
                    guid.phase = 'entry';
                    continue;
                end

            case 'entry'
                % FNPEG longitudinal control
                bank = run_fnpeg(bankGuess, t, xDyn.sph, t + tF, lRef,...
                    tRef, planet, veh, guid);
                % check bank reversal / deadband
                [direction, deadband, azErr] = try_reversal(prevDir, xDyn.sph, vRef, guid);
                bankProfile = guid.bankProfile; % save bank profile
                bankGuess = bank;   %reset for next iter
                prevDir = direction;   %reset for next iter
                % Enable termination for final energy
                opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
                    'Events', @(t,x) check_termination(t, x, lRef, tRef, 0, eF, planet, integ));
        end %guidance
    end
    % Save output (last step = false)
    dat = store_dat(false, dat, iter, t * tRef, xDyn, bank, direction, eCur, deadband, guid, azErr);
    
    %% Propagate dynamics forward
    [ttCur, xxDyn, ~, ~, ie] = ode45(@(t, x) lme_eom(t, x, lRef, tRef, bank, ...
        guid.sigmaF, direction, eCur, eF, bankProfile, planet, veh, ... 
        guid), [t, t + dt], xDyn.cart, opts);

    % Advance timestep
    t = ttCur(end);
    xDyn.cart = xxDyn(end,:);
    iter = iter + 1;
    
    % Transform to mixed-spherical coordinates
    xDyn.sph(1:6) = cart2mixsph(xDyn.cart(1:6));          % [r, lon, lat, v, fpa, hea]
    xDyn.sph(7) = distance(xDyn.sph(3), xDyn.sph(2), latT_FNPEG, lonT_FNPEG,'radians'); %range
    xDyn.sph(8) = xDyn.cart(7); % mass
    
    %Check exit conditions
    if ~isempty(ie)
        % Terminate integration
        break;
    end
    
end

% save last data point
dat = store_dat(true, dat, iter, t * tRef, xDyn, bank, direction, eCur, deadband, guid, azErr);

end %run_lme.m

