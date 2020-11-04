%% run_lme.m
%   run lift modulation problem with FNPEG guidance
% 
function dat = run_lme(ICs, planet, veh, guid, integ, lRef, tRef)
% function out = run_lme(ICs, planet, veh, guid, nav, integ, ...
%% Initialization
% Normalize planetary angular velocity
omega = planet.Om * tRef * pi/180;

% Additional reference quantities
vRef = lRef/tRef;
gRef = vRef/tRef;
mRef = veh.m0;

% Normalize thrust and specific impulse
% TMax = veh.TMax / (veh.m0 * gRef);
mDry = veh.mDry / veh.m0;

% Initial range-to-go (deg)
sCur_deg = calc_range(ICs.lat0, ICs.lon0, guid.latT, guid.lonT);

% Pack initial state vector and normalize
x0.sph = [( ICs.h0 + planet.r )/lRef, ICs.lon0 * pi/180, ICs.lat0 * pi/180, ...
    ICs.v0 / (lRef/tRef), ICs.fpa0 * pi/180, ICs.psi0 * pi/180, sCur_deg * pi/180];
x0.cart = mixsph2cart(x0.sph(1:6));
xDyn.sph  = [x0.sph, veh.m0/mRef];    % [r, lon, lat, v, fpa, hea, range, mass]
xDyn.cart = [x0.cart, veh.m0/mRef];   % [x, y, z, vx, vy, vz, mass]
xNav = xDyn;

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

% Set bank profile
bankProfile = guid.bankProfile;

% Normalize and initialize times and time steps
tF = integ.tF / tRef;
dt = guid.dt/tRef;
tCur = 0;

% Initialize outputs
nMax = 1000;
ldat1 = nan(nMax, 1);
dat.t = ldat1;
dat.m = ldat1;
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

% Initialize exit criterion, step counter
iter = 1;
% ie = 0;
done = false;

% Integration options
opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
    'Events', @(t,x) exit_conditions(t, x, lRef, tRef, 0, eF, mDry, planet, integ));

%% Run Simulation Loop
while (~done)
    %% NAVIGATION
    xNav = xDyn;
    
    % Current energy
    eCur = -rv2xi(1, xDyn.sph(1), xDyn.sph(4));
    
    %% GUIDANCE
    switch lower(guid.phase)
        case 'ballistic'
            % Trim bank angle
            bank = guid.bankTrim * pi/180; direction = 1;
            bankProfile = 'constant';
%             
%             % No thrust
%             a_Thrust = zeros(3,1);
%             
%             % Quantities for output
            course = azimuth(xDyn.sph(3)*180/pi, xDyn.sph(2)*180/pi,...
                guid.latT, guid.lonT);
            azOff = wrapTo180(xDyn.sph(6)*180/pi - course);
            deadband = guid.c1 * xDyn.sph(4)*vRef + guid.c0;
            if (tCur >= 0.065)
                % Ensure correct bank profile
                bankProfile = guid.bankProfile;
                
                % Enable termination for final energy
                opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
                'Events', @(t,x) exit_conditions(t, x, lRef, tRef, 0, eF,...
                mDry, planet, integ));
            
                % Switch to entry
                guid.phase = 'entry';
                continue
                
            end
        
        case 'entry'
            % FNPEG longitudinal control
            bank = run_fnpeg(bankGuess, tCur, xNav.sph, tCur + tF, lRef,...
                tRef, planet, veh, guid);
            % check bank reversal / deadband
            [direction, deadband] = try_reversal(prevDir, xNav.sph, vRef, guid);
            % save bank profile
            bankProfile = guid.bankProfile;
            % Reset guesses for next iteration
            bankGuess = bank;
            prevDir = direction;
            % Enable termination for final energy
            opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
                'Events', @(t,x) exit_conditions(t, x, lRef, tRef, 0, eF,...
                mDry, planet, integ));            
    end    
    % Save output (last step = false)
    dat = store_dat(false, dat, iter, tCur, xDyn, bank, direction, eCur, deadband, guid);
    
    %% Propagate dynamics forward
    [ttCur, xxDyn, ~, ~, ie] = ode45(...
    @(t, x) lme_eom(t, x, lRef, tRef, bank, ...
    guid.sigmaF, direction, eCur, eF, bankProfile,...
    omega, planet, veh, guid, integ.perts), [tCur, tCur + dt], ...
        xDyn.cart, opts);

    % Advance timestep
    tCur = ttCur(end);
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
dat = store_dat(true, dat, iter, tCur, xDyn, bank, direction, eCur, deadband, guid);

end %run_lme.m

