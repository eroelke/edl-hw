%% run_lme.m
%   run lift modulation problem with FNPEG guidance
% 
function [solution, controls, ie, events, navErrs, meas, atmo] = run_lme(ICs, planet, veh, guid, nav, integ, lRef, tRef)
% function out = run_lme(ICs, planet, veh, guid, nav, integ, ...
%% Initialization
% Normalize planetary angular velocity
omega = planet.Om * tRef * pi/180;

% Additional reference quantities
vRef = lRef/tRef;
gRef = vRef/tRef;
mRef = veh.m0;

% Normalize thrust and specific impulse
TMax = veh.TMax / (veh.m0 * gRef);
TMin = veh.TMin / (veh.m0 * gRef);
mDry = veh.mDry / veh.m0;
Isp  = veh.Isp / tRef;

% Initial range-to-go (deg)
sCur_deg = calc_range(ICs.lat0, ICs.lon0, guid.FNPEG.latT, guid.FNPEG.lonT);

% Pack initial state vector and normalize
x0.sph = [( ICs.h0 + planet.r )/lRef, ICs.lon0 * pi/180, ICs.lat0 * pi/180, ...
    ICs.v0 / (lRef/tRef), ICs.fpa0 * pi/180, ICs.psi0 * pi/180, sCur_deg * pi/180];
x0.cart = mixsph2cart(x0.sph(1:6));
xDyn.sph  = [x0.sph, veh.m0/mRef];    % [r, lon, lat, v, fpa, hea, range, mass]
xDyn.cart = [x0.cart, veh.m0/mRef];   % [x, y, z, vx, vy, vz, mass]
xNav = xDyn;

% Normalize final energy
eF   = guid.FNPEG.eT / (lRef/tRef)^2;

% Normalize FNPEG and landing target latitude, longitude
lonT_FNPEG = guid.FNPEG.lonT * pi/180;
latT_FNPEG = guid.FNPEG.latT * pi/180;

% Initial bank angle and direction guesses
bankGuess = guid.FNPEG.bank0;
prevDir = 1;

% Initial bank angle command
bank = guid.bankTrim * pi/180;
direction = 1;

% Set bank profile
bankProfile = guid.FNPEG.bankProfile;

% Set thrust = 0 during entry
a_Thrust = zeros(3,1);

% Normalize and initialize times and time steps
t0 = 0;
tF = integ.tF / tRef;
dt = guid.FNPEG.dt/tRef;
tCur = t0;

% Initialize outputs
nMax = 1000;
aux.cart  = zeros(nMax,8); aux.sph  = zeros(nMax,9);
aux.contr = zeros(nMax,11);
aux.nav   = zeros(nMax,7);
aux.meas  = zeros(nMax,7);
aux.atmo  = zeros(nMax,5);
% aux.pred  = cell(nMax,1);

% Initialize exit criterion, step counter
steps = 1;
ie = 0;
done = false;

% Integration options
hF = integ.hF/lRef;
opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
    'Events', @(t,x) exit_conditions(t, x, lRef, tRef, hF, eF, mDry, planet, integ));

% Events
events = struct();

%% Run Simulation Loop
while not(done)
    %% NAVIGATION

    % Add nav errors when above hF + hLim
    hCur = xDyn.sph(1) - 1;
    if hCur - hF > nav.hLim/lRef && nav.enabled
        xNav = navigation(xDyn,xNav,nav.Rbias,nav.Rvar,nav.Vbias,nav.Vvar,...
            latT_FNPEG,lonT_FNPEG,lRef,vRef);
        
    else
        xNav = xDyn;
    
    end
    
    % Current energy
    eCur = -rv2energy(1, xDyn.sph(1), xDyn.sph(4));
    
    % Predict density if desired
%     if strcmpi(guid.densType,'predicted')
%         acc = measurements(tCur, xDyn.cart, lRef, tRef, bank,...
%             guid.FNPEG.sigmaF, direction, eCur, eF, bankProfile, a_Thrust,...
%             TMax, planet, veh, guid);
%         
%         % Call to TMEP
%         [pred.logrho, pred.wE, pred.wN, featureSeq, TMEP] = ...
%             predict_TMEP_total(xNav.sph,acc,bank,featureSeq,TMEP,lRef,tRef);
%         
%     end
    
    %% GUIDANCE
    switch lower(guid.phase)
        case 'ballistic'
            % Trim bank angle
            bank = guid.bankTrim * pi/180; direction = 1;
            bankProfile = 'constant';
            
            % No thrust
            a_Thrust = zeros(3,1);
            
            % Quantities for output
            course = azimuth(xDyn.sph(3)*180/pi, xDyn.sph(2)*180/pi,...
                guid.FNPEG.latT, guid.FNPEG.lonT);
            azOff = wrapTo180(xDyn.sph(6)*180/pi - course);
            deadband = guid.FNPEG.c1 * xDyn.sph(4)*vRef + guid.FNPEG.c0;
            
            % Disable termination for final energy
            opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
                'Events', @(t,x) exit_conditions(t, x, lRef, tRef, hF, Inf,...
                mDry, planet, integ));

            % Check acceleration threshold
            acc = measurements(tCur, xDyn.cart, lRef, tRef, bank,...
            guid.FNPEG.sigmaF, direction, eCur, eF, bankProfile, a_Thrust,...
            TMax, planet, veh, guid);
            aMag_ms2 = norm(acc);
            
            if (tCur >= 0.065)
                % Ensure correct bank profile
                bankProfile = guid.FNPEG.bankProfile;
                
                % Enable termination for final energy
                opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
                'Events', @(t,x) exit_conditions(t, x, lRef, tRef, hF, eF,...
                mDry, planet, integ));
            
                % Switch to entry
                guid.phase = 'entry';
                
                continue
                
            end
        
        case 'entry'
            % Aerodynamic entry through FNPEG
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% START OF FNPEG GUIDANCE %%%

            % FNPEG longitudinal control
            [bank, k] = FNPEG_long(bankGuess, tCur, xNav.sph, tCur + tF, lRef,...
                tRef, planet, veh, guid);

            % FNPEG lateral control
            [direction, course, azOff, deadband] = FNPEG_lat(prevDir, xNav.sph, vRef, guid);
            
            % Ensure bank profile as desired
            bankProfile = guid.FNPEG.bankProfile;

            % Reset guesses for next iteration
            bankGuess = bank;
            prevDir = direction;
            
            % Enable termination for final energy
            opts = odeset('RelTol', integ.tol, 'AbsTol', integ.tol, ...
                'Events', @(t,x) exit_conditions(t, x, lRef, tRef, hF, eF,...
                mDry, planet, integ));

            %%% END OF FNPEG GUIDANCE %%%
            
    end    
    
    %% OUTPUT
    % Save output (last step = false)
    aux = saveAux(false, aux, steps, tCur, xDyn, xNav, lRef, tRef, bank, ...
    direction, eCur, eF, bankProfile, course, azOff, deadband, ...
    a_Thrust, TMax, planet, veh, guid);

    %% DYNAMICS
    
    % Propagate dynamics forward (Cartesian equations)
    [ttCur, xxDyn, ~, ~, ie] = ode45(...
    @(t, x) EOM_unified(t, x, lRef, tRef, bank, ...
    guid.FNPEG.sigmaF, direction, eCur, eF, bankProfile, a_Thrust, TMax, Isp,...
    omega, planet, veh, guid, integ.perts), [tCur, tCur + dt], ...
        xDyn.cart, opts);

    % Advance timestep
    tCur = ttCur(end);
    xDyn.cart = xxDyn(end,:);
    steps = steps + 1;
    
    % Transform to mixed-spherical coordinates
    xDyn.sph(1:6) = cart2mixsph(xDyn.cart(1:6));          % [r, lon, lat, v, fpa, hea]
    xDyn.sph(7) = distance(xDyn.sph(3), xDyn.sph(2), ...  % range-to-go
        latT_FNPEG, lonT_FNPEG,'radians');
    xDyn.sph(8) = xDyn.cart(7);                           % mass
    
    %Check exit conditions
    if ~isempty(ie)
        % Terminate integration
        break;
    end
    
end

% Save output (last step = true)
aux = saveAux(true, aux, steps, tCur, xDyn, xNav, lRef, tRef, bank, ...
    direction, eCur, eF, bankProfile, course, azOff, deadband, ...
    a_Thrust, TMax, planet, veh, guid);

%% FORMAT OUTPUT

% Cartesian trajectory
aux2 = aux.cart(1:steps,:);
solution.cart = array2table([aux2(:,1) .* tRef, aux2(:,2:4) .* lRef,...
    aux2(:,5:7) .* vRef, aux2(:,8) .* veh.m0]);
solution.cart.Properties.VariableNames = {'t', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'm'};
solution.cart.Properties.VariableUnits = {'s', 'm', 'm', 'm', 'm/s', 'm/s', 'm/s', 'kg'};

% Mixed-spherical trajectory
aux3 = aux.sph(1:steps,:);
solution.sph = array2table([aux3(:,1) .* tRef, aux3(:,2) .* lRef,...
    aux3(:,3:4) .* 180/pi, aux3(:,5) .* vRef, aux3(:,6:7) .* 180/pi, ...
    aux3(:,8) .* planet.r * 1e-3, aux3(:,8) .* veh.m0]);
solution.sph.Properties.VariableNames = {'t', 'r', 'lon', 'lat', 'v', 'fpa', 'hea', 'range', 'm'};
solution.sph.Properties.VariableUnits = {'s', 'm', 'deg', 'deg', 'm/s', 'deg', 'deg', 'km', 'kg'};

% Control history
aux4 = aux.contr(1:steps,:);
controls = array2table([aux4(:,1) .* tRef, aux4(:,2), aux4(:,3:6).*180/pi,...
    aux4(:,7).*veh.m0.*gRef, aux4(:,8).*veh.m0.*gRef, aux4(:,9:11)]);
controls.Properties.VariableNames = {'t','e','bank','course','dpsi','deadb_psi','TMag_cmd','TMag_eff','Tx','Ty','Tz'};
controls.Properties.VariableUnits = {'s', '-', 'deg', 'deg', 'deg', 'deg', 'N', 'N', '-','-','-'};

% Navigation errors
aux5 = aux.nav(1:steps,:);
navErrs = array2table([aux5(:,1) .* tRef, aux5(:,2:4) .* lRef, aux5(:,5:7) .* vRef]);
navErrs.Properties.VariableNames = {'t', 'dx', 'dy', 'dz', 'dvx', 'dvy', 'dvz'};
navErrs.Properties.VariableUnits = {'s', 'm', 'm', 'm', 'm/s', 'm/s', 'm/s'};

% Measurements
aux6 = aux.meas(1:steps,:);
meas = array2table([aux6(:,1) .* tRef, aux6(:,2:7)]);
meas.Properties.VariableNames = {'t', 'ax', 'ay', 'az', 'qInfty', 'pInfty', 'MInfty'};
meas.Properties.VariableUnits = {'s', 'm/s2', 'm/s2', 'm/s2', 'Pa', 'Pa', 'MInfty'};

% Density and winds
aux7 = aux.atmo(1:steps,:);
atmo = array2table([aux7(:,1) .* tRef, aux7(:,2:5)]);
atmo.Properties.VariableNames = {'t', 'rho', 'wE', 'wN', 'w'};
atmo.Properties.VariableUnits = {'s', 'kg/m3', 'm/s', 'm/s', 'm/s'};

% TMEP predictions
% aux8 = aux.pred(1:steps);
% predictions = {};
% predictions{1} = array2table([pred.alts; ...
%     log(10) .* (log(guid.atmo.rho0) - (pred.alts - 0)./guid.atmo.scaleH);
%     zeros(1,length(pred.alts)); zeros(1,length(pred.alts))]'); %%FIXME
% predictions{1}.Properties.VariableNames = {'h', 'logrho', 'wE', 'wN'};
% predictions{1}.Properties.VariableUnits = {'km', 'log10kgm3', 'm/s', 'm/s'};

end %run_lme.m

%% OUTPUT AUXILIARY FUNCTIONS

function aux = saveAux(last, aux, steps, tCur, xDyn, xNav, lRef, tRef, bank, ...
    direction, eCur, eF, bankProfile, course, azOff, deadband, ...
    a_Thrust, TMax, planet, veh, guid)
    %%%SAVEAUX Save output variables into 'aux' structure. Called at
    %%%each timestep.
    %
    % INPUTS:
    %   * last, logical, flag indicating whether this is the last step.
    %
    
    %%% SIMULATED MEASUREMENTS %%%
    [acc, q, pres, Mach] = measurements(tCur, xDyn.cart, lRef, tRef, bank,...
    guid.FNPEG.sigmaF, direction, eCur, eF, bankProfile, a_Thrust, TMax,...
    planet, veh, guid);
    aux.meas(steps,:) = [tCur, acc, q, pres, Mach];

    %%% TRAJECTORY %%%
    aux.cart(steps,:) = [tCur, xDyn.cart(1:6), xDyn.cart(7)];
    aux.sph(steps,:)  = [tCur, xDyn.sph];

    %%% CONTROLS %%%
    TMag_cmd = xDyn.sph(end) * norm(a_Thrust);
    TMag_eff = min([TMag_cmd,TMax]);
    if strcmpi(guid.phase,'entry') && ~last
        aux.contr(steps,:) = [tCur,eCur,direction*bank,course*pi/180,azOff*pi/180,...
        deadband*pi/180,TMag_cmd,TMag_eff,zeros(1,3)];

    else
        % During descent, zero out all FNPEG-related quantities
        aux.contr(steps,:) = [tCur,eCur,direction*bank,0,0,0,TMag_cmd,TMag_eff,zeros(1,3)];
        if TMag_eff >= 10*eps
            aux.contr(steps,:) = [tCur,eCur,direction*bank,0,0,0,TMag_cmd,TMag_eff,a_Thrust'./norm(a_Thrust)];

        end

    end

    %%% NAVIGATION %%%
    deltaR = xNav.cart(1:3) - xDyn.cart(1:3);
    deltaV = xNav.cart(4:6) - xDyn.cart(4:6);
    if ~last
        aux.nav(steps,:) = [tCur,deltaR,deltaV];

    else
        %There is no call to nav at the last step of integration,
        %therefore do not update the nav errors.
        aux.nav(steps,:) = [tCur, zeros(1,6)];

    end
    
    %%% DENSITY AND WINDS %%%
    h_km = (xDyn.sph(1) * lRef - planet.r) * 1E-3;
    dens = atm_exponential(h_km,planet.rho0,0,planet.H);
    wE = 0; wN = 0; w = 0;
    aux.atmo(steps,:) = [tCur, dens, wE, wN, w];
    
    %%% TMEP PREDICTIONS %%%
%     aux.pred{steps} = pred;
    
end

function screenOutput(last, ie, tCur, xDyn, tRef, lRef, vRef, gRef, mRef, ...
    eCur, bank, direction, azOff, a_Thrust, TMax, planet, ICs)
    %%%SCREENOUTPUT Print output to user on screen.
    %
    % INPUTS:
    %   * last, logical, flag indicating whether this is the last step.
    %
    
    TMag_cmd = xDyn.sph(end) * norm(a_Thrust);
    TMag_eff = min([TMag_cmd,TMax]);
    if ~last
        % Print output to user
        fprintf([repmat('%10.4g ',1,12),'%10.4g\n'],...
        tCur*tRef, (xDyn.sph(1)-1)*lRef, xDyn.sph(2)*180/pi,...
        xDyn.sph(3)*180/pi, xDyn.sph(4)*vRef, xDyn.sph(5)*180/pi,...
        xDyn.sph(6)*180/pi, xDyn.sph(7)*planet.r, eCur,...
        direction*bank*180/pi,azOff,TMag_eff*mRef*gRef*1e-3,...
        xDyn.sph(end)*mRef*1e-3);

    else
        % Final output to user
        fprintf([repmat('%10.4g ',1,12),'%10.4g\n'],...
                tCur*tRef, (xDyn.sph(1)-1)*lRef, xDyn.sph(2)*180/pi,...
                xDyn.sph(3)*180/pi, xDyn.sph(4)*vRef, xDyn.sph(5)*180/pi,...
                xDyn.sph(6)*180/pi, xDyn.sph(7)*planet.r, eCur,...
                direction*bank*180/pi,azOff,TMag_eff*mRef*gRef*1e-3,...
                xDyn.sph(end)*mRef*1e-3);
        if ie == 1
            fprintf('Nominal termination - reached final altitude.\n');

        elseif ie == 2
            fprintf('Nominal termination - reached final energy.\n');

        elseif ie == 3
            fprintf('Anomalous termination - reached final time.\n');

        elseif ie == 4
            fprintf('Anomalous termination - reached negative mass.\n');

        elseif ie == 5
            fprintf('Anomalous termination - NaNs in state vector.\n');

        end

        % Output useful stats
        fprintf('Range: %g deg\n',calc_range(ICs.lat0, ICs.lon0, xDyn.sph(3), sph(2)*180/pi));
        
    end

end