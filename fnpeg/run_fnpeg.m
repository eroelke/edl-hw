function sigmaCur = run_fnpeg(sigmaGuess, tCur, xCurSph, tF, lRef,...
    tRef, planet, vehicle, guid)

sigmaCur = sigmaGuess;
sigmaF   = guid.sigmaF;
epsSigma = 1 * pi/180;           % Epsilon to compute d(f)/d(\sigma) at zero-th/first iterations
epsExit  = guid.GNTol;     % Tolerance for the Gauss-Newton exit criterion
tolInt   = guid.tolInt;    % Integration tolerance

% Initialize energy, time
eCur  = -rv2xi(1,xCurSph(1),xCurSph(4));
eF    = guid.eT / (lRef/tRef)^2;
t0    = tCur;

opts = odeset('AbsTol',tolInt,'RelTol',tolInt,'Events',@(t, xSph) final_energy(t,xSph,eF));

%% First iteration
% First iteration. Propagate trajectory once and compute error function
[~, xxSph] = ode45(@(t, xSph) eom_fnpeg(t, xSph, sigmaCur, eCur,...
    sigmaF, eF, lRef, tRef, planet, vehicle, guid), [t0 tF], xCurSph, opts);
% check_eF(xxSph(end,1),xxSph(end,4))
z1 = xxSph(end,7);  %range, s1
f1 = 0.5 * z1^2;    %f(s1) = s1^2/2;

% Perturb initial \sigma_0 to compute finite difference
[~, xxSph] = ode45(@(t, xSph) eom_fnpeg(t, xSph, sigmaCur + epsSigma, eCur,...
    sigmaF, eF, lRef, tRef, planet, vehicle, guid), [t0 tF], xCurSph, opts);
% check_eF(xxSph(end,1),xxSph(end,4))
z1Pert = xxSph(end,7); %range, sPert
dzds = (z1Pert - z1) / epsSigma;

% If dz/ds = 0, we're likely close to the end of entry and we can assume
% sigma_0 = the initial guess
% if abs(dzds) <= eps
%     return
% end

% Initialize exit criterion, counter
exitFlag = false;
iter = 1;

% Set maximum iterations
iter_max = 15;

% Initial value of residual, error function
z = z1;
f = f1;

%% Gauss-Newton iterations
while(~exitFlag && iter <= iter_max)
    % Save past values of bank angle, residual, error function
    sigmaPastIter = sigmaCur;
    zPastIter = z;
    fPastIter = f;
    
    % Initialize step size, inner counter
    lambda = 1;
    j = 1;
%     while ((j == 1 || f > fPastIter) && j <= 10)
        % Compute value of the current bank angle using eq. (25). The sign of
        % sigma should be decided by the bank angle logic.
        sigmaCur = sigmaCur - lambda * z / dzds;
        sigmaCur = mod(sigmaCur, 2*pi);
        
        % Compute residual, error function, dzds
        [~, xxSph] = ode45(@(t, xSph) eom_fnpeg(t, xSph, sigmaCur, eCur,...
            sigmaF, eF, lRef, tRef, planet, vehicle, guid), [t0 tF], xCurSph, opts);
%         check_eF(xxSph(end,1),xxSph(end,4))
        z = xxSph(end,7);
        f = 0.5 * z^2;
        dzds = (z - zPastIter) / (sigmaCur - sigmaPastIter);
        
        % Reduce step-size and increase counter
        lambda = 0.5 * lambda;
        j = j + 1;
%     end
    exitFlag = abs(z * dzds) <= epsExit;    %exit conditions
    iter = iter + 1;    %increment loop counter
end

% Output magnitude only
sigmaCur = abs(sigmaCur);

    function check_eF(rF,vF)
        eF_reached = -rv2xi(1,rF,vF);
        if (abs(eF_reached - eF) > eps * 1e2)
            message = sprintf(['FNPEG_long: At tCur = %g s, the final energy reached '...
                'after the integration of the EoMs\n is different than the prescribed '...
                'final energy. If this is the first\n time step, try increasing the initial '...
                'bank angle guess.'],tCur*tRef);
            warning(message);
        end
        
    end

end