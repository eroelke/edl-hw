    function [sigmaCur, k] = FNPEG_long(sigmaGuess, tCur, xCurSph, tF, lRef,...
    tRef, planet, vehicle, guid)
%%FNPEG_long FNPEG longitudinal guidance
%
%%TODO
% All quantities are dimensionless.
%
% AUTHOR:
%   Davide Amato, CU Boulder, davide.amato@colorado.edu
% 
%% Initializations

sigmaCur = sigmaGuess;
sigmaF   = guid.FNPEG.sigmaF;
epsSigma = 1 * pi/180;           % Epsilon to compute d(f)/d(\sigma) at zero-th/first iterations
epsExit  = guid.FNPEG.GNTol;     % Tolerance for the Gauss-Newton exit criterion
tolInt   = guid.FNPEG.tolInt;    % Integration tolerance

% Initialize energy, time
eCur  = -rv2energy(1,xCurSph(1),xCurSph(4));
eF    = guid.FNPEG.eT / (lRef/tRef)^2;
t0    = tCur;

opts = odeset('AbsTol',tolInt,'RelTol',tolInt,'Events',@(t, xSph) final_energy(t,xSph,eF));

%% First iteration

% First iteration. Propagate trajectory once and compute error function
[~, xxSph] = ode45(@(t, xSph) EOM_FNPEG_long(t, xSph, sigmaCur, eCur,...
    sigmaF, eF, lRef, tRef, planet, vehicle, guid), [t0 tF], xCurSph, opts);
check_eF(xxSph(end,1),xxSph(end,4))
z1 = xxSph(end,7);
f1 = 0.5 * z1^2;

% Perturb initial \sigma_0 to compute finite difference
[~, xxSph] = ode45(@(t, xSph) EOM_FNPEG_long(t, xSph, sigmaCur + epsSigma, eCur,...
    sigmaF, eF, lRef, tRef, planet, vehicle, guid), [t0 tF], xCurSph, opts);
check_eF(xxSph(end,1),xxSph(end,4))
z1Pert = xxSph(end,7);
dzds = (z1Pert - z1) / epsSigma;

% If dz/ds = 0, we're likely close to the end of entry and we can assume
% sigma_0 = the initial guess
if abs(dzds) <= eps
    k = 0;
    return
    
end

% Initialize exit criterion, counter
exitCrit = false;
k = 1;

% Set maximum iterations
kMax = 15;

% Initial value of residual, error function
z = z1;
f = f1;

%% Gauss-Newton iterations

while(not(exitCrit) && (k <= kMax))
    % Save past values of bank angle, residual, error function
    sigmaPastIter = sigmaCur; zPastIter = z; fPastIter = f;
    
    % Initialize step size, inner counter
    lambda = 1;
    j = 1;
    
    while ((j == 1 || f > fPastIter) && j <= 10)
        % Compute value of the current bank angle using eq. (25). The sign of
        % sigma should be decided by the bank angle logic.
        sigmaCur = sigmaCur - lambda * z / dzds;
        sigmaCur = mod(sigmaCur, 2*pi);
        
        % Compute residual, error function, dzds
        [~, xxSph] = ode45(@(t, xSph) EOM_FNPEG_long(t, xSph, sigmaCur, eCur,...
            sigmaF, eF, lRef, tRef, planet, vehicle, guid), [t0 tF], xCurSph, opts);
        check_eF(xxSph(end,1),xxSph(end,4))
        z = xxSph(end,7);
        f = 0.5 * z^2;
        dzds = (z - zPastIter) / (sigmaCur - sigmaPastIter);
        
        % Reduce step-size and increase counter
        lambda = 0.5 * lambda;
        j = j + 1;
        
    end
    % Check exit condition
    exitCrit = abs(z * dzds) <= epsExit;
    
    % Advance counter
    k = k + 1;
end

% Output magnitude only
sigmaCur = abs(sigmaCur);

    function check_eF(rF,vF)
        eF_reached = -rv2energy(1,rF,vF);
        if (abs(eF_reached - eF) > eps * 1e2)
            message = sprintf(['FNPEG_long: At tCur = %g s, the final energy reached '...
                'after the integration of the EoMs\n is different than the prescribed '...
                'final energy. If this is the first\n time step, try increasing the initial '...
                'bank angle guess.'],tCur*tRef);
            warning(message);
        end
        
    end

end