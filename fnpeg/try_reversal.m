%% try_reversal.m
%   logic for FNPEG bank reversals
% 
function [direction, deadband, azErr] = try_reversal(prevDir, xCur, vRef, guid)
theta_deg = xCur(2) * 180/pi; phi_deg = xCur(3) * 180/pi;
v_kms = xCur(4) * vRef; psi_deg = xCur(6) * 180/pi;

% Velocity-dependent deadband (deg)
deadband = guid.c1 * v_kms + guid.c0;

% Azimuth of target site and offset (deg)
course  = azimuth(phi_deg, theta_deg, guid.latT, guid.lonT);
% Shortest arc between psi_deg and azT
azErr = wrapTo180(psi_deg - course);

if (abs(azErr) <= deadband)
    % Within the deadband, keep previous sign
    direction = prevDir;
else
    % Outside the deadband, change sign to minimize offset
    direction = -sign(azErr);
end

end