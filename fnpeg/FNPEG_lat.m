function [direction, course, azOff, deadband] = FNPEG_lat(prevDir, xCur, vRef, guid)
%%FNPEG_LAT FNPEG lateral guidance
%
%%TODO
%
% AUTHOR:
%   Davide Amato, CU Boulder, davide.amato@colorado.edu
% 
%% Unpack
theta_deg = xCur(2) * 180/pi; phi_deg = xCur(3) * 180/pi;
v_kms = xCur(4) * vRef; psi_deg = xCur(6) * 180/pi;

%% Sign rule for bank angle

% Velocity-dependent deadband (deg)
deadband = guid.FNPEG.c1 * v_kms + guid.FNPEG.c0;

% Azimuth of target site and offset (deg)
course   = azimuth(phi_deg, theta_deg, guid.FNPEG.latT, guid.FNPEG.lonT);
% Shortest arc between psi_deg and azT
azOff = wrapTo180(psi_deg - course);

if (abs(azOff) <= deadband)
    % Within the deadband, keep previous sign
    direction = prevDir;
    
else
    % Outside the deadband, change sign to minimize offset
    direction = -sign(azOff);
    
end

end