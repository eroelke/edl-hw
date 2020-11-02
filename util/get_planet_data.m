function planet = get_planet_data(name)

switch (name)
    case 'Earth'
        planet.H = 8500;  %atm scale height, m
        planet.r = 6.3568e6;  %radius, m
        planet.mu = 3.9860e14; %gravitational param, m^3/s^2
        planet.rho0 = 1.217;  %surface density, kg/m3
        planet.g = 9.798;   %grav accel, m/s2
        planet.hatm = 125;  %altitude of atmospheric interface
    case 'Mars'
        planet.mu   = 4.305E13;                      % m^3/s^-2 (gravitational parameter)
        planet.r    = 3397.2E3;                      % m (mean radius)
        planet.TRot = 1.02595675;                    % day (rotation period)
        planet.Om   = 360 / (planet.TRot * 86400);   % deg/s (angular velocity)
        planet.g    = planet.mu / planet.r^2;        % m/s^2 (surface acceleration)
        planet.gammaR = 1.2 * 191;                   % J/(kg K) (specific heat ratio * gas constant)
        planet.rho0 = 0.0263;
        planet.H = 10.1536;
        planet.hatm = 100;  %altitude of atmospheric interface
    otherwise
        error('Bad Planet');
end
end