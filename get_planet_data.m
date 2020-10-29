function [H, rP, mu, rho0, g, h0] = get_planet_data(planet)

switch (planet)
    case 'Venus'
        H = 15900;
        rP = 6.0518e6;
        mu = 3.24858592079e14;
        g = 8.870;
        rho0 = 65;
        h0 = 150;
    case 'Earth'
        H = 8500;  %atm scale height, m
        rP = 6.3568e6;  %radius, m
        mu = 3.9860e14; %gravitational param, m^3/s^2
        rho0 = 1.217;  %surface density, kg/m3
        g = 9.798;   %grav accel, m/s2
        h0 = 125;
    case 'Mars'
        rP = 3.3962e6;
        mu = 4.2830e13;
        H = 11100;
        rho0 = 0.02;
        g = 3.710;
        h0 = 100;
    otherwise
        error('Bad Planet');
end
end