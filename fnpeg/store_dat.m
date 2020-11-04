%% store_dat.m
%   store lift modulation simulation parameters
% 
function dat = store_dat(last, dat, iter, tCur, xDyn, bank, direction, eCur, deadband, guid)
    dat.t(iter) = tCur;
    dat.m(iter) = xDyn.sph(8);
    dat.r(iter) = xDyn.sph(1);
    dat.lon(iter) = xDyn.sph(2) * 180/pi;
    dat.lat(iter) = xDyn.sph(3) * 180/pi;
    dat.v(iter) = xDyn.sph(4);
    dat.fpa(iter) = xDyn.sph(5) * 180/pi;
    dat.heading(iter) = xDyn.sph(6) * 180/pi;
    dat.s(iter) = xDyn.sph(7);
    
    if strcmpi(guid.phase,'entry') && ~last
        dat.e(iter) = eCur;
        dat.bank(iter) = direction * bank;
        dat.deadband(iter) = deadband * 180/pi;
    end
    
    if (last)
        dat.idxend = iter;
    end
end
