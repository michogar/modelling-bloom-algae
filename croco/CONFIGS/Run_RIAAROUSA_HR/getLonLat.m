function [lon, lat] = getLonLat(nc)
% Returns lon and lat from xi and eta parameters
% Accepts string with the netcdf path or netcdf

    if (isstring(nc) || ischar(nc))
        ncdf = netcdf(nc);
    else
        ncdf = nc;
    end
    eta = ncdf{'eta_rho'}(:);
    xi = ncdf{'xi_rho'}(:);
    lon = ncdf{'lon_rho'}(eta, xi);
    lat = ncdf{'lat_rho'}(eta, xi);
end

