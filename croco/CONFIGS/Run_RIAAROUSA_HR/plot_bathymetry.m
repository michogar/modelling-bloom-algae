function ncdf = plot_bathymetry(nc)
% Plots the bathymetry of the netcdf
% Accepts string with the netcdf path or netcdf

    if (isstring(nc) || ischar(nc))
        ncdf = netcdf(nc);
    else
        ncdf = nc;
    end
    [lon, lat] = getLonLat(nc);
    h = ncdf{'h'}(:);
    surf(lon, lat, -h);
end