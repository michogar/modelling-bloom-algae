function ncdf = plot_surface(nc)
% Plots the surface of the netcdf
% Accepts string with the netcdf path or netcdf

    if (isstring(nc) || ischar(nc))
        ncdf = netcdf(nc);
    else
        ncdf = nc;
    end
    [lon, lat] = getLonLat(nc);
    h = ncdf{'h'}(:); h(h<100)=NaN;
    pcolor(lon, lat, -h); shading flat;
end