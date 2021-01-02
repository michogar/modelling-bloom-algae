function ncdf = plot_latitudinal_section(nc)
%
%
    if (isstring(nc) || ischar(nc))
        ncdf = netcdf(nc);
    else
        ncdf = nc;
    end
    
    [~, lat] = getLonLat(ncdf);
    h=ncdf{'h'}(:, floor(end/2)); [m, ~]=size(lat);
    time=ncdf{'time'}(:); [maxt, ~]=size(time);

    rgb=tempColorbar();
    colormap(rgb);

    for t=1:maxt
        zeta=ncdf{'zeta'}(t, :, floor(end/2));
        depth=squeeze(zlevs(h, zeta, 6, 0, 10, 32, 'r', 1)); x=depth;

        temp=squeeze(ncdf{'temp'}(t, :, :, floor(end/2))); temp(temp<-10)=NaN;
        for i=1:m
            y(i, :)=lat(i, 1);
        end
        pcolor(-y, depth, temp); colorbar; shading flat; title(strcat("Day", " - ", num2str(ncdf{'time'}(t:t)/86400))); pause(0.5);
    end
end