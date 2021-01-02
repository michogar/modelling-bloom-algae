function ncdf = plot_longitudinal_section(nc)
%
%
    if (isstring(nc) || ischar(nc))
        ncdf = netcdf(nc);
    else
        ncdf = nc;
    end
    
    [lon, ~] = getLonLat(ncdf); x=lon(floor(end/2), :);
    h=ncdf{'h'}(floor(end/2), :); [~, m]=size(lon);
    time=ncdf{'time'}(:); [maxt, ~]=size(time);

    rgb=tempColorbar();
    colormap(rgb);

    for t=1:maxt
        zeta=ncdf{'zeta'}(t, floor(end/2), :);
        depth=squeeze(zlevs(h, zeta, 6, 0, 10, 32, 'r', 1)); x=depth;
        for i=1:m
            x(:, i)=lon(1, i);
        end
        temp=squeeze(ncdf{'temp'}(t, :, floor(end/2), :)); temp(temp==0)=NaN;
        pcolor(x, depth, temp); colorbar; title(strcat("Day", " - ", num2str(ncdf{'time'}(t:t)/86400))); pause(0.5);
    end
end