function ncdf = plot_surface_to_gif(nc)
%
%
    if (isstring(nc) || ischar(nc))
        ncdf = netcdf(nc);
    else
        ncdf = nc;
    end
    
    [lon, lat] = getLonLat(ncdf);
    
    time=ncdf{'time'}(:); [maxt, ~]=size(time);

    rgb=tempColorbar();
    colormap(rgb);

    initLoop=1;

    for t=initLoop:maxt

        temp=ncdf{'temp'}(t, 32, :, :); temp(temp==0)=NaN;
        pcolor(lon, lat, temp); colorbar; shading flat; title(strcat("Day", " - ", num2str(ncdf{'time'}(t:t)/86400), " - ",  "temp ºC")); pause(0.2);

        % gif utilities
        set(gcf,'color','w'); % set figure background to white
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        outfile = 'surface.gif';

        % On the first loop, create the file. In subsequent loops, append.
        if t==initLoop
            imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
        else
            imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
        end

    end
    
end