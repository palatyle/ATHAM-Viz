function ATHAM_viz_multi(tracer_name,isovalue,domain_flux,upper_dir,vent_diam,lat,density_overlay, quiver_overlay, dep_calc)
mkdir(strcat(upper_dir,'direction_out'))
for k = 1:length(lat)
    for j = 1:length(vent_diam)
        tic
        % tracer_name = 'ash1';
        % isovalue = .001;
        % domain_flux = false;
        % small_grid = false;
        % upper_dir = '/Volumes/(M)ATHAM/';
        % vent_diam = '75m';
        % lat = 'polar';
        dir_info=dir(fullfile(strcat(upper_dir,lat(k),'/',vent_diam(j)),'**/atham_netCDF_MOV.nc'));
        dirs = extractfield(dir_info,'folder')';

        output = {'Vent speed (m/s)','Wind Speed (m/s)','stability mean','stability med','stability SD','Max plume height (km)','Neutral Buoyancy Height (km)', 'NBH err (km)'};
        for i = 1:length(dirs)
            fn = dirs{i};
            if vent_diam(j) == "20m" || vent_diam(j) == "30m"
                plane_offset = 1;
            elseif vent_diam(j) == "75m" || vent_diam(j) == "127_5m"
                plane_offset = 0;
            elseif vent_diam(j) == "303m" 
                plane_offset = -1;
            end 
                
            [flux_ratio_mean, flux_ratio_med, flux_ratio_SD, max_plume_height, NBH, NBH_err,gif_str] = ATHAM_viz_ts(fn, tracer_name, isovalue, domain_flux, plane_offset, density_overlay, quiver_overlay, dep_calc);
            split_dir = regexp(fn,filesep,'split');
            output_fn = split_dir{end};

            split_fn = regexp(output_fn,'_','split');
            vent_speed = split_dir{4};

            if split_dir{3} == "127_5m"
                wind_speed = split_fn{5};
            else
                wind_speed = split_fn{4};
            end
            output{i+1,1} = vent_speed;
            output{i+1,2} = wind_speed;
            output{i+1,3} = flux_ratio_mean;
            output{i+1,4} = flux_ratio_med;
            output{i+1,5} = flux_ratio_SD;
            output{i+1,6} = max_plume_height;
            output{i+1,7} = NBH;
            output{i+1,8} = NBH_err;
            
            copyfile(strcat(gif_str,'.mat'),strcat(upper_dir,'direction_out'))
        end
        close all
        cd(upper_dir)
        writecell(output,strcat(lat(k),'_',vent_diam(j),'.txt'))
        
%         setpref('Internet','E_mail','palatyle@isu.edu');
%         setpref('Internet','SMTP_Server','smtp.gmail.com');
%         setpref('Internet','SMTP_Username','palatyle@isu.edu');
%         setpref('Internet','SMTP_Password',passwd);
%         
%         props = java.lang.System.getProperties;
%         props.setProperty('mail.smtp.auth','true');
%         props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
%         props.setProperty('mail.smtp.socketFactory.port','465');
%         
%         sendmail('palatyle@isu.edu',strcat(vent_diam(j),' Done'),'Script Finished')
        disp(strcat('Done visualizing vent ',vent_diam(j),'in:'))
        toc
    end
end

end

