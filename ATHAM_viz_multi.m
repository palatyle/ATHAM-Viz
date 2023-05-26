function ATHAM_viz_multi(tracer_name,isovalue,domain_flux,upper_dir,vent_diam,lat,density_overlay, quiver_overlay, dep_calc)
    % ATHAM_viz_muilt.m  runs main visualization code for multiple runs depending on latitude and vent diameter

    % INPUTS
    % Tracer_name is the name of the tracer you want to visualize (":ash1":, "ash2", etc)
    % isovalue is a vector of the isovalues you want to visualize
    % domain_flux is a boolean that determines whether or not you want to visualize the flux through the domain edges of the simulation
    % upper_dir is the directory that contains the runs you want to visualize. This code assumes the directory structure is as such:
        % upper_dir:
        %     lat:
        %         vent_diam:
        %             vent_velocity:
        %                 lat_vent_diam_vent_velocity_wind_speed:
        %                      atham_netCDF_MOV.nc
    % vent_diam is a vector of the vent diameters you want to visualize
    % lat is a vector of the latitudes you want to visualize
    % density_overlay is a boolean that will color an isosurface by the bulk density at that point
    % quiver_overlay is a boolean that will plot the motion vectors of the plume on the isosurface
    % dep_calc will is a boolean that will calculate the total ash being deposited on the ground over the course of the simulation

    % OUTPUTS
    % This code will output a text file with the following columns:
        % Vent speed (m/s)
        % Wind Speed (m/s)
        % stability mean
        % stability med
        % stability SD
        % Max plume height (km)
        % Neutral Buoyancy Height (km)
        % NBH err (km)
    % as well as a .mat file with the PDC directionality for each sumulation

    % EXAMPLES:
    % The following will visualize all tropical runs in the directory /Volumes/MATHAM3 with a vent diameter of 22.5 m. The isovalues are set to 1 and .001
    % ATHAM_viz_multi('ash3',[1,.001],false,'/Volumes/MATHAM3/',["22_5m"],["tropical"],false, false, false)

    % The following will visualize all tropical runs in the directory /Volumes/MATHAM3 with a vent diameter of 75m and 135m. The isovalue is set to .1
    % ATHAM_viz_multi('ash3',[.1],false,'/Volumes/MATHAM3/',["75m","135m"],["tropical"],false, false, false)

    % Make output directory for PDC directionality output
    mkdir(strcat(upper_dir,'direction_out'))
    % Loop through latitude vector
    for k = 1:length(lat)
        % Loop through vent diameter vector
        for j = 1:length(vent_diam)
            tic
            % Get full directories for netcdf file in vent diameter folder
            dir_info=dir(fullfile(strcat(upper_dir,lat(k),'/',vent_diam(j)),'**/atham_netCDF_MOV.nc'));
            dirs = extractfield(dir_info,'folder')';

            % Initialize output cell array labels
            output = {'Vent speed (m/s)','Wind Speed (m/s)','stability mean','stability med','stability SD','Max plume height (km)','Neutral Buoyancy Height (km)', 'NBH err (km)'};
            % Loop through directories found above
            for i = 1:length(dirs)
                fn = dirs{i};
                % Set stability plane calculation height based on vent diameter
                if vent_diam(j) == "15m" || vent_diam(j) == "22_5m"
                    plane_offset = 1;
                elseif vent_diam(j) == "75m" || vent_diam(j) == "135m" || vent_diam(j) == "315m" 
                    plane_offset = 0;
                end 
                
                % Call ATHAM_viz_ts to visualize the simulation and calculate stabiltiy + other stats
                [flux_ratio_mean, flux_ratio_med, flux_ratio_SD, max_plume_height, NBH, NBH_err,gif_str] = ATHAM_viz_ts(fn, tracer_name, isovalue, domain_flux, plane_offset, density_overlay, quiver_overlay, dep_calc);

                % Get output file name
                split_dir = regexp(fn,filesep,'split');
                output_fn = split_dir{end};

                % Get vent speed
                split_fn = regexp(output_fn,'_','split');
                vent_speed = split_dir{6};
                
                % Get wind speed
                if split_dir{3} == "22_5m"
                    wind_speed = split_fn{5};
                else
                    wind_speed = split_fn{4};
                end

                % Add stats to output cell array
                output{i+1,1} = vent_speed;
                output{i+1,2} = wind_speed;
                output{i+1,3} = flux_ratio_mean;
                output{i+1,4} = flux_ratio_med;
                output{i+1,5} = flux_ratio_SD;
                output{i+1,6} = max_plume_height;
                output{i+1,7} = NBH;
                output{i+1,8} = NBH_err;
                
                % Copy PDC directionality output to dirwextion_out directory made earlier
                copyfile(strcat(gif_str,'.mat'),strcat(upper_dir,'direction_out'))
            end
            close all
            % Go back to upper directory
            cd(upper_dir)
            
            % Write output cell array to text file
            writecell(output,strcat(lat(k),'_',vent_diam(j),'.txt'))

            disp(strcat('Done visualizing vent ',vent_diam(j),'in:'))
            toc
        end
    end
end

