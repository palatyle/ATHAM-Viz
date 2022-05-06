function ATHAM_viz_multi(tracer_name,isovalue,domain_flux,upper_dir,vent_diam,lat,small_grid)

% tracer_name = 'ash1';
% isovalue = .001;
% domain_flux = false;
% small_grid = false;
% upper_dir = '/Volumes/(M)ATHAM/';
% vent_diam = '75m';
% lat = 'polar';
dir_info=dir(fullfile(strcat(upper_dir,lat,'/',vent_diam),'**/atham_netCDF_MOV.nc'));
dirs = extractfield(dir_info,'folder')';

output = {'Vent speed (m/s)','Wind Speed (m/s)','stability mean','stability med','stability SD','Max plume height (km)','Neutral Buoyancy Height (km)', 'NBH err (km)'};
for i = 1:length(dirs)
    fn = dirs{i};
    [flux_ratio_mean, flux_ratio_med, flux_ratio_SD, max_plume_height, NBH, NBH_err] = ATHAM_viz_ts(fn, tracer_name, isovalue, domain_flux, small_grid);
    split_dir = regexp(fn,filesep,'split');
    output_fn = split_dir{end};
    
    split_fn = regexp(output_fn,'_','split');
    vent_speed = split_fn{3};
    wind_speed = split_fn{4};
    
    output{i+1,1} = vent_speed;
    output{i+1,2} = wind_speed;
    output{i+1,3} = flux_ratio_mean;
    output{i+1,4} = flux_ratio_med;
    output{i+1,5} = flux_ratio_SD;
    output{i+1,6} = max_plume_height;
    output{i+1,7} = NBH;
    output{i+1,8} = NBH_err;
end
cd(upper_dir)
writecell(output,strcat(lat,'_',vent_diam,'.txt'))
end

