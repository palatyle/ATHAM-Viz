function [flux_ratio_mean, flux_ratio_med, flux_ratio_SD, max_plume_height, NBH, NBH_err,gif_str] = ATHAM_viz_ts(fn, tracer_name, isovalue, domain_flux, plane_offset, density_overlay, quiver_overlay, dep_calc)
% ATHAM_viz_ts.m  Visualize ATHAM output as a 3D animation. 

% INPUTS
% fn = directory of ATHAM output file (string)
% tracer_name = name of tracer to visualize (string)
% isovalue = isovalue of tracer to visualize (double)
% domain_flux is a boolean that determines whether or not you want to visualize the flux through the domain edges of the simulation
% plane_offset = index offset of plane from vent
% density_overlay is a boolean that will color an isosurface by the bulk density at that point
% quiver_overlay is a boolean that will plot the motion vectors of the plume on the isosurface
% dep_calc will is a boolean that will calculate the total ash being deposited on the ground over the course of the simulation

% OUTPUTS
% flux_ratio_mean = mean stability flux ratio over the course of the simulation
% flux_ratio_med = median stability flux ratio over the course of the simulation
% flux_ratio_SD = standard deviation of stability flux ratio over the course of the simulation
% max_plume_height = maximum plume height over the course of the simulation (km)
% NBH = Neutral Buoyancy Height of the plume over the course of the simulation (km) NOTE: currently unreliable)_
% NBH_err = error in NBH (km) NOTE: currently unreliable
% gif_str = string of the name of the animated gif created by this function

% EXAMPLES:
% The following will visualize the ATHAM netcdf file located at /Volumes/MATHAM3/tropical/30m/50ms/tropical_30m_50ms_0ms.
% It will visualzie the ash3 traceer at two isosurface (.1 and .001 g ash/kg air) and offsets the stabiltiy plane by 1 index value.
% [flux_ratio_mean, flux_ratio_med, flux_ratio_SD, max_plume_height, NBH, NBH_err] = ATHAM_viz_ts("/Volumes/MATHAM3/tropical/30m/50ms/tropical_30m_50ms_0ms", 'ash3', [.1,.001], false, 1, false, false, false) 

%% before loop!

tic
cd(fn)
% Assuming all ATHAM output files are named in the following format:
fn_data = 'atham_netCDF_MOV.nc';
% fn_data = 'MATHAM_spring_netCDF_MOV.nc';
fn_dep = 'atham_netCDF_PIC.nc';

% Setting up the figure
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);

xlabel([ax1 ax2], 'X Distance [km]')
ylabel([ax1 ax2], 'Y Distance [km]')
zlabel([ax1 ax2], 'Z Distance [km]')
grid([ax1 ax2], 'on')
grid([ax1 ax2], 'minor')

% Get tracer names from netcdf file
ash_names = find_num_trac(fn_data);

% Read in ash tracers in g/kg, multiply by scale factor and remove nans 
if sum(matches(ash_names, 'ash1')) == 1
    ash1 = read_datafile(fn_data,'ash1')*1000;
    ash1 = remove_nans(ash1);
    if dep_calc
        ash1_dep = read_datafile(fn_dep,'d_ash1')*1000;
        ash1_dep = remove_nans(ash1_dep);
    end
end
    
if sum(matches(ash_names, 'ash2')) == 1
    ash2 = read_datafile(fn_data,'ash2')*1000;
    ash2 = remove_nans(ash2);
    if dep_calc
        ash2_dep = read_datafile(fn_dep,'d_ash2')*1000;
        ash2_dep = remove_nans(ash2_dep);
    end
else
    ash2 = zeros(size(ash1));
end

if sum(matches(ash_names, 'ash3')) == 1
    ash3 = read_datafile(fn_data,'ash3')*1000;
    ash3 = remove_nans(ash3);
    if dep_calc
        ash3_dep = read_datafile(fn_dep,'d_ash3')*1000;
        ash3_dep = remove_nans(ash3_dep);
    end
else
    ash3 = zeros(size(ash1));
end

if sum(matches(ash_names, 'ash4')) == 1
    ash4 = read_datafile(fn_data,'ash4')*1000;
    ash4 = remove_nans(ash4);
    if dep_calc
        ash4_dep = read_datafile(fn_dep,'d_ash4')*1000;
        ash4_dep = remove_nans(ash4_dep);
    end
else
    ash4 = zeros(size(ash1));
end

% Make new variable for visualization
ash_iso = get_viz_trac(tracer_name);

% If guiver overlay is turned on, read in u and v vectors
if quiver_overlay
    u_vector = read_datafile(fn_data,'u'); % m/s
    v_vector = read_datafile(fn_data,'v'); % m/s
end

% Read in w vector for use in stability calculation
w_vector = read_datafile(fn_data,'w'); % m/s

% Read in density for use in stability calculation
den_full = read_datafile(fn_data,'density'); % kg/m^3

% Read in pressure and temperature for use in stability calculation
pnew = read_datafile(fn_data,'pnew'); % Pa
tempnew = read_datafile(fn_data,'tempnew'); % K

% Calculate air density for entire 4D array
p_air = calc_air_den(pnew,tempnew); % kg/m^3

% Calculate ash mass per unit volume [kg ash/m^3 air]
ash1_air_den = ash_air_den_calc(ash1,p_air);
ash2_air_den = ash_air_den_calc(ash2,p_air);
ash3_air_den = ash_air_den_calc(ash3,p_air);
ash4_air_den = ash_air_den_calc(ash4,p_air);

% Grab first timestep of density array. To be used in visualization of ground surface
den = den_full(:,:,:,1);

% Remove nans from entire density array
den_full = remove_nans(den_full);

% Read in x,y,z coordinates
[x,y,z] = read_geo(fn_data); % m

% Read in center of domain
[x_zoom_loc, y_zoom_loc] = read_zoom_loc(fn_data);

% Create mesh grid out of x,y,z coordinates
[xmg,ymg,zmg] = create_grid(x,y,z);

% Read in time array and conver to seconds
time_arr = read_datafile(fn_data,'time') * 60; %convert mins to seconds

% Find timestep value
time_num = find_timesteps(ash1);

% Find max and min of density and quiver overlays 
if density_overlay || quiver_overlay
    idxs = find(ash_iso(:) < isovalue+.01 & ash_iso(:) > isovalue-.01);
        iso_clr_max = max(max(max(max(den_full(idxs)))));
        iso_clr_min = min(min(min(min(den_full(idxs)))));


    color_min_max = [iso_clr_min,iso_clr_max];
end

% Get extent of stability plane 
[row_x,row_y] = get_plane_extent(x,y);

% Calcualte area of every grid cell
area_plane = area_calc(x,y,row_x,row_y); % m^2

% Find height index of plane to calcualte stabiltiy at
plane_height = find_plane_height(den,x,y,z,xmg,ymg,plane_offset);

% Find height index of plane to track PDC directionality at
lower_plane = round(plane_height*(7/8));

% Get boolean array of all points within radius of crater rim
rad_dist_bool = get_rad_array(den,x,y,z,xmg,ymg,plane_height);

% Constrain boolean array to only points within plane extent
rad_dist_bool = rad_dist_bool(row_x,row_y);

% Create isosurface for volcano
p_den = create_isosurf(ax1, xmg, ymg, zmg, remove_nans(den), .0005);
change_patch_props(p_den, true);

% Zoomed isosurface
p_den_z = create_isosurf(ax2, xmg, ymg, zmg, remove_nans(den), .0005);
change_patch_props(p_den_z, true);

% creeate str to be used in gif filename
gif_str = get_gif_str(fn, isovalue);

% initialize matrix to store direction of PDC
dir_matrix = zeros(size(area_plane));

disp(strcat('Currently visualizing',{' '}, gif_str))

%% In loop
% Loop through timesteps
for i = 1:time_num
    % Create new 3D array of ash concentration for each timestep
    ash_viz = ash_iso(:,:,:,i);

    % Create another, but just for ash threshold calculatiomns
    ash_threshold = ash_viz;
    ash_threshold(ash_threshold < isovalue(1)) = 0;

    % If quiver overlay is turned on, 3D array of u and v vectors at timestep
    if quiver_overlay
        u_i = u_vector(:,:,:,i);
        v_i = v_vector(:,:,:,i);
    end

    % Grab 3D array of w vector at timestep
    w_i = w_vector(:,:,:,i);

    % If quiver overlay is turned on, find quiver indexes amd plot quiver on those indices
    if quiver_overlay
        quiv_idxs = find(ash_iso(:,:,:,i) < isovalue+.01 & ash_iso(:,:,:,i) > isovalue-.01);
        hold(ax2,'on')
        q = quiver3(ax2, xmg(quiv_idxs),ymg(quiv_idxs),zmg(quiv_idxs),u_i(quiv_idxs),v_i(quiv_idxs),w_i(quiv_idxs),2,'LineWidth',1);
        q.Color = [0 0.4470 0.7410];
    end

    % If density overlay is turned on, color isosurface by density at timestep
    if density_overlay
        p_ash = create_isosurf(ax1, xmg, ymg, zmg, ash_viz, isovalue, den_full(:,:,:,i));
        p_ash = change_patch_props(p_ash, false, color_min_max, ax1);
        hold(ax1, 'on');
        % Ash 1 zoomed 
        p_ash_z = create_isosurf(ax2, xmg, ymg, zmg, ash_viz, isovalue, den_full(:,:,:,i));
        p_ash_z = change_patch_props(p_ash_z, false, color_min_max, ax2);
        hold(ax2, 'on');
    % Otherwise, color isosurface standard grey,. 
    else
        p_ash = create_isosurf(ax1, xmg, ymg, zmg, ash_viz, isovalue);
        p_ash = change_patch_props(p_ash, false);
        hold(ax1, 'on');
        % Ash 1 zoomed 
        p_ash_z = create_isosurf(ax2, xmg, ymg, zmg, ash_viz, isovalue);
        p_ash_z = change_patch_props(p_ash_z, false);
        hold(ax2, 'on');
    end


    % Main plot lighting and camera setup
    set_cam(ax1, false, x, y, z, x_zoom_loc, y_zoom_loc)
    set_lighting(ax1)
    % Zoomed plot lighting and camera setup
    set_lighting(ax2)
    set_cam(ax2, true, x, y, z, x_zoom_loc, y_zoom_loc)

    % Find current plume height
    plume_height_ts(i) = find_plume_height(ash_threshold, z);

    % Check to see if current plume height is greater than previous max
    if i ~= 1 && all(plume_height_ts(i) >= plume_height_ts(:))
        max_plume_height = plume_height_ts(i);
    elseif i == 1 
        max_plume_height = 0;
    end

    % Calclulate mass flux through each grid point at vent plane
    grid_mass_flux = mass_flux(x,y,row_x, row_y, plane_height, den_full, area_plane, w_vector, i);
    grid_mass_flux_ash = ash_mass_flux(ash1_air_den,area_plane,w_vector,row_x,row_y,plane_height,i) ...
        + ash_mass_flux(ash2_air_den,area_plane,w_vector,row_x,row_y,plane_height,i) ...
        + ash_mass_flux(ash3_air_den,area_plane,w_vector,row_x,row_y,plane_height,i) ...
        + ash_mass_flux(ash4_air_den,area_plane,w_vector,row_x,row_y,plane_height,i);

    % Calclulate mass flux through each grid point at lower plane
    grid_mass_flux_lower = mass_flux(x,y,row_x, row_y, lower_plane, den_full, area_plane, w_vector, i);
    grid_mass_flux_lower_ash = ash_mass_flux(ash1_air_den,area_plane,w_vector,row_x,row_y,lower_plane,i) ...
        + ash_mass_flux(ash2_air_den,area_plane,w_vector,row_x,row_y,lower_plane,i) ...
        + ash_mass_flux(ash3_air_den,area_plane,w_vector,row_x,row_y,lower_plane,i) ...
        + ash_mass_flux(ash4_air_den,area_plane,w_vector,row_x,row_y,lower_plane,i);
    
    % Ignore positive values and sum into dir_matrix for each timestep
    grid_mass_flux_lower_ash(grid_mass_flux_lower_ash>0) = 0; 
    dir_matrix = dir_matrix + grid_mass_flux_lower_ash;

    if dep_calc
        mass_total(i) = grid_mass_sum(grid_mass_flux); 
    end

    % Calculate stability of plume
    [flux_ratio(i), stability(i)] = stability_calc(grid_mass_flux_ash, ash1, ash2, ash3, ash4, row_x, row_y, plane_height, rad_dist_bool, i);
    
    % Draw text objects 
    text_objs = draw_text(ax1, x, y, z, i, time_arr, flux_ratio(i), max_plume_height);
    
    % Set figure size and fontsize
    set_figure_props(18)
    
    shg

    create_gif(i, strcat(gif_str,'.gif'))

    % Delete all patches and lights so that they don't overlap on next timestep
    delete(findall(gcf,'Type','light'))
    if quiver_overlay
        patches = {p_ash p_ash_z q};
    else
        patches = {p_ash p_ash_z};
    end
    
    for j = 1:length(patches)
        delete(patches{j})
    end
    for j = 1:length(text_objs)
        delete(text_objs{j})
    end
end

%% out of loop

% Calculate domain flux if true
if domain_flux == true
    domain_flux_calc(ash1, ash2, ash3, ash4, time_num, z)
end

% Calculate statistics of flux ratio
[flux_ratio_mean, flux_ratio_med, flux_ratio_SD] = stat_calc(flux_ratio);

% Calculate neutral buoyancy height (unreliable)
[NBH, NBH_err] = NBH_calc(ash_threshold, x, y, z);

% Calculate deposited ash if set to true
if dep_calc
    % Sum up final ash deposition
    ash_dep_mass = planar_density2mass(ash1_dep(row_x,row_y),area_plane) + planar_density2mass(ash2_dep(row_x,row_y),area_plane) + planar_density2mass(ash3_dep(row_x,row_y),area_plane) + planar_density2mass(ash4_dep(row_x,row_y),area_plane);
    
    % Sum up total amount of ash passing through stabiltiy calc plane
    
    mass_plane_final = sum(mass_total); 

    mass_plane_final/(mass_plane_final + ash_dep_mass)
end

% Save PDC directionality matrix 
save(strcat(gif_str,'.mat'),'dir_matrix')

disp('Done visualizing in:')
toc

function data = read_datafile(fn, tracername)
    %{
    Returns array of the specified tracername 

    Parameters
    ----------
    fn (str): full filename of netcdf data file
    tracername (str): tracer name of interest found in netcdf file

    Returns
    ----------
    data (array): data array likely containing ash concentration or some other
    parameter
    %}
    ncid=netcdf.open(fn,'NC_NOWRITE');
    varid=netcdf.inqVarID(ncid,tracername);
    if strcmp(tracername,'x') || strcmp(tracername,'y') || strcmp(tracername,'z')
        data = netcdf.getVar(ncid,varid,'single');
    else
        % Read in as double and set no data value to NaNs
        data = netcdf.getVar(ncid,varid,'double');
        data(data==-9.899999922296324e+32) = NaN;

    end
    netcdf.close(ncid)
end

function ash_names = find_num_trac(fn)
    %{
    Returns array of names of ash tracers

    Parameters
    ----------
    fn (str): full filename of netcdf data file

    Returns
    ----------
    ash_names (array): data array of ash tracer names
    %}
    file_info = ncinfo(fn);
    for i_func = 1:length(file_info.Variables)
        var_names{i_func} = file_info.Variables(i_func).Name;
    end
    
    var_names = string(var_names);
    ash_idx = matches(var_names,["ash1","ash2","ash3","ash4"]);
    ash_names = var_names(ash_idx);
    
end

function [x,y,z] = read_geo(fn)
    %{
    Returns x, y, and z vectors of grid 

    Parameters
    ----------
    fn (str): full filename of netcdf data file

    Returns
    ----------
    x,y,z (1D vectors): x, y, and z values containing grid information 
    %}

    x = read_datafile(fn,'x')/1000;
    y = read_datafile(fn,'y')/1000;
    z = read_datafile(fn,'z')/1000;
end

function [x_zoom_loc, y_zoom_loc] = read_zoom_loc(fn)
    %{
    Returns the zoom location in x and y (the (rough) center of the grid)

    Parameters
    ----------
    fn (str): full filename of netcdf data file

    Returns
    ----------
    x_zoom_loc, y_zoom_loc (1D vectors): zoom locations 
    %}
    x_zoom_loc = ncreadatt(fn,'/','x_zoom_location')/1000;
    y_zoom_loc = ncreadatt(fn,'/','y_zoom_location')/1000;
end

function data = remove_nans(data)
    %{
    Returns data array with nans repalced with zeros 

    Parameters
    ----------
    data (array): data array (w/nans)

    Returns
    ----------
    data (array): data array (no nans)
    %}

    data(isnan(data)) = 0;
end

function [xmg,ymg,zmg] = create_grid(x,y,z)
    %{
    Returns meshgrid varaibles in x, y, and z. For use in plotting

    Parameters
    ----------
    x, y, z (vectors): Input x, y, and z variables that will define grid
    space

    Returns
    ----------
    xmg, ymg, zmg (arrays): meshgrids for each dimension
    %}

    [xmg,ymg,zmg] = meshgrid(x,y,z);
end

function timestep_num = find_timesteps(data)
    %{
    Returns the amount of timesteps present in data file

    Parameters
    ----------
    data (array): data array 

    Returns
    ----------
    timestep_num (int): integer of number of timesteps present in data file
    %}

    data_dim = size(data);
    timestep_num = data_dim(4);
end

function ash_ts = get_viz_trac(tracer_name)

    if strcmp(tracer_name, 'ash1')
        ash_ts = ash1(:,:,:,:);
    elseif strcmp(tracer_name, 'ash2')
        ash_ts = ash2(:,:,:,:);
    elseif strcmp(tracer_name, 'ash3')
        ash_ts = ash3(:,:,:,:);
    elseif strcmp(tracer_name, 'ash4')
        ash_ts = ash4(:,:,:,:);
    end 
end

function p = create_isosurf(axnum, xmg, ymg, zmg, data, isovalue, den_ts)
    %{
    Returns a patch object of isosurface. Plots data onto axes axnum at
    the corresponding isovalue

    Parameters
    ----------
    axnum (obj): axes object
    xmg, ymg, zmg (arrays): meshgrids for each dimension
    data (array): data array
    isovalue (num): Concetration value to create isosurface at

    Returns
    ----------
    p (obj): Patch object of isosurface
    %}     
    for i_func = 1:length(isovalue)
        iso_surf(i_func) = isosurface(xmg,ymg,zmg,data,isovalue(i_func));
        p(i_func) = patch(axnum,iso_surf(i_func));
    
        % calculate normals
        isonormals(xmg,ymg,zmg,data,p(i_func))
    
        if exist('den_ts','var')
            isocolors(xmg,ymg,zmg,den_ts,p(i_func));
        end
    end
end

function p = change_patch_props(p, ground, color, axnum)
    %{
    Returns a patch object of isosurface. Changes properties of patch
    dependign on whether object is of ground or ash

    Parameters
    ----------
    p (obj): Patch object of isosurface
    ground (bool): Boolean controlling whether patch is of ground or ash

    Returns
    ----------
    p (obj): Patch object of isosurface
    %}
    for i_func = 1:length(p)
        if ground == true 
            p(i_func).FaceColor = [uint8(158) uint8(111) uint8(57)]; %brown
            p(i_func).EdgeColor = 'none';
            p(i_func).FaceAlpha = 1;
        elseif ground == false
            if exist('color','var')
                p(i_func).FaceColor='interp';
                p(i_func).EdgeColor='none';
                if quiver_overlay
                    p(i_func).FaceAlpha = 0.5;
                end
                colorbar(axnum)
                caxis(axnum,color)
                colormap(viridis())
            else
                p(i_func).FaceColor = [.3 .3 .3]; %grey
                p(i_func).EdgeColor = 'none';
                if quiver_overlay
                    p(i_func).FaceAlpha = 0.5;
                else
                    p(i_func).FaceAlpha = 1;
                end
                p(i_func).SpecularStrength = .1;
                p(i_func).DiffuseStrength = .1;  
            if i_func > 1
                p(i_func).FaceAlpha = 0.3;
            end
            end
        end
    end
end

function set_lighting(axnum)
    %{
    Sets lighting for isosurface on current aces

    Parameters
    ----------
    axnum (obj): axes object
    ground (bool): Boolean controlling whether patch is of ground or ash

    Returns
    ----------
    %}
    camlight(axnum) 
    lighting(axnum, 'gouraud')
    material(axnum, 'dull')
end

function set_cam(axnum, zoom, x, y, z, x_zoom_loc, y_zoom_loc)
    %{
    Sets camera for zoomed and full domain cases 

    Parameters
    ----------
    axnum (obj): axes object
    zoom (bool): Boolean controlling whether camera is on full domain or
    zoomed in to vent
    x,y,z (1D vectors): x, y, and z values containing grid information
    x_zoom_loc, y_zoom_loc (1D vectors): zoom locations 

    Returns
    ----------
    %}
    if zoom == true
        view(axnum,73,25)
        
        xlim(axnum,[x_zoom_loc-5 x_zoom_loc+5])
        ylim(axnum,[y_zoom_loc-5 y_zoom_loc+5])
        zlim(axnum,[0 4.3])
    elseif zoom == false
        view(axnum,90,0)
        % view(-62.907905123533659,14.341768859374485)
        
        xlim(axnum,[min(x) max(x)])
        ylim(axnum,[min(y) max(y)])
        zlim(axnum,[min(z) max(z)])  
    end
    axis(axnum, 'square')
end

function set_figure_props(fontsize)
    %{
    Sets figure properties for final output gif, including figure size and
    text size

    Parameters
    ----------
    fontsize (int): fontsize 

    Returns
    ----------
    %}
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    set(findall(gcf,'type','text'),'FontSize',fontsize)
    
    ax_num = findall(gcf, 'type', 'axes');
    for i_func = 1:length(ax_num)
        set(ax_num(i_func), 'Fontsize', fontsize)
    end
end

function rad_dist_bool = get_rad_array(density,x,y,z,xmg,ymg,plane_height)
    %{
    Finds the radial distance from the vent to the crater rim edge. Returns a 2D array

    Parameters
    ----------
    density (arr): 3D density array
    x (arr): x values
    y (arr): y values
    z (arr): z values
    xmg, ymg (arr): meshgrids for x and y
    plane_height (int): height of plane to find crater rim edge

    Returns
    ----------
    rad_dist_bool (arr): 2D array of radial distance from vent to crater rim edge
    %}

    k=0;
    % initialize variables
    all_nans_found = false;
    [upper_x_bool, lower_x_bool, upper_y_bool, lower_y_bool] = deal(true);
    rad_dist_bool = false(size(density(:,:,plane_height)));

    % find x and y values stepping out from center until nans are found
    while all_nans_found == false
        k = k+1;
        dist_away_from_cent = x((length(x)/2)+k) - x((length(x)/2));

        if isnan(density((length(x)/2)+k,length(y)/2,plane_height)) && upper_x_bool == true
            upper_x = (length(x)/2)+(k-1);
            upper_x_bool=false;
        elseif dist_away_from_cent >= .3 && upper_x_bool == true
            upper_x = (length(x)/2)+(k+1);
            upper_x_bool=false;
        end
        if isnan(density((length(x)/2)-k,length(y)/2,plane_height)) && lower_x_bool == true
            lower_x = (length(x)/2)-(k-1);
            lower_x_bool=false;
        elseif dist_away_from_cent >= .3 && lower_x_bool == true
            lower_x = (length(x)/2)-(k+1);
            lower_x_bool=false;
        end
        if isnan(density(length(x)/2,(length(y)/2)+k,plane_height)) && upper_y_bool == true
            upper_y = (length(y)/2)+(k-1);
            upper_y_bool=false;
        elseif dist_away_from_cent >= .3 && upper_y_bool == true
            upper_y = (length(y)/2)+(k+1);
            upper_y_bool=false;
        end
        if isnan(density(length(x)/2,(length(y)/2)-k,plane_height)) && lower_y_bool == true
            lower_y = (length(y)/2)-(k-1);
            lower_y_bool=false;
        elseif dist_away_from_cent >= .3 && lower_y_bool == true
            lower_y = (length(y)/2)-(k+1);
            lower_y_bool=false;            
        end

        if upper_x_bool == false && lower_x_bool == false && upper_y_bool == false && lower_y_bool == false
            all_nans_found = true;
            rad_dist_bool(lower_x:upper_x,lower_y:upper_y) = true;
        end
    end
        % If nan is found to early, something is wrong (in our case, this
        % only happens in the 303 m case. Instead of looking for nans,
        % we're going to look for lack of nans in a very similar way to
        % above. 
        if k == 1 && all_nans_found == true
            more_nans_needed = false;
            [upper_x_bool, lower_x_bool, upper_y_bool, lower_y_bool] = deal(true);
            while more_nans_needed == false
                k = k+1;
                if ~isnan(density((length(x)/2)+k,length(y)/2,plane_height)) && upper_x_bool == true
                    upper_x = (length(x)/2)+(k-1);
                    upper_x_bool=false;
                end
                if ~isnan(density((length(x)/2)-k,length(y)/2,plane_height)) && lower_x_bool == true
                    lower_x = (length(x)/2)-(k-1);
                    lower_x_bool=false;
                end
                if ~isnan(density(length(x)/2,(length(y)/2)+k,plane_height)) && upper_y_bool == true
                    upper_y = (length(y)/2)+(k-1);
                    upper_y_bool=false;
                end
                if ~isnan(density(length(x)/2,(length(y)/2)-k,plane_height)) && lower_y_bool == true
                    lower_y = (length(y)/2)-(k-1);
                    lower_y_bool=false;          
                end
                if upper_x_bool == false && lower_x_bool == false && upper_y_bool == false && lower_y_bool == false
                    more_nans_needed = true;
                end
            end
            rad_dist_bool(lower_x+5:upper_x-5,lower_y+5:upper_y-5) = true;
        end
end

function plane_height = find_plane_height(density,x,y,z,xmg,ymg,plane_offset)
    %{
    Finds height point just above vent to place stabiltiy calculation vent.
    Counts up from bottom of domain untill reaching a height where there
    are no nans(i.e. where the volcano ground stops)

    Parameters
    ----------
    density (array): 3D density array (at one timestep
    z (1D vector): z values containing vertical grid information

    Returns
    ----------
    plane_height (int): index of stabiltiy plane. 
    %}
    for i_func = 1:length(z)
%         if sum(sum(isnan(density(:,:,i_func)))) == 0
        if squeeze(isnan(density(length(x)/2,length(y)/2,i_func))) == 0
            plane_height = i_func + plane_offset;
            break
        end
    end
end

function [row_x,row_y] = get_plane_extent(x,y)
    %{
    Gets indices of a plane 1 index away from edges

    Parameters
    ----------
    x,y (1D vector): x and y values containing horizontal grid information

    Returns
    ----------
    row_x,row_y (1D vectors): indices of plane in both x and y dimensions
    %}
    row_x = 2:length(x)-1;
    row_y = 2:length(y)-1;
end

function area_plane = area_calc(x,y,row_x,row_y)
    %{
    Calculates area of each grid cell in the plane defiend by row_x and
    row_y

    Parameters
    ----------    
    x,y (1D vector): x and y values containing horizontal grid information
    row_x,row_y (1D vectors): indices of plane in both x and y dimensions
    
    Returns
    ----------
    area_plane (Array): array containing area of each chell in plane
    %}
    for i_func = 1:length(row_x)
        for j_func = 1:length(row_y)
            area_plane(i_func,j_func) = ((x(row_x(i_func)-1)/2 - x(row_x(i_func)+1)/2)*1000) * ((y(row_y(j_func)-1)/2 - y(row_y(j_func)+1)/2)*1000);
        end
    end

end


function p_air = calc_air_den(pnew,tempnew)
    %{
    Calculates air density based on pressure and temperature using the ideal gas law. 

    Parameters
    ----------
    pnew (4D array): pressure values
    tempnew (4D array): temperature values

    Returns
    ----------
    p_air (4D array): air density values
    %}
    R = 287.05; % J/kg*K (https://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html) for air
    p_air = pnew./(R.*tempnew); % kg/m^3
end

function ash_air_den = ash_air_den_calc(ash_trac,air_den)
    %{
    Convert ash concentration to ash mass per volume of air.

    Parameters
    ----------
    ash_trac (4D array): Ash concentration values [g ash/kg air]
    air_den (4D array): Air density values [kg/m^3]

    Returns
    ----------
    ash_air_den (4D array): ash concentration per volume of air [kg ash/m^3 air]
    %}

    % factor of 1000 to convert from g to kg
    ash_air_den = (ash_trac./1000) .* air_den; 
end

function grid_mass_flux_ash = ash_mass_flux(ash_air_den,area_plane,w,row_x,row_y,plane_height,i)
    %{
    Calculates ash mass flux based on ash mass per unit volume, area of plane, and upward velocity vector

    Parameters
    ----------
    ash_air_den (4D array): ash concentration per volume of air [kg ash/m^3 air]
    area_plane (4D array): array containing area of each cell in plane [m^2]
    w (4D array): upward velocity vector [m/s]
    row_x,row_y (1D vectors): indices of plane in both x and y dimensions
    plane_height (int): index representing height above the domain where vent exists. 
    i (int): index of timestep
    
    Returns
    ----------
    grid_mass_flux_ash (Array): array containing ash mass flux of each chell in plane [kg/s]
    %}
    grid_mass_flux_ash = ash_air_den(row_x,row_y,plane_height,i) .* area_plane .* w(row_x,row_y,plane_height,i);
end

function grid_mass_flux = mass_flux(x,y,row_x, row_y,plane_height, density, area_plane, w, i)
    %{
    Calculates mass flux through each grid cell

    Parameters
    ----------
    row_x,row_y (1D vectors): indices of plane in both x and y dimensions
    plane_height (int): index representing height above the domain where
    vent exists. 
    density (4D array): density array through time
    area_plane (Array): array containing area of each chell in plane
    w (4D array): Vertical motion vector through time
    i (int): current timestep index
    Returns
    ----------
    grid_mass_flux (Array): array containing mass flux through each cell in
    plane
    %}
    
    density_plane = density(row_x,row_y,plane_height,i);
    w_plane = w(row_x,row_y,plane_height,i);
    
    grid_mass_flux = density_plane .* w_plane .* area_plane;
%     sum(sum(grid_mass_flux))
%     h1 = figure;
%     pcolor(grid_mass_flux)
%     close(h1)
end
function mass_sum_ts = grid_mass_sum(grid_mass_flux)
    %{
    Sums up positive maass flux in grid. 

    Parameters
    ----------
    % grid_mass_flux (Array): array containing mass flux through each cell in vent plane
    Returns
    ----------
    % mass_sum_ts (float): sum of mass flux through vent plane [kg]
    %}
    mass_sum_ts = sum(grid_mass_flux(grid_mass_flux>0),'omitnan'); % mass sum in kgs
end

function final_dep_mass = planar_density2mass(dep,area_plane)
    %{
    calculates deposited ash mass based on planar density and area of plane 

    Parameters
    ----------
    % dep (Array): array containing planar density of ash in each cell of plane [g/m^2]
    % area_plane (Array): array containing area of each chell in plane [m^2]
    Returns
    ----------
    % final_dep_mass (float): sum of mass flux through vent plane [kg]
    %}

    % multiply planar density by area then multiply by 1000 to convert from g/m^2 to kg
    final_dep_mass_plane = (dep .* area_plane) / 1000; 
    final_dep_mass = sum(final_dep_mass_plane(:),'omitnan');
end

function flux_ratio = stability_calc(grid_mass_flux, ash1, ash2, ash3, ash4, row_x, row_y, plane_height, rad_dist_bool, i)
    %{
    Calculates stability based on mass flux through the plane along with
    ash concentrations

    Parameters
    ----------
    grid_mass_flux (Array): array containing mass flux through each cell in
    plane
    ash1 (4D array): ash1 array through time
    ash2 (4D array): ash2 array through time
    ash3 (4D array): ash3 array through time
    ash4 (4D array): ash4 array through time
    row_x,row_y (1D vectors): indices of plane in both x and y dimensions
    plane_height (int): index representing height above the domain where
    vent exists. 
    rad_dist_bool (Array): array containing boolean values for within crater rim vs outside crater rim
    i (int): current timestep index
    Returns
    ----------
    flux_ratio (float): Calculated flux ratio for current timestep
    %}
    
    % Inner mass flux
    inner_GMF = sum(grid_mass_flux(rad_dist_bool==1 & grid_mass_flux > 0),'omitnan');

    % Outer mass flux
    outer_GMF = sum(grid_mass_flux(rad_dist_bool==0 & grid_mass_flux < 0),'omitnan');

    % Ratio of inner to outer mass flux
    flux_ratio = inner_GMF/(abs(outer_GMF)+inner_GMF);

end

function plume_height = find_plume_height(ash_threshold,z)
    %{
    Finds the max plume height based on the isosurface being visualized 

    Parameters
    ----------
    ash_threshold (3D array): Current timestep ash array with value below isovalue set to 0 
    z (1D vector): z values containing vertical grid information
    
    Returns
    ----------
    max_plume_height (float): maximum plume height
    %}
    
    % Loop downwards until ash threshold is not 0. This is the max plume height
    for j_func = length(z):-1:1
        if sum(sum(ash_threshold(:,:,j_func))) ~= 0
            plume_height = z(j_func);
            break
        end
    end

    % If ash threshold is 0 for all values, set plume height to 0
    if max(max(max(ash_threshold))) == 0
        plume_height = 0;
    end
end

function text_objs = draw_text(ax1, x, y, z, i, time_arr, flux_ratio, max_plume_height)
    %{
    Places text objects in 3D space on the plot. Text objects include
    current timestep, stability ratio, and max plume height

    Parameters
    ----------    
    ax1 (obj): Axes object to place the text onto
    x,y,z (1D vector): x, y, and z values containing  grid information
    i (int): current index
    time_arr (1D vector): vector containing timestep information
    flux_ratio (num): stabiltiy flux ratio
    max_plume_height (num): Maximum plume height at timestep
    
    Returns
    ----------
    text_objs (obj cell array): cell array containign each text object 
    %}
    
    % current timestep
    text_objs{1} = text(ax1, double(max(x)), double(min(y)) + 1.5 , double(max(z)) - 1.5, strcat(num2str(round(time_arr(i))),' s'));
    % Stability ratio
    text_objs{2} = text(ax1, double(max(x)), double(min(y)) + 1.5 , double(max(z)) - 4.5, strcat('stability ratio: ',' ', num2str(flux_ratio)));
    % Maximum plume height
    text_objs{3} = text(ax1, double(max(x)), double(max(y)) - 32, double(max(z)) - 1.5, strcat('Plume height: ', num2str(max_plume_height),' km'));
end

function domain_flux_calc(ash1, ash2, ash3, ash4, time_num, z)
    %{
    Calculates the flux of ash leaving the model domain at every timestep.
    Also sums up any remaining ash in domain at last timestep. sums data
    into 1D vertical profile for input into the LMD GCM

    Parameters
    ----------    
    ash1 (4D array): ash1 array through time
    ash2 (4D array): ash2 array through time
    ash3 (4D array): ash3 array through time
    ash4 (4D array): ash4 array through time
    time_num (num): total number of timesteps
    z (1D vector): z values containing vertical grid information
    
    Returns
    ----------
    %}
    
    for i_func = 1:time_num-1
        ash1_row_ts(:,i_func) = row_sum_ts(ash1(:,:,:,i_func));
        ash2_row_ts(:,i_func) = row_sum_ts(ash2(:,:,:,i_func));
        ash3_row_ts(:,i_func) = row_sum_ts(ash3(:,:,:,i_func));
        ash4_row_ts(:,i_func) = row_sum_ts(ash4(:,:,:,i_func));
    end
    total_ash1_ts = sum(ash1_row_ts,2);
    total_ash2_ts = sum(ash2_row_ts,2);
    total_ash3_ts = sum(ash3_row_ts,2);
    total_ash4_ts = sum(ash4_row_ts,2);
    
    ash1_domain = sum_remaining_ash(ash1(:,:,:,end), z);
    ash2_domain = sum_remaining_ash(ash2(:,:,:,end), z);
    ash3_domain = sum_remaining_ash(ash3(:,:,:,end), z);
    ash4_domain = sum_remaining_ash(ash4(:,:,:,end), z);
    
    total_domain_flux_ash1 = total_ash1_ts + ash1_domain';
    total_domain_flux_ash2 = total_ash2_ts + ash2_domain';
    total_domain_flux_ash3 = total_ash3_ts + ash3_domain';
    total_domain_flux_ash4 = total_ash4_ts + ash4_domain';
    writematrix(['Height (km)','7_8 um ash conc (g/kg)','15_6 um ash conc (g/kg)','125 um ash conc (g/kg)','1 mm ash conc (g/kg)'],'mass_fluxes.txt','Delimiter',',')
    writematrix([z,total_domain_flux_ash1,total_domain_flux_ash2,total_domain_flux_ash3,total_domain_flux_ash4],'mass_fluxes.txt','WriteMode','append')
end

function last_step_domain_ash_var = sum_remaining_ash(ash_var, z)
    %{
    Sums up the remaining ash in the domain at the last timestep

    Parameters
    ----------    
    ash_var (4D array): ash array through time
    z (1D vector): z values containing vertical grid information
    
    Returns
    ----------
    last_step_domain_ash_var (1D vector): 1D vector containing the sum of ash in the domain at the last timestep
    %}
    for i_func = 1:length(z)
        last_step_domain_ash_var(i_func) = sum(sum(ash_var(:,:,i_func)));
    end
    
end

function row_sum = row_sum_ts(ash_var)
     %{
    Sums up the ash in the domain at each timestep row wise throughout the domain

    Parameters
    ----------    
    ash_var (4D array): ash array through time

    Returns
    ----------
    row_sum (1D vector): 1D vector containing the sum of ash in the domain at each timestep
    %}
    domain_exit1 = squeeze(ash_var(:,end,:)); % side
    domain_exit2 = squeeze(ash_var(:,1,:)); % side
    domain_exit3 = squeeze(ash_var(end,:,:)); % side
    domain_exit4 = squeeze(ash_var(1,:,:)); % side
    domain_exit5 = squeeze(ash_var(:,:,end)); % top

    % sum row wise 
    dom1_sum = sum(domain_exit1);
    dom2_sum = sum(domain_exit2);
    dom3_sum = sum(domain_exit3);
    dom4_sum = sum(domain_exit4);
    dom5_sum = sum(sum(domain_exit5)); %except here, just sum the whole shebang 
    
    
    row_sum = dom1_sum + dom2_sum + dom3_sum + dom4_sum;
    row_sum(end) = row_sum(end) + dom5_sum; % add in top domain to last row
    

end

function gif_str = get_gif_str(fn, isovalue)
    %{
    Creates a string to be used in the gif filename

    Parameters
    ----------    
    fn (str): filename of the netcdf ATHAM output file
    isovalue (1D vector): vector containing the isovalue(s) used in the visualization

    Returns
    ----------
    gif_str (str): string to be used in the gif filename
    %}

    % Split fn into distinct parts
    split_dir = regexp(fn,filesep,'split');
    gif_str = split_dir{end};
    if isempty(gif_str)
        gif_str = split_dir{end-1};
    end
    temp_str=[];
    for i_func = 1:length(isovalue)
        temp_str = strcat(temp_str,'_',num2str(isovalue(i_func)));
    end
    gif_str = strcat(gif_str,'_iso',temp_str);
end

function [flux_ratio_mean, flux_ratio_med, flux_ratio_SD] = stat_calc(flux_ratio)
    %{
    Calculates the mean, median, and standard deviation of the flux ratio vector
    
    Parameters
    ----------    
    flux_ratio (1D vector): vector containing the flux ratio values for every timestep

    Returns
    ----------
    flux_ratio_mean (float): mean of the flux ratio vector
    flux_ratio_med (float): median of the flux ratio vector
    flux_ratio_SD (float): standard deviation of the flux ratio vector
    %}
    flux_ratio_mean = mean(flux_ratio(1:end-1),'omitnan');
    flux_ratio_med = median(flux_ratio(1:end-1),'omitnan');
    flux_ratio_SD = std(flux_ratio(1:end-1),'omitnan');
end

function [NBH, NBH_err] = NBH_calc(ash_threshold,x,y,z)
    %{
    Calculates the Neutral Buoyancy Height (NBH) and NBH error from the ash threshold array. Currently unreliable
    
    Parameters
    ----------
    ash_threshold (3D array): ash threshold array
    x,y,z (1D vectors): x,y,z values containing grid information
    
    Returns
    ----------
    NBH (float): Nvalue
    NBH_err (float): NBH error
    %}

    % Go backwards through XZ plane (for some reason x and y have to be flipped here, god knows why. 
    % I'm so sorry to whoever has to use this code after me)
    for i_func = length(y):-1:1 
        if max(max(squeeze(ash_threshold(i_func,:,:)))) ~= 0
            break
        end
    end
    try
        NBH_vert_profile = find(ash_threshold(i_func-5,length(x)/2,:)); %Find nonzero indices in threshold profile
        if ~isempty(NBH_vert_profile) 
            NBH = z(ceil(median(NBH_vert_profile))); %median of profile (rounded up) gives index in z of NBH
            NBH_err = abs(z(min(NBH_vert_profile)) - z(max(NBH_vert_profile))); %Full span of nonzero numbers gives error
        else
            NBH = 0 ;
            NBH_err = 0;
        end
    catch 
        warning("NBH calc failed. Assigning NBH value to 0")
        NBH = 0;
        NBH_err = 0;
    end
             
end

function create_gif(i,gif_str)
    %{
    Creates a gif
    
    Parameters
    ----------
    i (int): timestep
    gif_str (str): string to be used in the gif filename
    
    Returns
    ----------
    None

    %}

    % Get frame from current figure and add it to the gif
    F = getframe(gcf);
    im = frame2im(F);
    [imind,cm] = rgb2ind(im,256);
    if i == 1 
      imwrite(imind,cm,gif_str,'gif', 'Loopcount',inf,'DelayTime', .1); 
    else
      imwrite(imind,cm,gif_str,'gif','WriteMode','append','DelayTime', .1); 
    end
end
end