function [flux_ratio_mean, flux_ratio_med, flux_ratio_SD, max_plume_height, NBH, NBH_err] = ATHAM_viz_ts(fn, tracer_name, isovalue, domain_flux, small_grid, density_overlay, quiver_overlay, dep_calc)
%% before loop!
% Setup plot windows
% h = figure()
tic
cd(fn)
fn_data = 'atham_netCDF_MOV.nc';
fn_dep = 'atham_netCDF_PIC.nc';

ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);


% timestep = 45;
xlabel([ax1 ax2], 'X Distance [km]')
ylabel([ax1 ax2], 'Y Distance [km]')
zlabel([ax1 ax2], 'Z Distance [km]')
grid([ax1 ax2], 'on')
grid([ax1 ax2], 'minor')


% fn = '/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/tropical_75m_50ms_35ms/atham_netCDF_MOV.nc';

ash_names = find_num_trac(fn_data);

% Read in ash tracers in g/kg, multiply by scale factor and remove nans \

   
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

ash_iso = get_viz_trac(tracer_name);

% Vertical motion vector
if quiver_overlay
    u_vector = read_datafile(fn_data,'u');
    v_vector = read_datafile(fn_data,'v');
end
w_vector = read_datafile(fn_data,'w');


den_full = read_datafile(fn_data,'density');

den = den_full(:,:,:,1);

% pres = read_datafile(fn_data,'pnew');

den_full = remove_nans(den_full);



[x,y,z] = read_geo(fn_data);
[x_zoom_loc, y_zoom_loc] = read_zoom_loc(fn_data);
[xmg,ymg,zmg] = create_grid(x,y,z);
time_arr = read_datafile(fn_data,'time') * 60; %convert mins to seconds
time_num = find_timesteps(ash1);

if density_overlay || quiver_overlay
    % Find indices of isovalue at each timestep
    idxs = find(ash_iso(:) < isovalue+.01 & ash_iso(:) > isovalue-.01);
        iso_clr_max = max(max(max(max(den_full(idxs)))));
        iso_clr_min = min(min(min(min(den_full(idxs)))));


    color_min_max = [iso_clr_min,iso_clr_max];
end

% Get extent of plane 
[row_x,row_y] = get_plane_extent(x,y);
% Calcualte area of every grid cell
area_plane = area_calc(x,y,row_x,row_y);
% Find index of plane to calcualte stabiltiy at
plane_height = find_plane_height(den,x,y,z,xmg,ymg,small_grid);
lower_plane = round(plane_height/4);
% Get boolean array at plane height of volcano vs air
rad_dist_bool = get_rad_array(den,x,y,z,xmg,ymg,plane_height);

% Get lower boolean array 
rad_dist_bool_lower = den(:,:,lower_plane);
rad_dist_bool_lower(isnan(rad_dist_bool_lower)) = 0;
rad_dist_bool_lower(rad_dist_bool_lower ~= 0) = 1;

% Get horizontal extent of volcano base
volc_base = den(:,:,2);
volc_base(isnan(volc_base))= 0;
volc_base(volc_base ~= 0)= 1;

rad_dist_bool = rad_dist_bool + volc_base*2;

rad_dist_bool_lower = rad_dist_bool_lower + volc_base;
rad_dist_bool_lower(rad_dist_bool_lower==2) = 0;



rad_dist_bool = rad_dist_bool(row_x,row_y);
rad_dist_bool_lower = rad_dist_bool_lower(row_x,row_y);
% Density 
p_den = create_isosurf(ax1, xmg, ymg, zmg, remove_nans(den), .0005);
p_den = change_patch_props(p_den, true);
% Density zoomed
p_den_z = create_isosurf(ax2, xmg, ymg, zmg, remove_nans(den), .0005);
p_den_z = change_patch_props(p_den_z, true);

gif_str = get_gif_str(fn, isovalue);

disp(strcat('Currently visualizing ',' ', gif_str))
% figure(2);



%% In loop
for i = 1:time_num
    ash_viz = ash_iso(:,:,:,i);

    ash_threshold = ash_viz;
    ash_threshold(ash_threshold < isovalue) = 0;
    % Ash 1
    % ash1_ts = ash1(:,:,:,timestep);

    u_i = u_vector(:,:,:,i);
    v_i = v_vector(:,:,:,i);
    w_i = w_vector(:,:,:,i);

    if quiver_overlay
        quiv_idxs = find(ash_iso(:,:,:,i) < isovalue+.01 & ash_iso(:,:,:,i) > isovalue-.01);
%         q1 = quiver3(ax1, xmg(quiv_idxs),ymg(quiv_idxs),zmg(quiv_idxs),u_i(quiv_idxs),v_i(quiv_idxs),w_i(quiv_idxs),4,'LineWidth',2.5);
        hold(ax2,'on')
%         streamline(ax2, xmg,ymg,zmg,u_i,v_i,w_i,xmg())
        q = quiver3(ax2, xmg(quiv_idxs),ymg(quiv_idxs),zmg(quiv_idxs),u_i(quiv_idxs),v_i(quiv_idxs),w_i(quiv_idxs),2,'LineWidth',1);
        q.Color = [0 0.4470 0.7410];
    end
    if density_overlay
        p_ash = create_isosurf(ax1, xmg, ymg, zmg, ash_viz, isovalue, den_full(:,:,:,i));
        p_ash = change_patch_props(p_ash, false, color_min_max, ax1);
        hold(ax1, 'on');
        % Ash 1 zoomed 
        p_ash_z = create_isosurf(ax2, xmg, ymg, zmg, ash_viz, isovalue, den_full(:,:,:,i));
        p_ash_z = change_patch_props(p_ash_z, false, color_min_max, ax2);
        hold(ax2, 'on');
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
    set_lighting(ax1)
    set_cam(ax1, false, x, y, z, x_zoom_loc, y_zoom_loc)

    % Zoomed plot lighting and camera setup
    set_lighting(ax2)
    set_cam(ax2, true, x, y, z, x_zoom_loc, y_zoom_loc)

    % Find current plume height
    plume_height_ts(i) = find_plume_height(ash_threshold, z);
%     if i ~= 1
%         ash_profile = squeeze(ash_viz(144/2,144/2,:))';
%         figure; plot(ash_profile,z)
%     end
%     new_find_plume_height(ash_viz, w_vector(:,:,:,i), z)
    % Find max plume height through time
    if i ~= 1 && all(plume_height_ts(i) >= plume_height_ts(:))
        max_plume_height = plume_height_ts(i);
    elseif i == 1 
        max_plume_height = 0;
    end

    % Calclulate mass flux through each grid point
    
    grid_mass_flux = mass_flux(x,y,row_x, row_y, plane_height, den_full, area_plane, w_vector, i);

    grid_mass_flux_lower = mass_flux(x,y,row_x, row_y, lower_plane, den_full, area_plane, w_vector, i);

    if dep_calc
        mass_total(i) = grid_mass_sum(grid_mass_flux); 
    end

    
    % Calculate stability of plume
    [flux_ratio(i), stability(i)] = stability_calc(x,y,grid_mass_flux, grid_mass_flux_lower, ash1, ash2, ash3, ash4, row_x, row_y, plane_height, rad_dist_bool, rad_dist_bool_lower, i);
    

%     pause
    text_objs = draw_text(ax1, x, y, z, i, time_arr, flux_ratio(i), max_plume_height);
    
    % Set figure size and fontsize
    set_figure_props(18)
    
    shg
%     pause
    create_gif(i, strcat(gif_str,'.gif'))
    delete(findall(gcf,'Type','light'))
    
    patches = {p_ash p_ash_z q};
    for j = 1:length(patches)
        delete(patches{j})
    end
    for j = 1:length(text_objs)
        delete(text_objs{j})
    end
end

%% out of loop
if domain_flux == true
    domain_flux_calc(ash1, ash2, ash3, ash4, time_num, z)
end

[flux_ratio_mean, flux_ratio_med, flux_ratio_SD] = stat_calc(flux_ratio);

% Calculate neutral buoyancy height
[NBH, NBH_err] = NBH_calc(ash_threshold, x, y, z);

if dep_calc
    % Sum up final ash deposition
    ash_dep_mass = planar_density2mass(ash1_dep(row_x,row_y),area_plane) + planar_density2mass(ash2_dep(row_x,row_y),area_plane) + planar_density2mass(ash3_dep(row_x,row_y),area_plane) + planar_density2mass(ash4_dep(row_x,row_y),area_plane);
    
    % Sum up total amount of ash passing through stabiltiy calc plane
    
    mass_plane_final = sum(mass_total); 

    mass_plane_final/(mass_plane_final + ash_dep_mass)
end

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
%     data = ncread(fn,tracername);
end

function ash_names = find_num_trac(fn)
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

function [ash_ts] = get_viz_trac(tracer_name)
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
    iso_surf = isosurface(xmg,ymg,zmg,data,isovalue);
    p = patch(axnum,iso_surf);

    % calculate normals
    isonormals(xmg,ymg,zmg,data,p)

    if exist('den_ts','var')
        isocolors(xmg,ymg,zmg,den_ts,p);
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
    if ground == true 
        p.FaceColor = [uint8(158) uint8(111) uint8(57)]; %brown
        p.EdgeColor = 'none';
        p.FaceAlpha = 1;
    elseif ground == false
        if exist('color','var')
            p.FaceColor='interp';
            p.EdgeColor='none';
            if quiver_overlay
                p.FaceAlpha = 0.5;
            end
            colorbar(axnum)
            caxis(axnum,color)
            colormap(viridis())
        else
            p.FaceColor = [.3 .3 .3]; %grey
            p.EdgeColor = 'none';
            if quiver_overlay
                p.FaceAlpha = 0.5;
            else
                p.FaceAlpha = 1;
            end
            p.SpecularStrength = .1;
            p.DiffuseStrength = .1;  
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
    % Find where volcano stops and air begins
    k=0;
    all_nans_found = false;
    [upper_x_bool, lower_x_bool, upper_y_bool, lower_y_bool] = deal(true);
    while all_nans_found == false
        k = k+1;
        if isnan(density((length(x)/2)+k,length(y)/2,plane_height)) && upper_x_bool == true
            upper_x = (length(x)/2)+(k-1);
            upper_x_bool=false;
        end
        if isnan(density((length(x)/2)-k,length(y)/2,plane_height)) && lower_x_bool == true
            lower_x = (length(x)/2)-(k-1);
            lower_x_bool=false;
        end
        if isnan(density(length(x)/2,(length(y)/2)+k,plane_height)) && upper_y_bool == true
            upper_y = (length(y)/2)+(k-1);
            upper_y_bool=false;
        end
        if isnan(density(length(x)/2,(length(y)/2)-k,plane_height)) && lower_y_bool == true
            lower_y = (length(y)/2)-(k-1);
            lower_y_bool=false;
        end

        if upper_x_bool == false && lower_x_bool == false && upper_y_bool == false && lower_y_bool == false
            all_nans_found = true;
        end

    end

    rad_dist_bool = false(size(density(:,:,plane_height)));
    rad_dist_bool(lower_x:upper_x,lower_y:upper_y) = true;
end

function plane_height = find_plane_height(density,x,y,z,xmg,ymg,small_grid)
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
            if small_grid == true
                plane_height = i_func+1;
            else
                plane_height = i_func; % if not 30 m
            end
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
    
    mass_sum_ts = sum(grid_mass_flux(grid_mass_flux>0),'omitnan'); % mass sum in kgs
end

function final_dep_mass = planar_density2mass(dep,area_plane)
    % multiply planar density by area then multiply by 1000 to convert from g/m^2 to kg
    final_dep_mass_plane = (dep .* area_plane) / 1000; 
    final_dep_mass = sum(final_dep_mass_plane(:),'omitnan');
end

function [flux_ratio, stability] = stability_calc(x,y,grid_mass_flux,grid_mass_flux_lower, ash1, ash2, ash3, ash4, row_x, row_y, plane_height, rad_dist_bool, rad_dist_bool_lower, i)
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
    i (int): current timestep index
    Returns
    ----------
    stability (int): array containing mass flux through each cell in
    plane
    %}
    ash_weighted_flux = grid_mass_flux .* (ash1(row_x,row_y,plane_height,i) ...
        + ash2(row_x,row_y,plane_height,i) ...
        + ash3(row_x,row_y,plane_height,i) ...
        + ash4(row_x,row_y,plane_height,i));
%     ash_weighted_flux = grid_mass_flux .* ash4(row_x,row_y,plane_height,i);
    
%     pos_flux = sum(ash_weighted_flux(ash_weighted_flux>0),'omitnan');
%     neg_flux = sum(ash_weighted_flux(ash_weighted_flux<0),'omitnan');
%     tot_flux = sum(sum(abs(ash_weighted_flux),'omitnan'));
%     flux_ratio = pos_flux/tot_flux;
%     flux_ratio = pos_flux/-neg_flux;


    inner = sum(ash_weighted_flux(rad_dist_bool==1 & ash_weighted_flux > 0),'omitnan');
%     inner = sum(ash_weighted_flux(ash_weighted_flux > 0),'omitnan');
    outer = sum(ash_weighted_flux(rad_dist_bool==0 & ash_weighted_flux < 0),'omitnan');
     
    inner_GMF = sum(grid_mass_flux(rad_dist_bool==1 & grid_mass_flux > 0),'omitnan');
    outer_GMF = sum(grid_mass_flux(rad_dist_bool==0 & grid_mass_flux < 0),'omitnan');

    outer_lower_GMF = sum(grid_mass_flux_lower(rad_dist_bool_lower==1 & grid_mass_flux_lower < 0),'omitnan');
    
    "upper, lower plane comparison"
    flux_ratio = inner_GMF/(abs(outer_lower_GMF)+inner_GMF)
    "upper plane comparison"
    inner_GMF/(abs(outer_GMF)+inner_GMF)
%     flux_ratio = abs(inner)/(abs(outer)+abs(inner));
    
    if flux_ratio > 0.80
        stability = 1; %stable
    elseif flux_ratio <= 0.80 && flux_ratio >= .2
        stability = 2; %partially  unstable
    elseif flux_ratio < .2 
        stability = 3; %fully unstable
    else
        stability = NaN;
    end

end

function p_height = new_find_plume_height(ash_viz, w_vector, z)
    for j_func = length(z):-1:1
        up_v = w_vector(:,:,j_func);
        [min_v,loc_min_v] = min(up_v(:));
        [max_v,loc_max_v] = max(up_v(:));
    end
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
    
    for j_func = length(z):-1:1
        if sum(sum(ash_threshold(:,:,j_func))) ~= 0
            plume_height = z(j_func);
            break
        end
    end
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
    
    % sum 
    for i_func = 1:time_num
        ash1_row_ts(:,i_func) = row_sum_ts(ash1(:,:,:,i_func));
        ash2_row_ts(:,i_func) = row_sum_ts(ash2(:,:,:,i_func));
        ash3_row_ts(:,i_func) = row_sum_ts(ash3(:,:,:,i_func));
        ash4_row_ts(:,i_func) = row_sum_ts(ash4(:,:,:,i_func));
    end
    total_ash1_ts = sum(ash1_row_ts,2);
    total_ash2_ts = sum(ash2_row_ts,2);
    total_ash3_ts = sum(ash3_row_ts,2);
    total_ash4_ts = sum(ash4_row_ts,2);
    
    ash1_domain = sum_remaining_ash(ash1, z);
    ash2_domain = sum_remaining_ash(ash2, z);
    ash3_domain = sum_remaining_ash(ash3, z);
    ash4_domain = sum_remaining_ash(ash4, z);
    
    total_domain_flux_ash1 = total_ash1_ts + ash1_domain';
    total_domain_flux_ash2 = total_ash2_ts + ash2_domain';
    total_domain_flux_ash3 = total_ash3_ts + ash3_domain';
    total_domain_flux_ash4 = total_ash4_ts + ash4_domain';
    writematrix(['Height (km)','37 um ash conc (g/kg)','150 um ash conc (g/kg)','900 um ash conc (g/kg)'],'mass_fluxes.txt','Delimiter',',')
    writematrix([z,total_domain_flux_ash1,total_domain_flux_ash2,total_domain_flux_ash3],'mass_fluxes.txt','WriteMode','append')
end

function last_step_domain_ash_var = sum_remaining_ash(ash_var, z)
    
    for i_func = 1:length(z)
        last_step_domain_ash_var(i_func) = sum(sum(ash_var(:,:,i_func)));
    end
    
end

function  row_sum = row_sum_ts(ash_var)
    
    
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
    split_dir = regexp(fn,filesep,'split');
    gif_str = split_dir{end};
    if isempty(gif_str)
        gif_str = split_dir{end-1};
    end
    gif_str = strcat(gif_str,'_iso_',num2str(isovalue));
end

function [flux_ratio_mean, flux_ratio_med, flux_ratio_SD] = stat_calc(flux_ratio)
    flux_ratio_mean = mean(flux_ratio(10:end-1),'omitnan');
    flux_ratio_med = median(flux_ratio(10:end-1),'omitnan');
    flux_ratio_SD = std(flux_ratio(10:end-1),'omitnan');
end

function [NBH, NBH_err] = NBH_calc(ash_threshold,x,y,z)
    for i_func = length(y):-1:1 %Go backwards through XZ plane (for some reason x and y have to be flipped here, god knows why. I'm so sorry to whoever has to use this code after me
        if max(max(squeeze(ash_threshold(i_func,:,:)))) ~= 0
            break
        end
    end
    
    NBH_vert_profile = find(ash_threshold(i_func-5,length(x)/2,:)); %Find nonzero indices in threshold profile
    if ~isempty(NBH_vert_profile) 
        NBH = z(ceil(median(NBH_vert_profile))); %median of profile (rounded up) gives index in z of NBH
        NBH_err = abs(z(min(NBH_vert_profile)) - z(max(NBH_vert_profile))); %Full span of nonzero numbers gives error
    else
        NBH = 0 ;
        NBH_err = 0;
    end
             
end

function create_gif(i,gif_str)
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