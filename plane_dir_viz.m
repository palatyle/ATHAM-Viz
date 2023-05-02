close all
clear all

% Load in fine and coarse grid geometry files
load('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/fine_geo.mat')
load('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/coarse_geo.mat')

x_coarse = x(row_x);
y_coarse = y(row_y);

x_fine = x(row_x);
y_fine = y(row_y);

% Get list of all files in directory
cd '/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/direction_out'
files = dir();
files = files(~ismember({files.name},{'.','..','.DS_Store'}));

% Loop through all files and load in data
for i = 1:length(files)
    split = regexp(files(i).name,'_','split');
    for j = 1:length(split)
        letter_ind{j} = isletter(split{j});
    end


    split{2}(letter_ind{2}) = [];
    split{3}(letter_ind{3}) = [];
    split{4}(letter_ind{4}) = [];
        
    
    temp = load(files(i).name);
    temp.dir_matrix(isnan(temp.dir_matrix)) = 0;
    temp.dir_matrix = abs(temp.dir_matrix);
    planes(i).data = temp.dir_matrix;
    planes(i).atmosphere = string(split{1});
    if str2double(split{2}) == 127
        split{5}(letter_ind{5}) = [];
        planes(i).vent = 127.5;
        planes(i).v_velocity = str2double(split{4});
        planes(i).windspeed = str2double(split{5});
    else
        planes(i).vent = str2double(split{2});
        planes(i).v_velocity = str2double(split{3});
        planes(i).windspeed = str2double(split{4});
    end

    
end

vents = [20, 30, 75, 127.5, 303];

% Plot each direction plamne
for i = 1:length(vents)
    input_plane = planes([planes.vent]==vents(i));
    if vents(i) == 20
        plane_plot(input_plane,vents(i),x_fine,y_fine)
    else
        plane_plot(input_plane,vents(i),x_coarse,y_coarse)
    end

end


function plane_plot(plane_filt,vent_size,X,Y)
    if vent_size == 20
        window_size = 15;
    else
        window_size = 8;
    end

    % Windspeed planes
    plane_0 = plane_filt([plane_filt.windspeed] == 0);
    plane_5 = plane_filt([plane_filt.windspeed] == 5);
    plane_10 = plane_filt([plane_filt.windspeed] == 10);
    plane_15 = plane_filt([plane_filt.windspeed] == 15);
    plane_20 = plane_filt([plane_filt.windspeed] == 20);
    plane_25 = plane_filt([plane_filt.windspeed] == 25);
    plane_30 = plane_filt([plane_filt.windspeed] == 30);
    plane_35 = plane_filt([plane_filt.windspeed] == 35);
    plane_40 = plane_filt([plane_filt.windspeed] == 40);
    plane_45 = plane_filt([plane_filt.windspeed] == 45);
    plane_50 = plane_filt([plane_filt.windspeed] == 50);
    
    plane_0_sum = sum(cat(3,plane_0.data),3);
    plane_5_sum = sum(cat(3,plane_5.data),3);
    plane_10_sum = sum(cat(3,plane_10.data),3);
    plane_15_sum = sum(cat(3,plane_15.data),3);
    plane_20_sum = sum(cat(3,plane_20.data),3);
    plane_25_sum = sum(cat(3,plane_25.data),3);
    plane_30_sum = sum(cat(3,plane_30.data),3);
    plane_35_sum = sum(cat(3,plane_35.data),3);
    plane_40_sum = sum(cat(3,plane_40.data),3);
    plane_45_sum = sum(cat(3,plane_45.data),3);
    plane_50_sum = sum(cat(3,plane_50.data),3);
    
    plane_0_sum(plane_0_sum==0) = NaN;
    plane_5_sum(plane_5_sum==0) = NaN;
    plane_10_sum(plane_10_sum==0) = NaN;
    plane_15_sum(plane_15_sum==0) = NaN;
    plane_20_sum(plane_20_sum==0) = NaN;
    plane_25_sum(plane_25_sum==0) = NaN;
    plane_30_sum(plane_30_sum==0) = NaN;
    plane_35_sum(plane_35_sum==0) = NaN;
    plane_40_sum(plane_40_sum==0) = NaN;
    plane_45_sum(plane_45_sum==0) = NaN;
    plane_50_sum(plane_50_sum==0) = NaN;
    
    plane_cell = {plane_0_sum, plane_5_sum plane_10_sum plane_15_sum plane_20_sum plane_25_sum plane_30_sum plane_35_sum plane_40_sum plane_45_sum plane_50_sum};
    
%     color_min = min(min([plane_cell{:}]));
    color_max = max(max([plane_cell{:}]));
    
    plane_cell = cellfun(@(x) x/color_max,plane_cell,'un',0);
    windspeeds = 0:5:50;
    fig = figure;
    set(gcf,'Position', [454,86,890,714])

    tl = tiledlayout(3,4,'TileSpacing','compact');
    for i = 1:length(plane_cell)
%         sp{i} = subplot(3,4,i);
        h(i) = nexttile(tl);
        hp = pcolor(h(i),X,Y,plane_cell{i});
%         axis equal
        axis(h(i),'equal')
        xlim([(X(end)/2)-window_size (X(end)/2)+window_size])
        ylim([(X(end)/2)-window_size (X(end)/2)+window_size])
        if vent_size ~= 20
            set(gca,'xtick',[15:5:30],'ytick',[15:5:30])
        else
            set(gca,'xtick',[35:5:60],'ytick',[35:5:60])
        end

        set(gca,'XTickLabel',[],'YTickLabel',[])
        set(hp, 'EdgeColor', 'none');
        grid on

        [t,s] = title(strcat(string(windspeeds(i))," m/s"));
        t.FontName = 'Myriad Pro';
        t.FontWeight = 'Bold';
        t.FontSize = 16;
        set(gca, 'Layer', 'top')
        colormap(flipud(magma))
    end

    c=colorbar(h(end));
    c.Layout.Tile='east';
    
    ylabel(c,'Normalized Ash Mass');
    caxis([0,1])
    t_sg = sgtitle(strcat(num2str(vent_size),' m Vent'));
    t_sg.FontName = 'Myriad Pro';
    t_sg.FontWeight = 'Bold';
    t_sg.FontSize = 28;
    c.FontName = 'Myriad Pro';
    c.FontSize = 14;
end
