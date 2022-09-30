close all
clear all
cd '/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/direction_out'

load('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/coarse_geo.mat')
x_coarse = x(row_x);
y_coarse = y(row_y);


load('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/fine_geo.mat')
x_fine = x(row_x);
y_fine = y(row_y);

files = dir();

files = files(~ismember({files.name},{'.','..','.DS_Store'}));


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
% 
% for i = 0:5:50
%     sum
% end

vents = [20, 30, 75, 127.5, 303];

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
        window_size = 7;
    end

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
    
    color_min = min(min([plane_cell{:}]));
    color_max = max(max([plane_cell{:}]));
    
    cellfun(@(x) x/color_max,plane_cell,'un',0);
    windspeeds = 0:5:50;
    fig = figure;
    for i = 1:length(plane_cell)
        sp{i} = subplot(3,4,i);
        hp = pcolor(X,Y,plane_cell{i});
        grid on
        xlim([(X(end)/2)-window_size (X(end)/2)+window_size])
        ylim([(Y(end)/2)-window_size (Y(end)/2)+window_size])
        set(hp, 'EdgeColor', 'none');
        axis square
        grid on
        title(strcat(string(windspeeds(i))," m/s"))
%         caxis(sp{i},[color_min,color_max]) %
    end
    h = axes(fig, 'visible','off');
    colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
%     caxis(h,[color_min,color_max])
    sgtitle(strcat(num2str(vent_size),' m Vent'))

end
