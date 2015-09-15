%% simple script used to plot how the camera angle changes per pixel in a specified column of the .ptx file

function camera_angle(input_file, col)

if nargin < 2
    col = 1;
end

fid = fopen(input_file);
width = str2num(fgets(fid));
height = str2num(fgets(fid));

disp(['Found ', num2str(height), ' pixels per column']);

for i = 1:8+(height*(col-1))
    fgets(fid);
end

x_coord = 1:height;
y_coord = zeros(1, height);
skipped_points = 0;
num_points = 0;

x_sum = 0;
y_sum = 0;

for i=1:height
    line = strsplit(fgets(fid), ' ');
    x = str2double(line(1));
    y = str2double(line(2));
    z = str2double(line(3));
    
    r = sqrt(x^2 + y^2 + z^2);
    
    if r < 1e-6
        skipped_points = skipped_points + 1;
        y_coord(i) = NaN;
        continue;
    end
    
    num_points = num_points + 1;
    y_coord(i) = asin(z/r);
    
    x_sum = x_sum + x_coord(i);
    y_sum = y_sum + y_coord(i);
end

if skipped_points > 0
    disp(['Warning: skipped ', num2str(skipped_points), ' points that were too close to origin']);
end

disp(['Angle of first pixel is ', num2str(y_coord(1)), ' rad, ', num2str(y_coord(1)*180/pi), ' deg']);
disp(['Angle of last pixel is ', num2str(y_coord(end)), ' rad, ', num2str(y_coord(end)*180/pi), ' deg']);

plot(x_coord, y_coord);
title([input_file, ' | Column ', num2str(col)]);

if num_points > 0
    x_avg = x_sum / num_points;
    y_avg = y_sum / num_points;
    
    x_coord(isnan(y_coord)) = [];
    y_coord(isnan(y_coord)) = [];
    
    slope = sum((x_coord - x_avg) .* (y_coord - y_avg)) / sum((x_coord - x_avg).^2);
    intercept = y_avg - slope*x_avg;
    
    x_fit = x_coord;
    y_fit = slope*x_fit + intercept;
    
    figure;
    plot(x_fit, y_fit);
    title([input_file, ' | Best fit line']);
    
    disp(['Best fit line (least squares): y = ', num2str(slope), 'x + ', num2str(intercept)]);
end