%% simple script used to plot how the camera angle changes per pixel in a specified column of the .ptx file

function camera_angle(input_file, col)

if nargin < 2
    col = 1;
end

fid = fopen(input_file);
width = str2num(fgets(fid));
height = str2num(fgets(fid));

for i = 1:8+(height*(col-1))
    fgets(fid);
end

y_coord = zeros(1, height);
skipped_points = 0;

for i=1:height
    line = strsplit(fgets(fid), ' ');
    x = str2double(line(1));
    y = str2double(line(2));
    z = str2double(line(3));
    
    r = sqrt(x^2 + y^2 + z^2);
    
    if r < 1e-6 
        skipped_points = skipped_points + 1;
        continue;
    end
    
    y_coord(i) = asin(z/r);
end

if skipped_points > 0
    disp(['Warning: found ', num2str(skipped_points), ' points that were too close to origin']);
end

x_coord = 1:height;
plot(x_coord, y_coord);

title([input_file, ' | Column ', num2str(col)]);