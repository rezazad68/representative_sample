function cloud = read_point_cloud_cv(filename, apply_flip)

data = readmatrix(filename);

% Select the columns second to fourth
data = data(:, 2:4);

if apply_flip
data(:,3) = 1-data(:,3); %flip z dimension
end
cloud = pointCloud(data);
pcshow(cloud); 
end