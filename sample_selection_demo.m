clear all
close all
clc

% Load configuration
config;

filtered_pointcloud = [];
%% Read sample pointcloud data
% cloud = read_point_cloud(raw_pc_path, true);    
cloud = read_point_cloud_cv(raw_pc_path, false); 

base_save = strcat('results/', exp_name);
if ~exist(base_save, 'dir')
    mkdir(base_save);
end


%% detect individual microgels:
[centers_x, centers_y, x_range, y_range] = detect_gels(cloud, 0);

%% Extract individual microgel and remove outlier 
for idx=1 : length(centers_x)
    gel1 = get_bounding_box(cloud, centers_x, centers_y, x_range, y_range, idx, 0, 0);
    filtered_point_cloud = removeOutliers(gel1.Location, 0.005);
    filtered_pointcloud{idx} = filtered_point_cloud;
    
end

noise_removed = [];
for idx=1 : length(filtered_pointcloud)
noise_removed = [noise_removed; filtered_pointcloud{idx}];
end
ptCloud_unnoisy = pointCloud(noise_removed);

%% Calculate correlation matrix
[sorted_corr_sum, sorted_corr_sum2, sortIdx]=calculate_corrolation_index(filtered_pointcloud, use_spherical_information);

%% Elbow strategy
% Find the change in correlation values
sorted_corr_sum_norm = sorted_corr_sum/max(sorted_corr_sum);
sorted_corr_sum_norm2 = sorted_corr_sum2/max(sorted_corr_sum2);
delta_corr = diff(sorted_corr_sum_norm);
% Calculate the second derivative to find the inflection point
second_derivative = diff(delta_corr);
% Find the index of the point with the maximum second derivative
[max_second_derivative, idx] = max(second_derivative);
% The index 'idx' corresponds to the elbow point
threshold_index = idx;


%% Visualization of the results
selected_color = [0, 255, 0];
dropped_color  = [255, 0, 0];
draw_color     = [0, 0, 255];


N = fix(length(filtered_pointcloud)*Percentage_select);
M = fix(length(filtered_pointcloud)*(1-Percentage_drop)); 

pcaggregate = [];
color_aggregated = [];
for idx=1 : length(filtered_pointcloud)
pcaggregate = [pcaggregate; filtered_pointcloud{sortIdx(idx)}];
if idx<=N
rgbList = repmat(selected_color, length(filtered_pointcloud{sortIdx(idx)}), 1);
elseif idx<M && idx>N
rgbList = repmat(draw_color, length(filtered_pointcloud{sortIdx(idx)}), 1);
else
rgbList = repmat(dropped_color, length(filtered_pointcloud{sortIdx(idx)}), 1);
end

color_aggregated = [color_aggregated; rgbList];

end
threshold_index = N;
ptCloud = pointCloud(pcaggregate, 'Color', color_aggregated/255);

figure(1)
pcshow(cloud)
set(gcf,'color','w');
set(gca,'color','w');
set(gca,'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
title('Raw microgel point cloud data from SRFM measurements');
savefig(strcat(base_save, '/Step1_raw_pointcloud.fig'));
saveas(gcf, strcat(base_save, '/Step1_raw_pointcloud.png'));

figure(2)
pcshow(ptCloud_unnoisy)
set(gcf,'color','w');
set(gca,'color','w');
set(gca,'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
title('Effect of removing noise');
savefig(strcat(base_save, '/Step2_denoized_raw_pointcloud.fig'));
saveas(gcf, strcat(base_save, '/Step2_denoized_raw_pointcloud.png'));

figure(3)
pcshow(ptCloud)
set(gcf,'color','w');
set(gca,'color','w');
set(gca,'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
title('Clustering result');
savefig(strcat(base_save, '/Step3_clustering_result.fig'));
saveas(gcf, strcat(base_save, '/Step3_clustering_result.png'));


%% Visualize the selected samples
selected_list = {};
for idx=1 : N
selected_list{idx} =  filtered_pointcloud{sortIdx(idx)};
end
totals = length(selected_list);
numRows = fix(sqrt(totals));
numCols = fix(sqrt(totals))+1;
figure(4);

% Loop through each image and display it in a subplot
for i = 1:min(numRows*numCols, numel(selected_list))
    subplot(numRows, numCols, i);
    x = selected_list{i};
    rgbList = repmat(selected_color, length(x), 1);
    x = x - mean(x);
    x = pointCloud(x, 'Color', rgbList/255);
    pcshow(x);
%     set(gcf,'color','w');
%     set(gca,'color','w');
%     set(gca,'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    zlim([-500 500])
    title(['Sample ', num2str(i)]);
    axis off;
end



%%
use_prompt = false
if use_prompt
% Prompt the user for input
userInput = inputdlg('Please input the indices of the microgels you wish to drop, separated by commas. If you do not want to drop any microgels, simply click OK', 'Input', [1 40]);

% Check if the user clicked cancel or entered an empty value
if isempty(userInput)
    disp('No input provided.');
    return;
end

% Split the input string into separate numbers using commas
inputNumbers = strsplit(userInput{1}, ',');

% Convert the input numbers to integers
integerList = str2double(inputNumbers);
if ~isnan(integerList)
selected_list(integerList) = [];
end
end

figure(5)
% Loop through each image and display it in a subplot
for i = 1:min(numRows*numCols, numel(selected_list))
    subplot(numRows, numCols, i);
    x = selected_list{i};
    rgbList = repmat(selected_color, length(x), 1);
    x = x - mean(x);
    x = pointCloud(x, 'Color', rgbList/255);
    pcshow(x);
%     set(gcf,'color','w');
%     set(gca,'color','w');
%     set(gca,'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    zlim([-500 500])
    title(['Sample ', num2str(i)]);
    axis off;
end
savefig(strcat(base_save, '/Step4_selected_samples.fig'));
saveas(gcf, strcat(base_save, '/Step4_selected_samples.png'));

save_path_selected = fullfile(base_save, 'plyfiles/');
if exist(save_path_selected, 'dir') == 0
    % Create the directory
    mkdir(save_path_selected);
end

for ids=1:numel(selected_list)
    % pccloud = pointCloud(selected_list{ids});
    % pcwrite(pccloud, strcat(save_path_selected, int2str(ids),'.ply'))
    xyz = selected_list{ids};
    filename = strcat(save_path_selected, int2str(ids),'.txt');
    dlmwrite(filename, xyz, ' ');
end

% Loop through each image and display it in a subplot
for i = 1:numel(selected_list)
    figure(20+i)
    PC = selected_list{i};
    ptCloud = pointCloud(PC);
    pcshow(ptCloud);
    savefig(strcat(save_path_selected, '/sample_',int2str(i),'_pointcloud.fig'));
    saveas(gcf, strcat(save_path_selected, '/sample_',int2str(i),'_pointcloud.png'));

end

%% Convexhull fitting
figure(6)
Agg_selected = [];
% Loop through each image and display it in a subplot
for i = 1:min(numRows*numCols, numel(selected_list))
    x_agg = selected_list{i};
    x_agg = x_agg - mean(x_agg);
    Agg_selected = [Agg_selected; x_agg];
end
Agg_selected = removeOutliers(Agg_selected, 0.005);
PC = Agg_selected;  % Replace with your point cloud data
ptCloud = pointCloud(PC);
pcshow(ptCloud);
% zlim([-500 500])
title('Representative sample');

% Save the figure as a .fig file
filename = strcat(save_path_selected, 'representative.txt');
dlmwrite(filename, Agg_selected, ' ');
savefig(strcat(base_save, '/Step5_representative_sample.fig'));
saveas(gcf, strcat(base_save, '/Step5_representative_sample.png'));
figure(7)

% Define the selected color with transparency for the point cloud
selected_color = [0, 255, 128] / 255;
% selected_color = [0.5, 0.5, 0.5]; % Gray color with transparency, modify as needed

% Create a pointCloud object with the specified color for all points
rgbaList = repmat(selected_color, length(PC), 1);
ptCloud = pointCloud(PC, 'Color', rgbaList);

% Calculate the centroid of the point cloud
centroid = mean(PC);

% Calculate the squared distance from each point to the centroid
distancesSquared = sum((PC - centroid).^2, 2);

% Calculate the radius of the sphere as the square root of the maximum squared distance
radius = sqrt(max(distancesSquared));

% Fit a 3D convex hull to the point cloud
k = convhull(PC);

% Extract the vertices of the convex hull
hullVertices = PC(unique(k(:)), :);

% Create a figure for visualization

hold on;

% Plot the convex hull vertices
scatter3(hullVertices(:, 1), hullVertices(:, 2), hullVertices(:, 3), 30, [1, 0.647, 0], 'filled');

% Plot the point cloud with transparency
pcshow(ptCloud);

% Define two points for the pointer line: the centroid and a point above the centroid
pointer_start = centroid;
pointer_end = centroid + [0, 0, radius*1.2]; % Extend the line upward

% Plot the pointer line with increased line width
plot3([pointer_start(1), pointer_end(1)], ...
      [pointer_start(2), pointer_end(2)], ...
      [pointer_start(3), pointer_end(3)], 'r', 'LineWidth', 4); % Increased line width

% Display the radius value on top of the pointer
text(pointer_end(1), pointer_end(2), pointer_end(3), ...
     sprintf('Radius: %.2f', radius), 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'center');

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Convexhull of the Microgel with Radiou');

% Set axis limits if needed
maxRange = max(radius, max(sqrt(sum(PC.^2, 2))));
xlim([-maxRange, maxRange]);
ylim([-maxRange, maxRange]);
zlim([-maxRange, maxRange]);

hold off;
set(gcf,'color','w');
set(gca,'color','w');
set(gca,'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
title('Representaive sample');
savefig(strcat(base_save, '/Step6_convex_fitting.fig'));
saveas(gcf, strcat(base_save, '/Step6_convex_fitting.png'));
%%
figure(8)
drop_sample = filtered_pointcloud{sortIdx(length(sortIdx)-1)};
Agg_selected = removeOutliers(drop_sample, 0.05);
Agg_selected = Agg_selected - mean(Agg_selected);
PC = Agg_selected;  % Replace with your point cloud data

% Define the selected color with transparency for the point cloud
selected_color = [0, 255, 128] / 255;
% selected_color = [0.5, 0.5, 0.5]; % Gray color with transparency, modify as needed

% Create a pointCloud object with the specified color for all points
rgbaList = repmat(selected_color, length(PC), 1);
ptCloud = pointCloud(PC, 'Color', rgbaList/255);

% Calculate the centroid of the point cloud
centroid = mean(PC);

% Calculate the squared distance from each point to the centroid
distancesSquared = sum((PC - centroid).^2, 2);

% Calculate the radius of the sphere as the square root of the maximum squared distance
radius = sqrt(max(distancesSquared));

% Fit a 3D convex hull to the point cloud
k = convhull(PC);

% Extract the vertices of the convex hull
hullVertices = PC(unique(k(:)), :);

% Create a figure for visualization

hold on;

% Plot the convex hull vertices
scatter3(hullVertices(:, 1), hullVertices(:, 2), hullVertices(:, 3), 30, [1, 0.647, 0], 'filled');

% Plot the point cloud with transparency
pcshow(ptCloud);

% Define two points for the pointer line: the centroid and a point above the centroid
pointer_start = centroid;
pointer_end = centroid + [0, 0, radius*1.2]; % Extend the line upward

% Plot the pointer line with increased line width
plot3([pointer_start(1), pointer_end(1)], ...
      [pointer_start(2), pointer_end(2)], ...
      [pointer_start(3), pointer_end(3)], 'r', 'LineWidth', 4); % Increased line width

% Display the radius value on top of the pointer
text(pointer_end(1), pointer_end(2), pointer_end(3), ...
     sprintf('Radius: %.2f', radius), 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'center');

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Convexhull of the Microgel with Radiou');

% Set axis limits if needed
maxRange = max(radius, max(sqrt(sum(PC.^2, 2))));
xlim([-maxRange, maxRange]);
ylim([-maxRange, maxRange]);
zlim([-maxRange, maxRange]);

hold off;
set(gcf,'color','w');
set(gca,'color','w');
set(gca,'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
title('dropped sample');

savefig(strcat(base_save, '/Step7_dropped_sample.fig'));
saveas(gcf, strcat(base_save, '/Step7_dropped_sample.png'));


text = ['The dropped sample has ' num2str(sorted_corr_sum_norm2(length(sortIdx)-1)) ' correlation with others'];


%%
% Create a pointCloud object with the specified color for all points
rgbaList = repmat(selected_color, length(PC), 1);
ptCloud = pointCloud(PC, 'Color', rgbaList/255);

% Fit a 3D convex hull to the point cloud
k = convhull(PC);

% Extract the vertices of the convex hull
hullVertices = PC(unique(k(:)), :);

% Calculate the centroid of the point cloud
centroid = mean(PC);

% Calculate the distance from each convex hull vertex to the centroid
distancesToCentroid = sqrt(sum((hullVertices - centroid).^2, 2));

% Calculate the radius of the sphere as the average distance from convex hull vertices to the centroid
radius = mean(distancesToCentroid);


hold on;

% Plot the convex hull vertices
scatter3(hullVertices(:, 1), hullVertices(:, 2), hullVertices(:, 3), 30, [1, 0.647, 0], 'filled');

% Plot the point cloud with transparency
pcshow(ptCloud);
% Define two points for the pointer line: the centroid and a point above the centroid
pointer_start = centroid;
pointer_end = centroid + [0, 0, radius*1.2]; % Extend the line upward


% Extract convex hull vertices in XY plane
hullXY = hullVertices(:, 1:2);

% Extract convex hull vertices in XZ plane
hullXZ = hullVertices(:, [1, 3]);

% Extract convex hull vertices in YZ plane
hullYZ = hullVertices(:, 2:3);

% Calculate diameters
diameterXY = mean( max(pdist2(hullXY, hullXY)));
diameterXZ = mean( max(pdist2(hullXZ, hullXZ)));
diameterYZ = mean( max(pdist2(hullYZ, hullYZ)));

% Calculate the spherical metric
sphericalMetric = min([diameterXY, diameterXZ, diameterYZ]) / max([diameterXY, diameterXZ, diameterYZ]);
formattedString = sprintf('%.2f', sphericalMetric);

text2 = ['The dropped sample has spherical score of ' num2str(sphericalMetric)];

fileID = fopen(strcat(base_save, '/report.txt'), 'w');
fprintf(fileID, '%s\n%s\n', text, text2);
fclose(fileID);

close all

function filteredPointCloud = removeOutliers(pointCloud, thresh)
    % Calculate the center of the point cloud
    center = mean(pointCloud);
    
    % Calculate the distances of each point to the center
    distances = sqrt(sum((pointCloud - center).^2, 2));
    
    % Sort the distances in descending order
    [~, sortedIndices] = sort(distances, 'descend');
    
    % Determine the number of points to remove (5% of the total points)
    numPointsToRemove = round(thresh * size(pointCloud, 1));
    
    % Remove the points with the highest distances
    filteredPointCloud = pointCloud;
    filteredPointCloud(sortedIndices(1:numPointsToRemove), :) = [];
end









