close all

% Loop through each image and display it in a subplot
for i = 1:numel(selected_list)
    PC = selected_list{i};
    ptCloud = pointCloud(PC);
    pcshow(ptCloud);
    savefig(strcat(base_save, '/sample_',int2str(i),'_pointcloud.fig'));
    saveas(gcf, strcat(base_save, '/sample_',int2str(i),'_pointcloud.png'));

end


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


