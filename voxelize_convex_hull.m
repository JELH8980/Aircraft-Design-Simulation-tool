function elements = voxelize_convex_hull(points, nx, ny, nz)
% VOXELIZE_CONVEX_HULL - Partitions a convex hull into cubic voxels and returns voxel properties.
%
% This function takes a set of 3D points defining a convex hull, discretizes the enclosed volume into a grid of cubic voxels, 
% and computes the centroids and properties of the voxels that lie within the hull. It uses a specified number of divisions 
% along each axis (X, Y, Z) to determine voxel size and placement, returning only those centroids inside the convex hull along 
% with voxel dimensions and volume.
%
% INPUTS:
%   points     - Nx3 matrix of convex hull vertices, where each row is [x, y, z] coordinates.
%   nx         - Number of voxels along the X-axis (integer).
%   ny         - Number of voxels along the Y-axis (integer).
%   nz         - Number of voxels along the Z-axis (integer).
%
% OUTPUTS:
%   elements   - A struct containing voxel properties:
%                - elements.centroids: Mx3 matrix of voxel centroids inside the convex hull, where each row is [x, y, z].
%                - elements.volume: Scalar value representing the volume of each voxel (dx*dy*dz).
%                - elements.dx: Voxel size along the X-axis.
%                - elements.dy: Voxel size along the Y-axis.
%                - elements.dz: Voxel size along the Z-axis.
%
% FUNCTIONALITY:
% - Computes the convex hull of the input points using the `convhull` function, returning vertex indices (K).
% - Determines the bounding box of the hull by finding the minimum and maximum coordinates along each axis.
% - Calculates voxel sizes (dx, dy, dz) by dividing the bounding box dimensions by the specified number of voxels (nx, ny, nz).
% - Generates a grid of candidate voxel centroids using `ndgrid`, placing centroids at the center of each voxel.
% - Filters the candidate centroids to keep only those inside the convex hull using the `inhull` function.
% - Stores the valid centroids, voxel volume, and individual voxel dimensions in the output struct.
%
% NOTES:
% - The function assumes the `inhull` function is available to test point inclusion within the convex hull.
% - Voxel centroids are offset by half the voxel size (e.g., dx/2) to position them at the center of each voxel.
% - The output `elements.centroids` will have M rows, where M <= nx*ny*nz, depending on how many voxels lie fully within the hull.
% - Voxel sizes (dx, dy, dz) are uniform within each dimension but may differ between dimensions if the bounding box is not cubic.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

% Compute convex hull
[K, ~] = convhull(points);

% Find bounding box of the hull
min_vals = min(points, [], 1);
max_vals = max(points, [], 1);

% Compute voxel size in each dimension
dx = (max_vals(1) - min_vals(1)) / nx;
dy = (max_vals(2) - min_vals(2)) / ny;
dz = (max_vals(3) - min_vals(3)) / nz;

% Generate voxel centroid coordinates
[x, y, z] = ndgrid(linspace(min_vals(1) + dx/2, max_vals(1) - dx/2, nx), ...
                    linspace(min_vals(2) + dy/2, max_vals(2) - dy/2, ny), ...
                    linspace(min_vals(3) + dz/2, max_vals(3) - dz/2, nz));
centroid_candidates = [x(:), y(:), z(:)];

% Check which centroids are inside the convex hull
in_hull = inhull(centroid_candidates, points, K);

% Store valid centroids
centroids = centroid_candidates(in_hull, :);

elements.centroids = centroids;
elements.volume = dx*dy*dz;
elements.dx= dx;
elements.dy = dy;
elements.dz = dz;


end

