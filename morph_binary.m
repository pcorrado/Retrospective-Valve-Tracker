function [interpMask] = morph_binary(mask1,mask2,n)
% MORPH_BINARY takes two 2D binary images and "morphs" one binary mask
% into the other (with n intervening steps) using interpolation. This 
% results in a 3D binary output composed of mask1, n interpolated slices,
% and mask2. This function may be useful when attempting to create a
% volume-of-interest on a 3D image - the user can create masks only on
% certain slices (for instance, every 10 slices) and then interpolate the
% intervening ones to create a full 3D binary volume.
% 
% Inputs: mask1 and mask2 - 2D logical mask representing "slices" of a 3D
%                           binary image stack. mask1 and mask2 must have
%                           the same array size.
%                       n - The number of intervening slices between mask1
%                           and mask2
% 
% Outputs:     interpMask - 3D output mask composed of mask1, mask2, and n
%                           intervening slices interpolated from mask1 and
%                           mask2. If mask1 and mask2 are of size [A,B],
%                           interpMask will be of size [A,B,n+2]
% 
% Example:
% 
%       % Create sample 2D masks 'A' (a square) and 'B' (a diamond)
%       A = false(100,100); A(25:75,25:75) = true;
%       B = B = poly2mask([10,50,90,50],[50,10,50,90],100,100);
%   
%       % Morph A into B using 18 intervening slices, creating 3D array 'C'
%       C = morph_binary(A,B,18); % C has size [100,100,20]
% 
% Copyright: Michael D. Newton 2017
%            Orthopaedic Research Laboratories
%            Beaumont Health
% 
% Contact: michael.newton@beaumont.org

%% Create bwperim of each mask
% Define the boundaries of the binary masks

mask1P = bwperim(mask1);
mask2P = bwperim(mask2);

%% Create distance maps of each mask
% Compute the distance of each pixel to the perimeter

mask1D = bwdist(mask1P);
mask2D = bwdist(mask2P);

%% Invert values of distance maps inside the perimeter
% Make all of the values inside the perimeter negative so that there are
% negative values inside, positive values outside, and zeros at the
% perimeter.

mask1I = mask1 & ~mask1P;
mask2I = mask2 & ~mask2P;

mask1D(mask1I == 1) = -mask1D(mask1I == 1);
mask2D(mask2I == 1) = -mask2D(mask2I == 1);

%% Build stack
% Create an image stack from the distance maps for interpolation

sz = size(mask1);
stackD = zeros([sz,2]);
stackD(:,:,1  ) = mask1D;
stackD(:,:,end) = mask2D;

%% Interpolate values of intervening slices
% Interpolate between the distance maps to generate maps for intervening
% slices

[x ,y ,z ] = size(stackD);
[xG,yG,zG] = meshgrid(1:y,1:x,1:z);
[xQ,yQ,zQ] = meshgrid(1:y,1:x,1:1/(n+1):2);

VQ = interp3(xG,yG,zG,stackD,xQ,yQ,zQ);

%% Determine borders
% Zero-crossings and negative values are considered part of the final mask,
% whereas positive values lie outside the mask

interpMask = (VQ + abs(VQ)) == 0;
