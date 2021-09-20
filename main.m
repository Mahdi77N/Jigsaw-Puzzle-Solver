clc;
clear;

% ************************* First, run cvx_setup *************************

% ************************ Enter the path manually ************************
patches_path = 'C:\Users\Mahdi\Desktop\Mahdi Naderi - 9636713\Puzzle_1_160\';
% *************************************************************************

patches = dir(patches_path);

proc_img = im2double(imread([patches_path patches(8).name]));
corner = imread([patches_path patches(3).name]);

H = size(proc_img, 1);	% height of output image
W = size(proc_img, 2);	% width of output image

h = size(corner, 1);	% height of patches
w = size(corner, 2);	% width of patches

rows = H / h;
cols = W / w;
blocksize = h;

cnt = 3;
for i=1:rows
	for j=1:cols
		tmp = imread([patches_path patches(cnt).name]);
		proc_img((i-1)*h+1:i*h, (j-1)*w+1:j*w, :) = im2double(tmp);
		cnt = cnt+1;
		if cnt == 7
			cnt = 9;
		end
	end
end

numOfPatches = rows * cols;
rgbPatches = zeros([blocksize, blocksize, 3, numOfPatches], class(proc_img));

for i = 1:numOfPatches
	rowStart = (ceil(i/cols)-1) * blocksize + 1;
	rowEnd = rowStart + (blocksize-1);
	colStart = mod(i-1, cols)  * blocksize + 1;
	colEnd = colStart + (blocksize-1);
	rgbPatches(:, :, :, i) = proc_img(rowStart:rowEnd, ...
		colStart:colEnd, :);
end

partsOrder = randperm(numOfPatches);
imageblocks = rgbPatches(:, :, :, partsOrder);

result = linearProg(imageblocks, rows, cols);

orig_img = im2double(imread([patches_path patches(7).name]));

fprintf('*************************************\n');
fprintf('MSE = %f\n', immse(result, orig_img));
fprintf('*************************************\n\n');
