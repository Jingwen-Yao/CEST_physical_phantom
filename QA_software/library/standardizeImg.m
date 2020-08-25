function [img_256,padNumber] = standardizeImg(img)
% resize and pad to 256 x 256

% remove dimensions of 1
img = squeeze(img);

% calculate scale based on max dimension
matSize = size(img);
maxSize = max(matSize);
scale = 256/maxSize;

% resize so that largest dimension is 256
img_resize = imresize(img,scale);

% calculate number of padding arrays
matSize = size(img_resize);
pad_row1 = ceil((256 - matSize(1))/2);
pad_row2 = floor((256 - matSize(1))/2);
pad_col1 = ceil((256 - matSize(2))/2);
pad_col2 = floor((256 - matSize(2))/2);

% pad the array
img_pad = padarray(img_resize,[pad_row1 pad_col1],0,'pre');
img_256 = padarray(img_pad,[pad_row2 pad_col2],0,'post');
padNumber = [pad_row1 pad_row2 pad_col1 pad_col2];

end