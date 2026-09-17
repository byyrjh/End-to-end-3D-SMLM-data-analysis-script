function img=threeDimgread(imgpath,stack_idx,num_slice)
obj = Tiff(imgpath,'r');
dimx = obj.getTag('ImageWidth');
dimy = obj.getTag('ImageLength');
% dimz = obj.getTag('PageNumber');   % # of Slices
img = zeros(dimy,dimx,num_slice,'uint16');
for i = 1:num_slice
    obj.setDirectory((stack_idx-1)*num_slice+i);
    img(:,:,i) = obj.read();
end
end