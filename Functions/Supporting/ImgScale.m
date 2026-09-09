%Function to show image scaled by percentage
%INPUTS - 2D non-complex image
%       - Percentage lower bound
%       - Percentage upper bound
function ImgScale(img, lower, upper)
    c_scale = [prctile(img,lower,'all'), prctile(img,upper,'all')];
    imshow(img, c_scale)
end
