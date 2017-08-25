clc;
close all;

result = zeros (256,1);

rgb_image   	= imread ('fingerprint.jpg');
image           = rgb2gray (rgb_image);
image_size      = size (image);
n_pixels        = image_size(1) * image_size(2);

for threshold_level = 1 : 256
    
    image_hist      = imhist (image);
    image_hist_p    = image_hist./n_pixels;
    pixel_max_level = size (image_hist_p,1);
    
    mau_t           = 0;
    mau_k           = 0;
    omega_k         = 0;
    
    for i = 1:threshold_level
        omega_k = omega_k + image_hist_p(i);
    end
    
    for i = 1:threshold_level
        mau_k = mau_k + (i*image_hist_p(i));
    end
    
    for i = 1:pixel_max_level
        mau_t = mau_t + (i*image_hist_p(i));
    end
    
    sigma_sqr_b_k   = (((mau_t*omega_k)-mau_k)^2)/(omega_k*(1-omega_k));
   
    result (threshold_level) = sigma_sqr_b_k;
end

plot (result);

max_thresh = max (result);

idx = mean(find(result == max_thresh));

level = (idx-1)/(pixel_max_level-1);

disp ('From MATLAB IN BUILT PROGRAM');
disp (graythresh (image));

disp ('From OUR PROGRAM');
disp (level);

BW_IMG = im2bw(image,level);
imshow (BW_IMG);
