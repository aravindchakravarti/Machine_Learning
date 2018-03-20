%% File Information                         
% File Name             : HistOfOrientGrad.m
%
% Description           : This file implements HOG (Histogram of Oriented 
%                         Gradients
%
% References            : https://www.learnopencv.com/histogram-of-oriented
%                        -gradients/
%
% Author(s)             : Aravind D. Chakravarti
%
% Version History       :
% Ver   Name        Change Description
% 1.0   Aravind D   Started with basics. Just a trial
%
%

function [features] = HistOfOrientGrad (I) %#ok<INUSD>
%% Initial Stuff; Delete this when function is done <TO-DO>
clc;
close all;
clear;

%% Inbuilt command HOG for trial
img = imread('cameraman.tif');
[~, features] = extractHOGFeatures(img); 
figure;
imshow(img); hold on;
plot(features);
title ('HOG by Built-in Function');

%% Let us start Manual implementation

% For now, read some dummy image
I_raw = double (imread('cameraman.tif'));
%[im_row,im_col] = size(I_raw);
I = applyGaussianFilter(I_raw);

% This point forward, I am sticking here with desription provided by 
% https://www.learnopencv.com/histogram-of-oriented-gradients/

%% STEP - I : Preprocessing
% This step is skipped as of now. Since we are not tracking any object here
% <TO-DO>

%% STEP - II : Calculate the Gradient Images
% We are going to run SOBEL operator on the image to calculate x and y
% gradients of intensities, and angle of gradients. 
[I_Grad, I_Directn] = runSobelOnImage (I);
% figure;
% imshow (uint8(I_Grad));
% title ('Sobel Operated Image');

%% STEP - III : Calculate Histogram of Gradients in 8×8 cells 
% In this step, the image is divided into 8×8 cells and a histogram of 
% gradients is calculated for each 8×8 cells.
% Note: The webpage talks about implementing HOG for RGB image. Hence he
% tells that image patch contains 8x8x3 = 192 pixels. Since we are
% implementing for Grayscale we will have 8x8 = 64 pixels


end % End of HistOfOrientGrad Function

function [I] = applyGaussianFilter(image)

% Run Guassian filter to remove un-wanted noise
G = [2,  4,   5,  4,  2;
    4,  9,  12,  9,  4;
    5, 12,  15, 12,  5;
    4,  9,  12,  9,  4;
    2,  4,   5,  4,  2];

G = (1/159).*G;

I = image;
[im_row,im_col] = size(I);

% Convolution of image with Guassian filter
for i = 3 : (im_row-2)
    for j = 3 : (im_col-2)
        test_location   = image(i-2:i+2,j-2:j+2);
        I(i,j)          = sum(sum(G.*test_location));
    end
end

end

function [I_Grad, I_Directn] = runSobelOnImage (I)
% Sobel X direction operator
G_x = [ -1,  0,  1;
        -2,  0,  2;
        -1,  0,  1];

% Sobel Y direction operator
G_y = [ 1,  2,  1;
        0,  0,  0;
       -1, -2, -1];
   
[im_row,im_col] = size(I);

% These two variables hold the magniture and direction of image intensities
I_Grad    = I;
I_Directn = uint8(zeros(size(I)));

for i = 2 : (im_row-1)
    for j = 2 : (im_col-1)
        % Calculating Sobel gradients' magnitude
        test_location       = I(i-1:i+1,j-1:j+1);
        x_grad              = sum(sum(G_x.*test_location));
        y_grad              = sum(sum(G_y.*test_location));
        
        % Magnitude
        I_Grad(i,j)         = sqrt(x_grad.^2 + y_grad.^2);
        % Angle (atan2(y_grad, x_grad))*(180/pi)
        I_Directn(i,j) 		= uint8(atan2(y_grad, x_grad)*(180/pi));
        if (I_Directn(i,j)  < 0)
            I_Directn(i,j)  = I_Directn(i,j) + 180;
        end
%         if (I_Directn(i,j) > 180)
%             I_Directn(i,j) = I_Directn(i,j) - 180;
%         end
        
    end
end

end


