%% Initialization
clc;
close all;
clear;

%% Read the data
gr_tr = load('groundtruth_rect.txt');
gr_tr(:,1:2) = gr_tr(:,1:2)-10; 
gr_tr(:,3:4) = gr_tr(:,3:4)+30; 

files = dir('./Training/*.jpg');

aspect_ratio_plot = gr_tr(:,4)./gr_tr(:,3);

aspect_ratio      = mean(aspect_ratio_plot);

% <to-do>
training_pos      = zeros(20, (146*62));
training_neg      = zeros(100, (146*62));
HOG_pos_ex        = zeros(20, 3672);
HOG_neg_ex        = zeros(100, 3672);

PCA_pos_ex        = zeros(20, 500);
PCA_neg_ex        = zeros(100, 500);

temp_index        = 1;
%% Training positive
for mat_idx = 1:5:size(files, 1)
    file_name = strcat('./Training/', files(mat_idx).name);
    image     = imread (file_name);
    if (size(image, 3) == 3)
        image   = rgb2gray(image);
    end
    
    image = imcrop(image, gr_tr(mat_idx,:));
    image = imresize(image, [146, 62]);
    HOG_pos_ex(temp_index,:) = extractHOGFeatures(image);
    temp_index = temp_index+1; 
    
%     training_pos(mat_idx,:) = image(:);
%     figure;
%     imshow (image);
%     imshow(image);
%     hold on;
%     rectangle('Position', gr_tr(mat_idx,:), 'EdgeColor', 'r');
%     hold off;
    
end

%% Train negative
file_name = strcat('./Training/', files(1).name);
image_for_neg_train  = rgb2gray(imread (file_name));

for mat_idx = 1:100
    rand_num_x = randperm((288-146),1);
    rand_num_y = randperm((384-62),1);
    image = imcrop(image_for_neg_train,[rand_num_y rand_num_x ... 
                            61 145]);
    HOG_neg_ex(mat_idx,:) = extractHOGFeatures(image);
                        
%     figure;
%     imshow (image);
%     training_neg(mat_idx,:) = image(:);
end


%% Testing begins
[training_labels{ 1:20}] = deal('human');
[training_labels{21:120}] = deal('no_human');

clubbed_train_data = vertcat(HOG_pos_ex, HOG_neg_ex);

[PCA_norm, ~, ~] = featureNormalize(clubbed_train_data);

[U, ~]    = pca_manual(PCA_norm);

PCA_train_data = projectData(clubbed_train_data, U, 100);

svm_model = fitcsvm(clubbed_train_data, training_labels);

%%  Let us test on the 1st image

for img_idx = 1: 100
    
    file_name = strcat('./Training/', files(img_idx).name);
    image_for_gen_test  = rgb2gray(imread (file_name));
    
    scale_factor = 1;
    found_human  = 0;
    
    disp (img_idx);
    
    while ((scale_factor < 1.20) && (found_human == 0))
        
        image_for_gen_test   = imresize(image_for_gen_test, scale_factor);
        [im_height,im_width] = size(image_for_gen_test);
        
        w_pt = 1;
        h_pt = 1;
        
        %store_val = [1 1 1 1];
        while (((h_pt+146) <= im_height) && (found_human == 0))
            w_pt = 1;
            while (((w_pt+62) <= im_width) && (found_human == 0))
                image_crop = imcrop(image_for_gen_test, [w_pt h_pt 61 145]);
                
                HOG_feature = extractHOGFeatures(image_crop);
                
                result = predict(svm_model, HOG_feature);
                
                %disp (result);
                
                if (strcmp(result, 'human'))
                    
                    store_val = [h_pt, w_pt, 62 146];
                    found_human = 1;
                    imshow (image_for_gen_test);
                    hold on;
                    rectangle ('position', store_val, 'EdgeColor', 'r', 'LineWidth', 4);
                    saveas(gcf,['.\results\' num2str(img_idx) '.png']);
                end
                
                w_pt = w_pt + 30;
            end
            h_pt = h_pt + 30;
        end
        
        if (found_human == 0)
            scale_factor = scale_factor + 0.02;
        end
        
    end
end

   
%% ONLY IF I NEED TO VERIFY SVM MODEL
% for mat_idx = 1: 100
%     result = predict(svm_model, HOG_neg_ex(mat_idx,:));
%     disp(result);
% end

