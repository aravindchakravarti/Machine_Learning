% --------- PARTICLE SWARM OPTIMIZATION FOR IMAGE SEGMENTATION ------------
%
% File Name     : PSO_Otsu.m
%
% Description   : Image segmention technique is primay technique used
%                 image processing. It 'segments' the foreground and
%                 background. To segment the image, we need to find the
%                 threshold value, which will suite optimally for an image.
%                 OTSU proposed such method (1970s) that says "The optimum
%                 threshold value can be found either by 'increasing the
%                 between class variance' or by 'decreasing the in class
%                 varience'". In the implementation, the former one
%                 is used. 
%
%                 Finding the one threshold is a uni-dimension problem, but
%                 having multiple threshold is a multi-dimension problem.
%                 It is a optimization (maximization) problem. This
%                 problem has been solved using 'Particle Swarm
%                 Optimization' paradigm of CI
%
% Reference     : PT - 2017 MIS501 slides and codes
%                 Author: Dr. R. V. Kulkarni
%                 M. S. Ramaiah University of Applied Sciences, Bengaluru
%
%                 A Threshold Selection Method from Gray-Level Histograms 
%                 IEEE Transactions on Systems, Man, and Cybernetics, 
%                 Vol.9, No. 1, 1979, pp. 62-66.
%                 Author: Noboyuki Otsu
%
% File History  : | Ver. |  Date     |              Changes               |
%                 | 1.0  | 03/AUG/17 | Base version of PSO                |
%                 |      |           |                                    |
%                 | 1.1  | 01/SEP/17 | Integrated OTSU fitenss fucntion   |
%                 |      |           |                                    |
%                 | 1.2  | 02/SEP/17 | Fixed errors in PSO for image proc.|
%                 |      |           |                                    |
%                 | 1.3  | 03/SEP/17 | Added stopping criterion           |
% -------------------------------------------------------------------------

% Clear the clutter
clc;            % Clean the screen
clear;          % Remove all run time variables
close all;      % Any window open? Close them!

% Initilizations 

% Our search boundary. This is same as min and max of thresholds
thresh_max = 256;
thresh_min = 001;

% Our maximum velocities
v_max = +120;
v_min = -120;

% Charecteristics of PSO
wt          = 0.5;     % These are standard values recommended
c_1         = 2.0;
c_2         = 2.0;
max_iter    = 400;
decrs_wt    = 1;

% Number of particles we need and there place in fitness landscape
no_particles    = 40;
thresh_levels   = 08;
no_of_runs      = 10;

% These variables are for our statistical analysis
thresh_stats     = zeros (no_of_runs, thresh_levels);
timer_stop       = zeros (1,no_of_runs);
% For plotting
gBest_fitness_plot  = zeros (no_of_runs, max_iter);

% Lets try to initialize our image related things
rgb_image   	= imread ('Lena_std.tif');
% Let us convert image to grayscale
image           = rgb2gray (rgb_image);
image_size      = size (image);
n_pixels        = image_size(1) * image_size(2);

% Hiistogram related information extraction
image_hist      = imhist (image);
image_hist_p    = image_hist./n_pixels;
pixel_max_level = size (image_hist_p,1);

% Total mean level of the image
mau_t           = 0;
for i = 1:pixel_max_level
    mau_t = mau_t + (i*image_hist_p(i));
end

% Our stopping criteria
last_fourty_fitness         = zeros (40,1);
enable_stopping_criterion   = 1;

% Number of trials we would like to run
for run=1:no_of_runs
    
    % To identify how much MATLAB takes in executing the code
    tic;
    
    % Randomly allocate them in space with different velocities
    % We are rounding them because, we cannot have decimal values in
    % threshold values. Threshold values are always discrete
    particle_pos  = round (thresh_min + ((thresh_max-thresh_min)*...
                    (rand(no_particles,thresh_levels))));
    particle_vel  = v_min + ...
                    ((v_max-v_min)*(rand(no_particles,thresh_levels)));
    
    % Because thresholding requires always sorted data, 
    % let sort the particles in assending order.
    particle_pos  = sort (particle_pos,2);
    
    % For the particles, current best fitness is
    pBest_fitnesses     = zeros (1,no_particles);
        
    % Call, Otsu function to determine the fitness of each particle to know
    % what is best 'between class variance' we get.
    for i = 1:no_particles
        pBest_fitnesses(i) = tryOtsu(image_hist_p, thresh_levels, ...
                                     particle_pos(i,:), mau_t);
    end
    
    % Out of all fitness, this particle had minimum fitness
    [gBest_fitness, index]  = max (pBest_fitnesses);
    pBest_positions         = particle_pos;
    gBest_position          = particle_pos(index,:);
    
    iteration               = 1;
    
    stop_this_exec          = 0;
    
    while ((iteration < max_iter) && (stop_this_exec ~= 1))
        
        iteration = iteration + 1;
        
        % Let us scan through all the particles
        for index = 1: no_particles
            
            % If decreasing weight is allowed
            if decrs_wt == 1
                wt = 0.5 - (iteration*(0.4/max_iter));
            end
            
            % This is the velocity influence of 'current' particle itself
            inertial_vel = wt * particle_vel (index,:);
            
            % This is the congignet velocity; 
            % Influence of 'particles own best position'
            cognitv_vel = c_1*rand*(pBest_positions(index,:)-...
                                    particle_pos(index,:));
            
            % This is the social velocity 
            % Influence of 'global best particle'
            social_vel  = c_2*rand*(gBest_position-particle_pos(index,:));
            
            % For current particle, net influence of velocity is 
            velocity_factor = inertial_vel + cognitv_vel + social_vel;
            
            % Throttling the velociy component
            velocity_factor (velocity_factor > v_max) = v_max;
            
            velocity_factor (velocity_factor < v_min) = v_min;
            
            % Update velocity factor for the particle
            particle_vel(index, :) = velocity_factor;
            
            % Current particles new position.
            % p(p1, p2, ... ) = p(p1, p2, ... ) + (velocity_factor * 1sec)
            particle_pos (index, :) = round (particle_pos(index, :) + ...
                                      particle_vel(index,:));
          
            % Throttling position components
            particle_pos (index,(particle_pos(index, :) > thresh_max))= ...
                                                                thresh_max;
            particle_pos (index,(particle_pos(index, :) < thresh_min))= ...
                                                                thresh_min;
            
            % Before we ESTIMATE the fitness, we need sort the values.
            % Beccause OTSU criteria. It is also important that we also 
            % exchange velocities associalted with each particle.
            
            % BUBBLE SORTING
            for sort_indx_1 = 1 : thresh_levels
                for sort_indx_2 = 1 : thresh_levels - sort_indx_1
                    % If particle heigher than next one, exchange it
                    if (particle_pos(index,sort_indx_2)) > ...
                                      (particle_pos(index,(sort_indx_2+1)))
                                  
                        temp_store_1 = particle_pos(index,sort_indx_2);
                        temp_store_2 = particle_pos(index,(sort_indx_2+1));
                        
                        particle_pos(index,sort_indx_2) = temp_store_2;
                        particle_pos(index,(sort_indx_2+1)) = temp_store_1;
                        
                        % If we exchanged the particle, then we need to
                        % exchange velocities associated with it too
                        temp_store_1 = particle_vel(index, sort_indx_2);
                        temp_store_2 = particle_vel(index,(sort_indx_2+1));
                        
                        particle_vel(index, sort_indx_2) = temp_store_2;
                        particle_vel(index,(sort_indx_2+1)) = temp_store_1;
                        
                    else
                        % Else do nothing
                    end
                end % End of sort_index_2
            end % End of sort_index_1
            
            % Let us calculate the fitness of particle now
            current_fitness = tryOtsu(image_hist_p, thresh_levels, ...
                                      particle_pos(index,:), mau_t);
            
            % Update pBest and gBest components. 
            if current_fitness > pBest_fitnesses(index)
                pBest_positions(index,:) = particle_pos (index,:);
                pBest_fitnesses(index)   = current_fitness;
            end
            
            if current_fitness > gBest_fitness
                gBest_position = particle_pos (index, :);
                gBest_fitness  = current_fitness;
            end
            
        end % End of for (number of particles)
        
        % This is useful for plotting
        gBest_fitness_plot(run, iteration) = gBest_fitness;
        
        % If there is no improvement in 40 iterations, let us not waste
        % time :)
        if (enable_stopping_criterion == 1)
            if (iteration > 40)
                last_fourty_fitness(iteration-40) = gBest_fitness;
                if all(last_fourty_fitness == last_fourty_fitness(1))
                    stop_this_exec = 1;
                    % Without this, graph looks ugly :)
                    gBest_fitness_plot(run,(iteration:max_iter))...
                                                           = gBest_fitness;
                end
            end % End of if (iteration > 40)
        end % End of enable_stopping_criterion
        
    end % End of while (iteration < max_iter)
          
    % fprintf ('The best position which this algorithm achieved is: \n');
    disp(gBest_position);
    
    thresh_stats(run,:) = gBest_position;
    
    timer_stop(run) = toc;
    
end % End of number of runs

disp ('Mean time taken for execution');
disp (mean (timer_stop));

disp ('Standard deviation for execution');
disp (std(timer_stop));

% If required, we can also plot fitness achieved in each position
figure;
plot (gBest_fitness_plot');
title ('Fitness function, Iteration Vs Number of trials');
xlabel ('Iterations');
ylabel ('Fitness function');

% Next section is to print the thresholded 'Lena'

thresh_image = image;

% If we have 'N' thresholds, then we will have N+1 levels in image
for i = 1:(thresh_levels+1)
    % What should be mmy thresold values?
    if i == 1
        start_thresh = 1;
        stop_thresh  = gBest_position(i);
        gray_level   = 0;
    % If it is last one
    elseif i == (thresh_levels+1)
        start_thresh = gBest_position(i-1)+1;
        stop_thresh  = 256;
        gray_level   = 256;
    else
        start_thresh = gBest_position(i-1)+1;
        stop_thresh  = gBest_position(i);
        gray_level   = (start_thresh+stop_thresh)/2;
    end
    
    % Thresholding
    for pixels = 1: n_pixels
        if (image(pixels) >= start_thresh) && (image(pixels) <= stop_thresh)
            thresh_image(pixels) = gray_level;
        end
    end
end

mean_of_trials = mean (thresh_stats);
std_of_trials  = std (thresh_stats);

if no_of_runs > 1 
    disp ('Means:');
    disp (mean_of_trials);
    disp ('Standard Deviation:');
    disp (std_of_trials);
end

% Printing result here
figure;
subplot(1,2,1);
imshow (image);
title ('Original Image');

subplot (1,2,2);
imshow (thresh_image);
title ('3 - Level Segmented Image');