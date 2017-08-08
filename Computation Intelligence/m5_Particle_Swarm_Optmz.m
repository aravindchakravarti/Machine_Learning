% -------------------- PARTICLE SWARM OPTIMIZATION ------------------------
%
% File Name     : m5_Particle_Swarm_Optmz.m
%
% Description   : This file implements PSO on sphere function. It is a
%                 minimization problem
%
% Author        : Aravind D. Chakravarti
%                 M. S. Ramaiah University of Applied Sciences, Bengaluru
%
% Reference     : PT - 2017 MIS501 slides and codes
%
% Version       : 1.0       03/AUG/2017
%                 Base Version
% -------------------------------------------------------------------------

% Clear the clutter
clc;
clear;
close all;

% Initilizations 

% Our search boundary
x_max = +5;
x_min = -5;

% Our maximum velocities
v_max = +1;
v_min = -1;

% 'Variables' of PSO
w  = 0.5;
c1 = 2.0;
c2 = 2.0;
max_iterations = 1000;
decrease_weight= 1;

% How many particles we need? And where they rest in search space
no_of_particles = 40;
dimension       = 10;

% Randomly allocate their space and speed (velocity!)
particle_positions = x_min + ((x_max-x_min)*(rand(no_of_particles,dimension)));
particle_velocities = v_min + ((v_max-v_min)*(rand(no_of_particles,dimension)));

% This is our fitness function
pBest_fitnesses = sum ((particle_positions.^2),2);

% Out of all fitness, this particle had minimum fitness
[gBest_fitness, index] = min (pBest_fitnesses);
pBest_positions = particle_positions;
gBest_position  = particle_positions(index,:);

% These variables are just used for plotting
pBest_plot = zeros (max_iterations,1);
gBest_plot = zeros (max_iterations,1);

iteration = 1;

while iteration < max_iterations
    
    iteration = iteration + 1;
    
    % Let us scan through all the particles
    for index = 1: no_of_particles
        
        % If decreasing weight is allowed
        if decrease_weight == 1
            w = 0.5 - (iteration*(0.4/max_iterations));           
        end
        
        % This is the velocity influence of 'current' particle itself
        inertial_vel = w * particle_velocities (index,:);
        
        % This is the velocity influence of 'current particles best
        % position'
        cognitv_vel = c1*rand*(pBest_positions(index,:)-particle_positions(index,:));
        
        % This is the velocity influence of 'global best particle'
        social_vel  = c2*rand*(gBest_position - particle_positions(index,:));
        
        % For current particle the net influence of velocity is given by
        velocity_factor = inertial_vel + cognitv_vel + social_vel;
        
        % Throttling the velociy component
        velocity_factor (velocity_factor > v_max) = v_max;
        
        velocity_factor (velocity_factor < v_min) = v_min;
        
        % Update velocity factor for the particle
        particle_velocities(index, :) = velocity_factor;
               
        % Current particles new position.
        % p(x1, x2, ... ) = p(x1, x2, ... ) + (velocity_factor * 1s)
        particle_positions (index, :) = particle_positions(index, :) + particle_velocities(index,:);
        
        % Throttling position components
        particle_positions ((particle_positions(index) > x_max),:) = x_max;
        
        particle_positions ((particle_positions(index) < x_min),:) = x_min;
        
        % This is my fitness function
        current_fitness = sum (particle_positions(index,:).^2);
        
        % Update pBest and gBest components
        if current_fitness < pBest_fitnesses(index)
            pBest_positions(index,:) = particle_positions (index,:);
            pBest_fitnesses(index)   = current_fitness;
        end
        
        if current_fitness < gBest_fitness
            gBest_position = particle_positions (index, :);
            gBest_fitness  = current_fitness;
        end
        
        % Updating plotting variables
        pBest_plot(iteration) = pBest_fitnesses(index);
        gBest_plot(iteration) = gBest_fitness;
               
    end
   
end


% Plot the particles for visualization. 
% -- CAUTION : HOLDS GOOD ONLY FOR TWO DIMENSIONS !!!! ----
plot (pBest_positions(:,1), pBest_positions(:,2), 'o');
        
figure;
% Ploting the results
subplot (2,1,1);
plot (pBest_plot);
title ('pBest positions');

subplot (2,1,2);
plot (gBest_plot);
title ('gBest positions');

% Display few things on command window
fprintf ('The best position which this algorithm achieved is: \n');
disp(gBest_position);
fprintf ('Fitness value %e\n', gBest_fitness);
