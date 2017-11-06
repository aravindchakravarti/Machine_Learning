%% -------------------- Function Approximation ----------------------------
%
% File Name     : Mackey_glass.m
%
% File Info.    : This file includes implementation of 5x5x3x1 neural
%                 network. Network is trained to solve Mackey-Glass series.
%                 Back-propagation algorithm is being used for training.
% 
% Version       : 1.0
%      Author   : Aravind D. Chakravarti           
%
% References    : Codes and Lecture Slides of MIS502
%      Author     Dr. R. V. Kulkarni
%                 MSRUAS, Bangalore.
% -------------------------------------------------------------------------

%% Neural Network
%{
            THIS IS HOW NEURAL NETWORK LOOKS LIKE!
            ___         ___          ___
           /   \       /   \        /   \
   X ----> |   |  ...  |   | .....  |   |
           \___/       \___/        \___/

            ___         ___          ___           ___       
           /   \       /   \        /   \         /   \      _______
   1 ----> |   |  ...  |   | .....  |   | ......  |   | ---->|  J   |---> e
           \___/       \___/        \___/         \___/      |______|
                                                                ^
                        ___                                     |
                       /   \                                    |
                       |   | .....                           desired
                       \___/                           
---------------------------------------------------------------------------
            a1     |    a2     |     a3       |     a4      |
                 theta_1    theta_2        theta_3
      X            | a2_i      |a3_i          |a4_i         |
               a1_o|       a2_o|         a3_o |        a4_o |
                   |           |              |             |
---------------------------------------------------------------------------
%}

%% Clear the clutter from MATLAB windows
clc; close all; clear;

%% Initialize the Neural Network paramaeters
% Initialize number of units in each layer, number of training set 
% Read Mackey Glass data and display them to visualize
inp_layer_units     =   4;
hid_lay_1_units     =   5;
hid_lay_2_units     =   5;
out_layer_units     =   1;

m                   =   30000; 
learning_gain       =   0.025;
momentum_gain       =   0.975;

mg_samples  = load ('mgdata.dat');
mg_data     = mg_samples(:,2);

% Random initialization of network weights. sqrt(6)/layers gives the best
% apporoximation for getting the results.
epsilon_1   = sqrt (6)/(sqrt((hid_lay_1_units)+(inp_layer_units)));
epsilon_2   = sqrt (6)/(sqrt((hid_lay_2_units)+(hid_lay_1_units)));
epsilon_3   = sqrt (6)/(sqrt((out_layer_units)+(hid_lay_2_units)));

%% Randomize the weights of synapses 
theta_1             = rand (hid_lay_1_units, inp_layer_units+1) * 2 ...
                      * epsilon_1  - epsilon_1;
theta_2             = rand (hid_lay_2_units, hid_lay_1_units) * 2 ...
                      * epsilon_2  - epsilon_2;
theta_3             = rand (out_layer_units, hid_lay_2_units) * 2 ...
                      * epsilon_3  - epsilon_3;

%% Neural Network Training Begins here. It has two steps; 
%  1. Forward propogation
%  2. Backward (Back) propogation

error_exemplar  = zeros (600,1);

error_training  = zeros (m,1);

input_X         = zeros (inp_layer_units+1, 1);

for i = 1: m
    
    delta_theta_2 = 0;
    delta_theta_1 = 0;
    delta_theta_3 = 0;
    
    fprintf ('Current iteration %d\n', i);
      
    rand_int = randi([19 600],1, 50);
    
    for j = 1: 50
              
        sel_rand_int   = rand_int(j);
        
        % Input Neurons; they just pass what they got.
        input_X(1)     = 1;
        input_X(2)     = mg_data(sel_rand_int);
		input_X(3)     = mg_data(sel_rand_int-6);
		input_X(4)     = mg_data(sel_rand_int-12);
		input_X(5)     = mg_data(sel_rand_int-18);
        
        a1_o        = input_X;

        % Hidden Neurons; inputs to them will be weighted sum of input
        % layer outputs and ouput will be sigmoid output 
        a2_i        = theta_1 * a1_o;
        a2_o        = 1./(1+exp(-a2_i));
        
        a3_i        = theta_2 * a2_o;
        a3_o        = 1./(1+exp(-a3_i));
   
        % Output Neuron; this is a linear Neuron. It just combines what it
        % got at its input
        a4_i 		= theta_3 * a3_o;
        a4_o		= a4_i;
         
        % For our cost function; the desired output what we expected is;
        d           = mg_data(sel_rand_int+6);
        
        error       = d-a4_o;
        
        error_exemplar(j) = error;
        
        % 2. Back propogation begins here
        
        % The formulae here come from partial differentiation and chain
        % rule; doo_* here is representing partial differentiation
        
        doo_thetha_3 = error*1*a3_o';
        
        doo_thetha_2 = ((error*1*theta_3').*(a3_o.*(1-a3_o)))*a2_o';
        
        doo_thetha_1 = ((((error*1*theta_3').*(a3_o.*(1-a3_o)))'*...
                         theta_2)'.*(a2_o.*(1-a2_o)))*a1_o';
          
        % We don't want to get stuck in local minima; hence we are going
        % for momentum based; where, momentum represents kind of
        % acceleration, so that it will pull us out of local minima if 
        % we get stuck there..
        delta_theta_3 = (learning_gain * doo_thetha_3) + ...
                        (momentum_gain * delta_theta_3);
        
        delta_theta_2 = (learning_gain * doo_thetha_2) + ...
                        (momentum_gain * delta_theta_2);
                    
        delta_theta_1 = (learning_gain * doo_thetha_1) + ...
                        (momentum_gain * delta_theta_1);
    end
    theta_3 = theta_3 + delta_theta_3;
    theta_2 = theta_2 + delta_theta_2;
    theta_1 = theta_1 + delta_theta_1;
    
    error_training(i) = mean(error_exemplar.^2);
end

%% Testing phase

% Let us have more samples than what we trained so far

mg_testing          = zeros (size(mg_data,1));
mg_testing          = mg_data (1:606,:);

test_input_X        = zeros (inp_layer_units, 1);

for time = 600 : (1201-6);

    test_input_X(1)     = 1;
    test_input_X(2)     = mg_testing(time);
	test_input_X(3)     = mg_testing(time-6);
	test_input_X(4)     = mg_testing(time-12);
	test_input_X(5)     = mg_testing(time-18);
  
    test_a1_o        = test_input_X;
    
    test_a2_i        = theta_1 * test_a1_o;
    test_a2_o        = 1./(1+exp(-test_a2_i));
    
    test_a3_i        = theta_2 * test_a2_o;
    test_a3_o        = 1./(1+exp(-test_a3_i));
    
    test_a4_i 		 = theta_3 * test_a3_o;
	test_a4_o		 = test_a4_i;
    
    mg_testing(time+6)= test_a4_o;
        
end

desired_out = mg_data;

% Plotting the results
figure;
plot (error_training);
xlabel ('Trial data');
ylabel ('MSE');
title  ('Error in each Epoch');

figure;
plot (desired_out)
hold on;
plot (mg_testing, 'r');
title ('Mackey-Glass Time Series');
xlabel ('Time');
ylabel ('f(x)');
legend ('Desired Output', 'ANN Output');
hold off;

time_series_prediction_err = desired_out - mg_testing;

fprintf ('The mean square error of predicted output is %f', ...
        mean (time_series_prediction_err.^2));

