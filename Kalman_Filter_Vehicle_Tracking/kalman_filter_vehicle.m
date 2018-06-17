%% Kalman Filter
%--------------------------------------------------------------------------
% Description: 
% This file implements the Kalman filter to estimate the position of the 
% vehicle moving on the track.
% A super simplistic model is built to understand & analyze the filter.
% 
% Kalman filter is has 2 steps
%                ______________________
%               ^                      |
%           prediction              Correction
%               |______________________|
%
% Predition step is the apriori estimate made by the filter, and 
% correction step is the correction applied to apriori estimations based
% on the observations.
% 
% In this program,
% Any apriori estimates (We call it as prediction) will begin with predct_
% Any estimation will begin with estimd_
%                                                    _
% predct_X = (A * estimd_X) + B                       |-> Prediction state
% predct_P = (A * estimd_P * A') + Q                  |
%                                                    _
% K = (C*predct_P)/{(C*predct_P*C')+R}                |
% estimd_X = predct_X + {K * [Z - (C * predct_X)]}    |-> Correction state
% estimd_P = (1-K*C) predict_P                        |
%
%
% Author    : Aravind D. Chakravarti
%             aravind dot chakravarti at gmail dot com
%
% References: (1) Notes of Dr. Hariharan Ramasangu 
%             (2) http://studentdavestutorials.weebly.com/
%             (2) Understanding the Basis of the Kalman Filter Via a Simple 
%                 and Intuitive Derivation - Ramsey Faragher
% 
% Version history:
% 
% Ver.    Date                             Changes
%--------------------------------------------------------------------------
% 1.0   14-June-2018    Started implementation
% 1.1   16-June-2018    Restarted everything again
% 
% 
%--------------------------------------------------------------------------

% Some basic things for MATLAB
clc; clear; close all;

% Time of flight: We are going to measure for 10 seconds
TOF =  10;
% Update rate
dt =   0.1;
% Time instances
time_inst = 0:0.1:TOF;

% A = State transition matrix; B = Control input matrix; 
% C = Transformation matrix
A = [1     dt; 0    1];
B = [(dt^2/2);     dt]; 
C = [1     0];

% U = Input matrix
U = 1.5;

proces_noise = 0.05;
measur_noise = 10;

%% Data Generation
% Below variables are declared which are used for data generation purpose
% only

ideal_travel    = zeros(length(time_inst), 2);
state_data      = zeros(length(time_inst), 2);
observd_data    = zeros(length(time_inst), 1);

idx = 1; 
for t = 0.1:0.1: TOF
    % Ideal data is used only for plotting purpose. Refer file header
    % section for more information on below equation
    ideal_travel(idx+1,:) = (A*ideal_travel(idx,:)') + (B*U);
    state_data(idx+1,:)   = (A*state_data(idx,:)') + (B*U) + ...
                                                proces_noise*[randn;randn];
    observd_data(idx+1,:) = C*state_data(idx+1,:)' + ...
                                                measur_noise*randn;
    idx = idx+1;
end

figure;
subplot (1, 2, 1);
plot (time_inst, ideal_travel(:,1), 'b', 'LineWidth', 2);
hold on;
plot (time_inst, observd_data, 'r');
title ('Synthetic Data Generation');
hold off;

%% Begin of Kalman filter
estimd_X    = [0; 0];
%estimd_P    = [randn, randn; randn, randn];
%Q           = [randn, randn; randn, randn];
%R           = randn;
Kalman_out  = zeros(length(time_inst), 2);
R           = measur_noise^2;
Q           = (proces_noise^2) * [dt^4/4, dt^3/2; dt^3/2, dt^2];
estimd_P    = Q;

idx = 1;
for t = 0:0.1: TOF
    % Prediction phase--------------------------------------------------
    predct_X = (A * estimd_X) + (B*U);                       
    predct_P = (A * estimd_P * A') + Q;
    
    % Correction phase--------------------------------------------------
    % Calculating Kalman Gain
    K = (predct_P*C')/((C*predct_P*C')+R);
    % How much do we need to correct our prediction?
    estimd_X = predct_X + (K * (observd_data(idx) - (C * predct_X)));
    estimd_P = (eye(2)-K*C)*predct_P;
    
    % For plotting purpose----------------------------------------------
    Kalman_out(idx, :) = estimd_X;
    idx = idx+1;
end

subplot (1, 2, 2);
plot (time_inst, ideal_travel(:,1), 'b', 'LineWidth', 2);
hold on;
plot (time_inst, Kalman_out(:,1), 'r');
hold off;
title ('Kalman Filter Output');

%% Things TO-DO
% Why would we relate Q with time? 
% Why to set estimated P with prior known value?
