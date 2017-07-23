% Implementing three unimodal and three multimodal benchmark functions
% 
% Author : Aravind D. Chakravarti
%          PT-2017
%          Machine Learning and Intelligent Systems
%          M S Ramaiah University of Applies Sciences
%          Bangalore  
%
% Ref.   : None
%
% Ver.   : 1.0
%          Base version
%_________________________________________________________________________

% Some initial 'arrangments'
clc;
close all;
clear;

% What should be our 'space'?
x1 = -5:0.4:5;
x2 = x1;

% Lets prepare for surface plot
[X1, X2] = meshgrid (x1, x2);

% A sphere function
func_1 = (X1.*X1) + (X2.*X2);
figure;
surf (X1, X2, func_1);
grid on;

% Function to plot BOHACHEVSKY FUNCTION
func_2 = (X1.*X1) + (2*(X2.*X2)) - (0.3*cos(3*pi*X1)) - (0.4*cos(4*pi*X2)) + 0.7; 
figure;
surf (X1, X2, func_2);
grid on;


% Function to plot EASOM FUNCTION
func_3 =(-1*cos(X1)).*(cos(X2)).*(exp((-1*((X1-pi).^2))-((X2-pi).^2)));
figure;
surf (X1,X2, func_3);