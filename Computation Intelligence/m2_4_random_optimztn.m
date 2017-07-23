% This m-file finds the minima of the given function using random search
% based algorithm
%
% Author            : Aravind D. Chakravarti
%                     PT-2017
%                     MLIS
%                     M S Ramaiah University of Applied Sciences, Bengaluru
%                    
% Version           : 1.0
%                     Base version

% Clear everything on command window
clc;

% Clear variables used
clear;

% close if any windows is open!
close all;

% Our search space is limited to
range_min = -5;
range_max =  5;

% Dimension of our world
dimension =  2;

% Let us initialiaze few things which are helpful
minimum_postn = [5 5];
minimum_fitns = 50; % 5^2 + 5^2 = 50

% Lets start searching
for index = 1:1000
    
    % Generate a random number in the range range_min to range_max
    rand_num = range_min + ((range_max - range_min)* (rand(1,dimension)));
    
    % Given fitness function is sphere function
    fitnes_of_fun = (rand_num(1).^2)+(rand_num(2).^2);
    
    % Is this function has minimum than previous one?
    if fitnes_of_fun < minimum_fitns
        minimum_fitns = fitnes_of_fun;
    end
    
end

% Display the result
disp (minimum_fitns);