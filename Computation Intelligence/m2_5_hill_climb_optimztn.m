% This m-file finds the minima of the given function using 'hill-climbing'
% algorithm. 
%
% Though we are decending a hill here, because mathematically
% hill climbing and decending is same except signs, we can still keep the
% same name :)
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

% I want to pop up figures as soon as run the program.
close all;

% Our search space is limited to
range_min =   -5;
range_max =    5;
step_size = 0.05;

% Dimension of our world
dimension =  2;

% Let us initialiaze few things which are helpful
minimum_fitns = 50; % 5^2 + 5^2 = 50

% Generate a random number in the range range_min to range_max
minima_coardntes = range_min + ((range_max - range_min)* (rand(1,dimension)));

% Lets start climbing the hill
for index = 1:1000
    
    %                    o      o
    %                     .    . 
    %                      .  .
    %                       P
    %                      . .
    %                     .   .
    %                    o     o
    %
    %
    
    co_ordinate_1 = [minima_coardntes(1)+step_size minima_coardntes(2)];
    co_ordinate_2 = [minima_coardntes(1)-step_size minima_coardntes(2)];
    
    co_ordinate_3 = [minima_coardntes(1) minima_coardntes(2)+step_size];
    co_ordinate_4 = [minima_coardntes(1) minima_coardntes(2)-step_size];
    
    % What is depth at positions near me?
    fitnes_1 = sphere_fun(co_ordinate_1);
    fitnes_2 = sphere_fun(co_ordinate_2);
    fitnes_3 = sphere_fun(co_ordinate_3);
    fitnes_4 = sphere_fun(co_ordinate_4);
    
    % Am I standing below the earlier one?
    if fitnes_1 < minimum_fitns
        minimum_fitns = fitnes_1;
        minima_coardntes = co_ordinate_1;
    end
    
    % Very cruid implementation. Do you know anything better?
    if fitnes_2 < minimum_fitns
        minimum_fitns = fitnes_2;
        minima_coardntes = co_ordinate_2;
    end
    
    if fitnes_3 < minimum_fitns
        minimum_fitns = fitnes_3;
        minima_coardntes = co_ordinate_3;
    end
    
    if fitnes_4 < minimum_fitns
        minimum_fitns = fitnes_4;
        minima_coardntes = co_ordinate_4;
    end

end

% Display the result
fprintf('You wanted to reach the hill\nUnfortunately you ended up in valley\n');
fprintf('Your current position is %d\t and co-ordinates are %d,%d\n',minimum_fitns, minima_coardntes(1), minima_coardntes(2));