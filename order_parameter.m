clear all
close all
clc
%% polar order parameter example
% an example of polar order parameter calculation, the three feet are 
% initially completey out of phase but reach the same phase at the end

theta = [0 120  240;      % a matrix containing the angles of each foot (column) at every time step (row)    
         23 110  200;
         49 104 175;
         67 101 1143;
         100 100 100];
     
N = [3, 3, 3, 3, 3]';     % a vector of number of attached feet at each time
     
theta = (pi/180)*theta;   % converts angles into radians (needed for measurments in degrees only)

% looping over the length of our data vector (time) and computing order
% parameter at every time step
for i =1:length(N)
    p(i) = polar_order(theta(i,:),N(i));
end


p

%%
function [ OP ] = polar_order(phase,num)
% takes in phases of the oscillators (an array) and the number of oscillators (scalar) 
% computes order parameter of the oscillators (scalar)

z = 0;  % initialize order parameter to zero 

for j = 1:num
    z = z + exp(1i*phase(j)); 
end

z = z/num;   % divide by number of oscillators
OP = abs(z); % compute absolute value which is the polar order paremter
end



