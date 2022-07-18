% File name: grad.m
% Usage: This file contains the gradient of the function that is to be minimized using the steepest
% descent algorithm.
% Currently, the gradient corresponds to the function from Example 8.1 in Chong and Zak.

function y = grad(x)
y(1) = 4*(x(1) - 4)^3;
y(2) = 2*(x(2) - 3);
y(3) = 16*(x(3) + 5)^3;
y = y';