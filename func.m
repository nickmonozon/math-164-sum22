% File name: func.m
% Usage: This file contains the function that is to be minimized using the steepest
% descent algorithm.
% Currently, the function is from Example 8.1 in Chong and Zak.

function y = func(x)
y = (x(1) - 4)^4 + (x(2) - 3)^2 + 4*(x(3) + 5)^4;