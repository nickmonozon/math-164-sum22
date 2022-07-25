% Implementation of the conjugate gradient algorithm from Chong and Zak,
% Introduction to Optimization (4th Edition).

% Steps to the algorithm can be found on page 183.

% File name: conjugate_gradient.m
% Usage: Script to find the minimizer of a function based on its gradient, 
% which can be found in grad.m.

%OPTIONS(1) controls how much display output is given; set
% to 1 for a tabular display of results, (default is no display: 0).

%OPTIONS(2) is a measure of the precision required for the final point.

%OPTIONS(3) is a measure of the precision required of the gradient.

%OPTIONS(5) specifies the formula for beta (0 is
% Powell, 1 is Fletcher-Reeves, 2 is Polak-Ribiere, 3 is
% Hestenes-Stiefel)

%OPTIONS(14) is the maximum number of iterations.

function [x,N] = conjugate_gradient(grad, xnew, options)

if nargin ~= 3
    options = [];
    if nargin ~= 2
        disp("Wrong number of arguments.");
        return;
    end
end

numvars = length(xnew);
if length(options) >= 14
    if options(14)==0
        options(14)=1000*numvars;
    end
else
    options(14)=1000*numvars;
end

clc;
format compact;
format short e;

options = foptions(options);
print = options(1);
epsilon_x = options(2);
epsilon_g = options(3);
max_iter=options(14);

g_curr=feval(grad,xnew);
if norm(g_curr) <= epsilon_g
    disp("Terminating: Norm of initial gradient less than");
    disp(epsilon_g);
    return;
end 

d=-g_curr;
reset_cnt = 0;
for k = 1:max_iter

    xcurr=xnew;
    alpha=linesearch_secant(grad,xcurr,d);
    % alpha=-(d *g_curr)/(d *Q*d);
    xnew = xcurr+alpha*d;

    if print
        disp("Iteration number k =")
        % Index of iteration (k)
        disp(k); 
        disp("alpha =");
        % Alpha
        disp(alpha); 
        disp("Gradient = ");
        % Gradient
        disp(g_curr); 
        disp("New point =");
        % New point
        disp(xnew); 
    end 

    if norm(xnew-xcurr) <= epsilon_x*norm(xcurr)
        disp("Terminating: Norm of difference between iterates less than");
        disp(epsilon_x);
        break;
    end 

    g_old=g_curr;
    g_curr=feval(grad,xnew);

    if norm(g_curr) <= epsilon_g
        disp("Terminating: Norm of gradient less than");
        disp(epsilon_g);
        break;
    end 

    reset_cnt = reset_cnt+1;
    if reset_cnt == 3*numvars
        d=-g_curr;
        reset_cnt = 0;
    else
        %Powell
        if options(5)==0 
            beta = max(0,(g_curr*(g_curr-g_old))/(g_old*g_old));
        % Fletcher-Reeves
        elseif options(5)==1
            beta = (g_curr*g_curr)/(g_old*g_old);
        % Polak-Ribiere
        elseif options(5)==2 
            beta = (g_curr*(g_curr-g_old))/(g_old*g_old);
        % Hestenes-Stiefel
        else 
            beta = (g_curr*(g_curr-g_old))/(d*(g_curr-g_old));
        end 
        d=-g_curr+beta*d;
    end

    if print
        disp("New beta =");
        disp(beta);
        disp("New d =");
        disp(d);
    end

    if k == max_iter
        disp("Terminating with maximum number of iterations");
    end
end 

if nargout >= 1
    x=xnew;
    if nargout == 2
        N=k;
    end
else
    disp("Final point =");
    disp(xnew);
    disp("Number of iterations =");
    disp(k);
end