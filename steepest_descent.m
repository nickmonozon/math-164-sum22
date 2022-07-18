% Implementation of the steepest descent algorithm from Chong and Zak,
% Introduction to Optimization (4th Edition)

% File name:  steepest_descent.m
% Usage: This MATLAB script implements the steepest descent algorithm.
% Stopping condition: tolerance (epsilon) supplied by user
% Note: The function value is read from the file "func.m".
% Note: The gradient value is read from the file "grad.m".

% Prompt user to supply inputs
  n = input('How many variables does the function have? ');
  x = input('Give an initial condition: ');
  epsilon = input('Supply a value of epsilon: ');

% Stepsize rules 
  sigma = .1;
  beta = .5;
  obj = func(x);
  g = grad(x);
  k = 0;                    % k = iteration counter
  nf = 1;					% nf = number of function evals

% Begin algorithm
  while norm(g) > epsilon   
    d = -g;                 % steepest descent direction (negative gradient)
    a = 1;
    newobj = func(x + a*d);
    nf = nf+1;
    while (newobj-obj)/a > sigma*g'*d
      a = a*beta;
      newobj = func(x + a*d);
      nf = nf+1;
    end
    x = x + a*d;
    obj = newobj;
    g = grad(x);
    k = k + 1;
  end

% Output number of iterations and approximation
fprintf('\nThe steepest descent algorithm converges within tolerance after %d iterations. The approximation is: \n', k);
x
