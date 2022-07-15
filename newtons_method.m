%% Newton's Method

% This code approximates a solution to f(x) = 0 given an initial
% approximation p0.

% INPUTS:
%   p0          = initial approximation
%   epsilon     = tolerance
%   max_iter    = maximum number of iterations

% OUTPUTS:
%   p           = 	approximate solution
%   or error message

%% Information and set up

f = @(x) exp(x) - 2;        % function whose root we want to approximate
fprime = @(x) exp(x);       % derivative of function whose root we want to approximate

p0 = 1;

tol = 1e-4;                 % tolerance, 1e-4 = 10^{-4}

max_iter = 30;              % max number of iterations

%% Newton's Method

i = 1;                                        % iteration count

fprintf('i\tp_i\t\tf(p_i)\n');                % for display
fprintf('%d\t%.10f\t%.10f\n',0,p0,f(p0));     % displays iteration i, p_i, f(p_i)

while( i <= max_iter)

    % get p_i 
    p = p0 - f(p0)/fprime(p0);    % p_i 
    
    % display information
    fprintf('%d\t%.10f\t%.10f\n',i,p,f(p));   % displays iteration i, p_i, f(p_i)
    
    % check stopping condition 
    if(abs(p - p0) < epsilon)
        break;
    end
    
    % increase iteration count
    i = i + 1;
    
    % prepare for next iteration
    p0 = p;
    
end

%% Display Information

if( i <= max_iter )         % successful
    fprintf('\nNewton''s Method approximated the solution p = %.10f after %d iterations.\n\n',p,i);
else                        % not successful 
    fprintf('\nNewton''s Method did not converge within the tolerance in %d iterations.\n\n',max_iter)
end
