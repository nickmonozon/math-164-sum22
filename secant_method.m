%% Secant Method

% This code approximates a solution to f(x) = 0 given initial approximations p0 and p1.

% INPUTS:
%   p0          =   initial approximation
%   p1          =   another initial approximation
%   epsilon     =   tolerance

% OUTPUTS:
%   p           = 	approximate solution

%% Information and set up

f = @(x) (2*x-1)^2 + 4*(4-1024*x)^2;        % function whose root we want to approximate
p0 = 0;                                     % first initial approximation
p1 = 1;                                     % second initial approximation
tol = 1e-5;                                 % tolerance, e.g. 1e-4 = 10^{-4}


%% Secant Method

i = 2;                                    % iteration count

q0 = f(p0);
q1 = f(p1);

while true

    % get p_i 
    p = p1 - q1*(p1 - p0)/(q1 - q0);      
    
    % check stopping condition 
    if(abs(p - p1) < tol)
        break;
    end
    
    % increase iteration count
    i = i + 1;
    
    % prepare for next iteration
    p0 = p1;
    q0 = q1;
    p1 = p;
    q1 = f(p);
    
end

%% Display Information
 
fprintf('\nSecant Method approximated the solution p = %.10f after %d iterations.\n\n',p,i);