%% Secant Method

% This code approximates a solution to f(x) = 0 given initial approximations p0 and p1.

% INPUTS:
%   pn1         =   initial approximation
%   p1          =   another initial approximation
%   epsilon     =   tolerance

% OUTPUTS:
%   p           = 	approximate solution

%% Information and set up

f = @(x) (2*x-1)^2 + 4*(4-1024*x)^4;        % function whose root we want to approximate
pn1 = 0;                                    % first initial approximation
p0 = 1;                                     % second initial approximation
epsilon = 1e-5;                             % tolerance, e.g. 1e-4 = 10^{-4}


%% Secant Method

qn1 = f(pn1);
q0 = f(p0);

while true

    % get p_i 
    p = p0 - ((p0 - pn1) * q0)/(q0 - qn1);
    
    % check stopping condition 
    if(abs(p - p0) < abs(p0)*epsilon)
        break;
    end
    
    
    % prepare for next iteration
    pn1 = p0;
    qn1 = q0;
    p0 = p;
    q0 = f(p);
    
end

%% Display Information
 
fprintf('\nSecant Method approximated the solution p = %.10f.\n\n',p);
