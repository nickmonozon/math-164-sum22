%% Line search using the Secant Method

% Information and set up

% INPUTS:
%   'grad'      =   .m file with gradient
%   x           =   starting line search point
%   d           =   search direction
%   epsilon     =   tolerance, e.g. 1e-4 = 10^{-4}
%   max         =   maximum number of iterations

% OUTPUTS:
%   alpha       =   value returned by function

alpha_curr = 0;
alpha = 0.001;
dphi_zero = feval(grad,x)'*d;
dphi_curr=dphi_zero;
epsilon = 1e-4;
max = 100;

i = 0;

while abs(dphi_curr) > epsilon*abs(dphi_zero),
    alpha_o = alpha_curr;
    alpha_curr = alpha;
    dphi_o = dphi_curr;
    dphi_curr = feval(grad, x+alpha_curr*d)'*d;
    alpha = (dphi_curr*alpha_o-dphi_o*alpha_curr)/(dphi_curr-dphi_o);
    i = i + 1;
    if (i >= max) & (abs(dphi_curr)>epsilon*abs(dphi_zero)),
        fprintf('\nLine search algorithm converges after %d iterations.\n\n',i);
        break;
    end
end
