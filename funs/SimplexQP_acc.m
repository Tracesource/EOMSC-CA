%  min  x'*A*x - x'*b
%  s.t. x'1=1, x>=0
function [x, obj]=SimplexQP_acc(A, b, x0)


NIter = 500;
NStop = 20;

[n] = size(A,1);
if nargin < 3
    x = 1/n*ones(n,1);
else
    x = x0;
end;

x1 = x;
t = 1;
t1 = 0;
r = 0.5;  % r=1/mu;
%obj = zeros(NIter,1);
for iter = 1:NIter
    p = (t1-1)/t;
    s = x + p*(x-x1);
    x1 = x;
    g = 2*A*s - b;
    ob1 = x'*A*x - x'*b;
    for it = 1:NStop
        z = s - r*g;
        z = EProjSimplex_new(z,1);  % z'1=1;z>=0; z=alpha;
        ob = z'*A*z - z'*b;
        if ob1 < ob
            r = 0.5*r; % rho=2;
        else
            break;
        end;
    end;
    if it == NStop
        obj(iter) = ob;
        %disp('not');
        break;
    end;
    x = z;
    t1 = t;
    t = (1+sqrt(1+4*t^2))/2;
    
    
    obj(iter) = ob;
end
   
1;



function [x, ft] = EProjSimplex_new(v, k)

%
%% Problem
%
%  min  1/2 || x - v||^2
%  s.t. x>=0, 1'x=k
%

if nargin < 2
    k = 1;
end;

ft=1;
n = length(v);

v0 = v-mean(v) + k/n;
%vmax = max(v0);
vmin = min(v0);
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = v0 - lambda_m;
        posidx = v1>0;
        npos = sum(posidx);
        g = -npos;
        f = sum(v1(posidx)) - k;
        lambda_m = lambda_m - f/g;
        ft=ft+1;
        if ft > 100
            x = max(v1,0);
            break;
        end;
    end;
    x = max(v1,0);

else
    x = v0;
end;    
