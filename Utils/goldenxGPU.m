% GOLDENX Computes local maximum of univariate function on interval via Golden Search
%   Vectorized version of Golden
% USAGE
%   [x,fval] = goldenx(f,a,b,P1,P2,...);
% INPUTS
%   f         : name of function of form fval=f(x)
%   a,b       : left, right endpoints of interval
%   P1,P2,... : optional additional arguments for f
% OUTPUTS
%   x       : local maximum of f
%   fval    : function value estimate
%
% USER OPTIONS (SET WITH OPSET)
%   tol     : convergence tolerance

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

% This function has been modified to
%   1. Check that bounds are actually correct, and return errors if not
%   2. Adjust for NaN's

%% Indexing is very slow on GPU modify calculation of x1 and f1
function [x1,f1] = goldenxGPU(f,a,b,varargin)

%%% make sure 
tol = optget('goldenx','tol',sqrt(eps));

alpha1 = (3-sqrt(5))/2;
alpha2 = 1-alpha1;
d  = b-a;
x1 = a+alpha1*d;
x2 = a+alpha2*d;
s  = ones(size(x1));
f1 = feval(f,x1,varargin{:});
f2 = feval(f,x2,varargin{:});
d = alpha1*alpha2*d;
while any(any(d>tol))
  f2(isnan(f2)) = -1e22;
  f1(isnan(f1)) = -1e22;
  i = f2>f1;
  x1 = x2.*i+x1.*(1-i);
% f1(i) = f2(i);
  f1 = f2.*i+f1.*(1-i);
  d = d*alpha2;
  x2 = x1+s.*(i-(~i)).*d;
  s = sign(x2-x1);
  f2 = feval(f,x2,varargin{:});
end

% Return the larger of the two
i = f2>f1; x1(i) = x2(i);  f1(i) = f2(i);
