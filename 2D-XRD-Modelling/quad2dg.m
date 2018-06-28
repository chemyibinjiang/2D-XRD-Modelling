function [int,numpts] = quad2dg(fun,xlow,xhigh,ylow,yhigh,tol,varargin)
%usage:  int = quad2dg('Fun',xlow,xhigh,ylow,yhigh)
%or
%	 int = quad2dg('Fun',xlow,xhigh,ylow,yhigh,tol)
%	 [int,numpts] = quad2dg('Fun',xlow,xhigh,ylow,yhigh,tol)
%or
%        int = quad2dg('Fun',xlow,xhigh,ylow,yhigh,tol,p1,...)
%        [int,numpts] = quad2dg('Fun',xlow,xhigh,ylow,yhigh,tol,p1,...)
%
%This function is similar to QUAD or QUAD8 for 2-dimensional integration,
%but it uses a Gaussian quadrature integration scheme.
%	int	-- value of the integral
%	numpts	-- order of finally integration quadrature
%	Fun	-- Fun(x,y) (function to be integrated)
%	xlow	-- lower x limit of integration
%	xhigh	-- upper x limit of integration
%	ylow	-- lower y limit of integration
%	yhigh	-- upper y limit of integration
%	tol	-- tolerance parameter (optional)
%       p1,...  -- optional parameters for the function fun
%Note that if there are discontinuities the region of integration
%should be broken up into separate pieces.  And if there are singularities,
%a more appropriate integration quadrature should be used
%(such as the Gauss-Chebyshev for a specific type of singularity).

%This routine could be optimized.

if nargin<6
  tol=1e-3;
elseif tol==[],
  tol=1e-3;
end

n=length(xlow);
nquad=2*ones(n,1);
[bpx,bpy,wfxy] = grule2d(2,2);
int_old=gquad2d(fun,xlow,xhigh,ylow,yhigh,bpx,bpy,wfxy,varargin{:});

converge='n';
for i=1:7,
  [bpx,bpy,wfxy] = grule2d(2^(i+1),2^(i+1));
  int=gquad2d(fun,xlow,xhigh,ylow,yhigh,bpx,bpy,wfxy,varargin{:});

  if (abs(int_old-int) < abs(tol*int) | abs(int_old-int) < abs(tol)^2),
    converge='y';
    break;
  end
  int_old=int;
end



if nargout==2,
  numpts=2^(i+1);
end
