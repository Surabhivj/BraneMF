function [alam,x,fval,exitflag,output] = linemin(fun, x0, d,options,N,M,L,D,Q,alpha,beta)
% LINEMIN Line minimizer
%
% [alam,x,fval,exitflag,output] = linemin(fun, x0, d,options)
%
% alam is the scalar (less than or equal to 1) amount stepped in
% direction d for x to get the resulting value fval. exitflag
% indicates the termination condition, and output is a struct with
% information about the optimization.
%
% Valid fields for options struct include:
%
% MaxIter	Maximum number of optimization iterations
% TolFun	Tolerance on the function value
% Display	off: 	No display information
%			iter: 	Display information at every iteration [default]
%
% Values of the exitflag indicate
%
%  -2     Internal error
%  -1     Slope is positive
%   0     Iteration count exceed
%   1     Sufficient decrease (normal termination)
%   2     Jump required is too small (alam returned is 0)

% Written by Jerod Weinman (jerod@acm.org)
%
% Ref:
% See Numerical Recipes in C, p.385. lnsrch. A simple backtracking line
% search. No attempt at accurately finding the true minimum is
% made. The goal is only to ensure that linemin will return a
% position of lower value.

  if nargin<4
	options = struct;
  end;
  
  if isfield(options,'MaxIter') & ~isempty(options.MaxIter)
	maxIter = options.MaxIter;
  else
	maxIter = 100;
  end;
  
  if isfield(options,'MaxStep') & ~isempty(options.MaxStep)
	stpmax = options.MaxStep;
  else
	stpmax = 100;
  end;
  
  if isfield(options,'TolFun') & ~isempty(options.TolFun)
	tolx = options.TolFun;
  else
	tolx = 1e-10;
  end;

  if isfield(options,'Display') & ~isempty(options.Display)
	switch options.Display
	 case 'off'
	  dspOpt = 0;
	 case 'iter'
	  dspOpt = 1;
	 otherwise
	  error('Unknown value for Display option');
	end;
  else
	dspOpt = 1;
  end;
  
  
  alf = 1e-4; % Ensures sufficient decrease in function value.
  
  stpmod = 1;
  
  alam2 = 0;
  tmplam = 0;

  x = x0;

  
  [fval,g] = feval(fun,x , N,M,L,D,Q,alpha,beta);
  
  fold = fval;
  f2 = fold;

  funcCount = 1;
  
  if (norm(d,2)>stpmax)
	if dspOpt
	  display('linemin: step too big: scaling direction');
	end;
	
	stpmod = stpmax/norm(d,2);
	d = d*stpmod;
  end;
  
  slope = g*d';
  

  if (slope>0)
	%%%%% Bad Slope 
	alam = 0;
	x = x0;
	fval = fold;
	msg = 'Line min: Slope is positive';
	exitflag = -1;

	output = struct('iterations',{0},...
					'funcCount',{funcCount},...
					'algorithm',{'backtrack line minimization'},...
					'message',{msg});
	if dspOpt
	  display(msg);
	end;
	return;
	%%%%%
  end;

  alamin = tolx/max(abs(d)./max(abs(x),1));
  %alam = 1;
  alam = 0.5;
  oldalam = 0;
  
  for k=1:maxIter
	
	%x = x + d.*(alam-oldalam);
    x = x + d;

	
	oldalam = alam;
	
	if (alam<alamin)
	  %%%% Small jump
	  alam = 0;
	  x = x0;
	  fval = fold;
	  msg = 'Line min: Jump too small';
	  exitflag = 2;
	  
	  output = struct('iterations',{k},...
					  'funcCount',{funcCount},...
					  'algorithm',{'backtrack line minimization'},...
					  'message',{msg});
	  if dspOpt
		display(msg);
	  end;
	  return;
	  %%%%%%
	end;

	funcCount= funcCount + 1;
	fval = feval(fun,x,N,M,L,D,Q,alpha,beta);

	if (fval <= fold + alf*alam*slope) 
	  %%%% Wolfe condition

	  alam = alam*stpmod;
	  
	  msg = 'Line min terminated: sufficient decrease';
	  exitflag = 1;
	  
	  output = struct('iterations',{k},...
					'funcCount',{funcCount},...
					'algorithm',{'backtrack line minimization'},...
					'message',{msg});
	  
	  if dspOpt
		display(msg);
	  end;
	  
	  return;
	  %%%%%
	elseif isinf(fval) | isinf(f2)
	  
	  if dspOpt
		display('Line min: Inf value. Scaling step size.');
	  end;
	  
	  tmplam = .1* alam;
	  
	else % Backtrack
	  if alam==1 % first time through
		tmplam = -slope./(2*(fval-fold-slope));
	  else
		rhs1 = fval-fold-alam*slope;
		rhs2 = f2-fold-alam2*slope;
		
		if (alam==alam2)
		  msg = 'Failure: dividing by alam-alam2=0';
		  exitflag = -2;
		  
		  output = struct('iterations',{k},...
						  'funcCount',{funcCount},...
						  'algorithm',{'backtrack line minimization'},...
						  'message',{msg});
		  
		  if dspOpt
			display(msg);
		  end;
		  return;
		end;
		
		a = (rhs1./(alam.*alam)-rhs2./(alam2.*alam2))./(alam-alam2);
		b = (-alam2.*rhs1./(alam.*alam)+alam.*rhs2./(alam2.*alam2))./(alam-alam2);

		if (a == 0.0) 
		  tmplam = -slope./(2.0*b);
		else
		  disc = b.*b-3.0.*a.*slope;
		  if (disc < 0.0)
			tmplam = .5 * alam;
		  elseif (b <= 0.0)
			tmplam = (-b+sqrt(disc))/(3.0*a);
		  else 
			tmplam = -slope./(b+sqrt(disc));
		  end;
		end;

		if (tmplam > .5*alam)
		  tmplam = .5*alam;    % lambda <= .5 lambda_1
		end;
	  end;
	end;
	
	alam2 = alam;
	f2 = fval;

	if dspOpt
	  display(sprintf('linemin\t%d\t%g\t%g\t%g',k,fval,alam,alamin));
	end;
	
	alam = max(tmplam, .1*alam);
	

  end;
  
  
  alam = 0;
  x = x0;
  fval = fold;
  msg = 'Line min: Too many iterations';
  exitflag = 0;
  output = struct('iterations',{k},...
				  'funcCount',{funcCount},...
				  'algorithm',{'backtrack line minimization'},...
				  'message',{msg});
  if dspOpt
	display(msg);
  end;
