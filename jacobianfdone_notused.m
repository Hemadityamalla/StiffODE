function Jest = jacobianfdone(estPoint,stepsize,func,depvar)
%Example
% f = @(x)[x(1)^2 + x(2)^2; x(1)^3.*x(2)^3];
%Point at which to estimate it
% x = [1;1];
%Step to take on each dimension (has to be small enough for precision)
% h = 1e-5*ones(size(x));
% jacobianfdone(x,h,f);
%------------------------------------------%
% Jacobian functor
J = @(x,h,F)(F(depvar,repmat(x,size(x'))+diag(h))-F(depvar,repmat(x,size(x'))))./h';
% Compute the jacobian
Jest = J(estPoint,stepsize,func);
end