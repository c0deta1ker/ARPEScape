function dispname(varargin)
% dispname(x [,...]) displays variables in coma-separated argument list with names
for i=1:nargin
  disp([inputname(i) '= ' num2str(varargin{i})])
end
% function dispname(name)
% % dispname(X) -- Display scalar variable(s) in argument list with name(s)
% disp([inputname(1) '= ' num2str(name)])