function x=indeks(y,in)
% x=INDEKS(y,in)
%
% INPUT:
%
% y         Some vector
% in        Some set of indices
%
% Extracts indexed positions out of simple matrices
%
% indeks([1 2],2) 
% indeks([1 2],':,2n')
% indeks([1 2],'end')
%
% Works for logical and numeric indices.
%
% Last modified by fjsimons-at-alum.mit.edu, 04/09/2007


defval('in',1)

if ~isstr(in)
  x=y(in);
else
  eval([ 'x=y(' in ');'])
end
