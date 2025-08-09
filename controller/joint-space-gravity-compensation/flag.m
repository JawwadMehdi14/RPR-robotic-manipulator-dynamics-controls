function map = flag(m)
if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end
% f = [red; white; blue; black]
f = [1 0 0; 1 1 1; 0 0 1; 0 0 0];
% Generate m/4 vertically stacked copies of f with Kronecker product.
e = ones(ceil(m/4),1);
map = kron(e,f);
map = map(1:m,:);
