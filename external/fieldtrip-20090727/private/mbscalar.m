function mbscalar(a);

if ~all(size(a)==1)
  error('Argument to mbscalar must be scalar');
end
