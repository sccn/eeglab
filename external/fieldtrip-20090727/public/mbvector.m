function mbvector(a)

if ndims(a) > 2 | (size(a, 1) > 1 & size(a, 2) > 1)
  error('Argument to mbvector must be a vector');
end
