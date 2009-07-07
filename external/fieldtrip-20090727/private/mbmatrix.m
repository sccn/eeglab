function bool = mbmatrix(x)

bool = (ndims(x) == 2) & (prod(size(x)) > 1);
