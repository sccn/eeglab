% MYFPVAL - parametric F tail used by statcond and teststat.
function p = myfpval(F, df1, df2)

p = ones(size(F), class(F));
p(isnan(F)) = NaN;
use = (F > 0);
if ~any(use(:))
    return
end
Fv = F(use);
lo = Fv < (df2 ./ df1);
pv = zeros(size(Fv), class(F));
if any(~lo(:))
    Ft = Fv(~lo);
    pv(~lo) = betainc(df2 ./ (df2 + df1 .* Ft), df2/2, df1/2);
end
if any(lo(:))
    Ft = Fv(lo);
    pv(lo) = betainc((df1 .* Ft) ./ (df2 + df1 .* Ft), df1/2, df2/2, 'upper');
end
p(use) = pv;
