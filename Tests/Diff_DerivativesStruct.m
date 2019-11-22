function err = Diff_DerivativesStruct(f, df_x, x, param)

fun   = f(x);
der1  = df_x;
s_all = 10.^(linspace(-10,-5,6));

err = Inf;
for i=1:length(s_all)
    
    s    = s_all(i);
    x.B_mu.a_w(param) = x.B_mu.a_w(param)+s;
    der2 = (f(x) - fun)/s;
    err  = min(err, norm(der1-der2)/max(norm(der1), norm(der2)));
end

end