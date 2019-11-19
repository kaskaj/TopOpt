function err = Diff_Derivatives(f, df_x, x, dir)

fun   = f(x);
der1  = df_x'*dir;
s_all = 10.^(linspace(-10,-5,6));

err = Inf;
for i=1:length(s_all)
    
    s    = s_all(i);
    der2 = (f(x+s*dir) - fun)/s;
    err  = min(err, norm(der1-der2)/max(norm(der1), norm(der2)));
end

end