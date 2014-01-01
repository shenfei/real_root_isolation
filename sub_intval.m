function Z = sub_intval(X,Y,n)

Z = intval(zeros(n,1));
for i = 1:n
    if in0(Y(i),X(i)) || in0(X(i),Y(i)) || Y(i) == X(i) || isnan(Y(i))
        Z(i) = X(i);
    elseif inf(X(i)) < inf(Y(i))
        Z(i) = infsup(inf(X(i)),inf(Y(i)));
    else
        Z(i) = infsup(sup(Y(i)),sup(X(i)));
    end
end
