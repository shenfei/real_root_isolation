%function [real_roots,real_roots_num] = narrowing(real_roots,real_roots_num,F,n,tau)
function  narrowing(F,n,tau)
global real_roots real_roots_num

for i = 1:real_roots_num
    X0 = real_roots{i,1};
    %while norm(rad(X0),inf) >= tau
    %不能使用无穷范数，否则可能出现最窄的区间变成nan的情况
    while min(rad(X0)) >= tau
        gx = gradientinit(X0);
        gy = F(gx);
        Y = mid(gy.dx);
        Xs = Krawczyk(F,X0,Y,gy,n);
        Xt = intersect(Xs,X0);
%         infsup(X0)
%         infsup(Xs)
%         infsup(Xt)
        
        if Xt == X0
            [M,idx] = max(rad(X0));
            X1 = X0;
            X2 = X0;
            X1(idx) = infsup(inf(X0(idx)),mid(X0(idx)));
            X2(idx) = infsup(mid(X0(idx)),sup(X0(idx)));
            [Xt,flag] = zero_check(F,X1,n);
            if ~flag
                [Xt,flag] = zero_check(F,X2,n);
            end
        end
        X0 = Xt;
    end
    real_roots{i,1} = X0;
end