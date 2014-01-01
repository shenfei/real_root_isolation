function [Xr,flag] = zero_check(F,X0,n)

%判断X0区域内是否有根，且假定初始区域最多只有一个零点
Xr = X0;
y = F(slopeinit(mid(X0),X0));
if ~all(in0(0,y.r))
    flag = false;
    return
end

gx = gradientinit(X0);
gy = F(gx);
Y = mid(gy.dx);
Xs = Krawczyk(F,X0,Y,gy,n);

if all(in0(Xs,X0))
    flag = true;
    Xr = Xs;
    return
end

Xt = intersect(Xs,X0);
if any(isnan(Xt))
    flag = false;
    Xr = X0;
    return
end

Xt = Xs;
Xs = X0;
while ~all(in0(Xt,Xs))
    X = intersect(Xt,Xs);
    if X == Xs
        Xr = X;
        flag = true;
        return
    end
    
    Y0 = Y;
    gx = gradientinit(X);
    gy = F(gx);
    gx0 = gradientinit(Xs);
    gy0 = F(gx0);
    Y = mid(gy.dx);
    if norm(eye(n) - Y \ gy.dx) > norm(eye(n) - Y0 \ gy0.dx)
        Y = Y0;
    end
    Xs = X;
    Xt = Krawczyk(F,Xs,Y,gy,n);
    if any(isnan(intersect(Xt,Xs)))
        flag = false;
        return
    end
end

flag = true;
%Xr = Xt;
Xr = X0;

