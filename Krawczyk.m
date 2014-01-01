function Xt = Krawczyk(F,X,Y,gy,n)

%gy is the interval Jacobian matrix

%Y is the inverse of ordinary nonsigular Y
% gx = gradientinit(X);
% gy = F(gx);
% 
% Xt = mid(X) - Y \ F(mid(X));
% Xt = Xt + abs(Y) \ rad(gy.dx) * rad(X)*(1/4)*midrad(0,1);

%Y = mid(gy.dx);
Xt = mid(X) - Y \ F(mid(X));
Xt = Xt + (eye(n) - Y \ gy.dx) * rad(X) * midrad(0,1);
