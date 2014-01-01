function [r1,r2] = init_width7(F,x,n)

% residue = norm(F(x));
% if residue > 0.5
%     r = -1;
%     return
% end

gx = gradientinit(x);
gy = F(gx);
J = gy.dx;
eta = norm(gy.dx \ F(x),inf);
B = norm(inv(J),inf);
w = 2 * eta;
X0 = midrad(x,w);

% hx = hessianinit(X0);
% hy = F(hx);
% K = -1;
% for j = 1:n
%     ct = sum(max(sup(abs(hy(j).hx))));
%     if ct > K
%         K = ct;
%     end
% end
gK = F(gradientinit(X0));
K = norm(sup(abs(gK.dx)),inf);

h = n * K * B * eta;

while h > 1/2
    x = x - gy.dx \ F(x);
    gx = gradientinit(x);
    gy = F(gx);
    J = gy.dx;
    eta = norm(gy.dx \ F(x),inf);
    B = norm(inv(J),inf);
    w = 2 * eta;
    X0 = midrad(x,w);

%     hx = hessianinit(X0);
%     hy = F(hx);
%     K = -1;
%     for j = 1:n
%         ct = sum(max(sup(abs(hy(j).hx))));
%         if ct > K
%             K = ct;
%         end
%     end
gK = F(gradientinit(X0));
K = norm(sup(abs(gK.dx)),inf);
    h = n * K * B * eta;
end

% B
% K
% eta
% h
r1 = (1 - sqrt(1-2*h)) / h * eta;
r2 = (1 + sqrt(1-2*h)) / h * eta;
