function r = init_width5(F,x,n)

%residue = norm(F(x),inf);
residue = norm(F(x));
if residue > 0.5
    r = -1;
    return
end

gx = gradientinit(x);
gy = F(gx);
J = gy.dx;
delta = gy.dx \ F(x);

hx = hessianinit(x);
hy = F(hx);
C = -1;
for j = 1:n
    ct = sum(max(abs(hy(j).hx)));
    if ct > C
        C = ct;
    end
end

%invJ = norm(inv(J),inf);
invJ = norm(inv(J));
%r = invJ * residue / (1 - sqrt(n) * C * invJ * norm(delta,2));
r = sqrt(n) * invJ * residue / (1 - n * C * invJ * norm(delta));
