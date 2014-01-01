%function [real_roots,real_roots_num,flag] = special_judge(F,X0,n,real_roots,real_roots_num,k)
function flag = special_judge(F,X0,n,k)
global real_roots real_roots_num

if max(rad(X0)) <= 0.5 * 1e-12
    [Xt,flag] = zero_check(F,X0,n);
    if flag
        real_roots_num = real_roots_num + 1;
        real_roots{real_roots_num,1} = Xt;
        real_roots{real_roots_num,2} = k;
    end
    return
end

y = F(slopeinit(mid(X0),X0));
if ~all(in0(0,y.r))
    %Xt = X0;
    flag = false;
    return
end

%if rad(y.r) < 1e-20
%    Xt = X0;
%    flag = true;
%    return
%end

gx = gradientinit(X0);
gy = F(gx);
Jac = gy.dx;
% find the plane to cut
singular_flag = any(inf(abs(Jac)));
if ~all(singular_flag)
    idx = find(singular_flag == 0);
else
    [m,idx] = max(rad(X0));
end
% if all(sigular_flag)
%    idx = 1;
%    temp = min(inf(abs(Jac(:,1))));
%    for j = 2:length(X0)
%        if min(inf(abs(Jac(:,j)))) < temp
%            idx = j;
%            temp = min(inf(abs(Jac(:,j))));
%        end
%    end
% else
%     idx = find(sigular_flag == 0);
% end
% 
% if rad(X0(idx)) < 1e-10
%     [m,idx] = max(rad(X0));
% end


X1 = X0;
X2 = X0;
X1(idx) = infsup(inf(X0(idx)),mid(X0(idx)));
X2(idx) = infsup(mid(X0(idx)),sup(X0(idx)));
%X1 = infsup(inf(X0),mid(X0));
%X2 = infsup(mid(X0),sup(X0));
% [real_roots,real_roots_num,s1] = Krawczyk_Moore(F,X1,n,real_roots,real_roots_num,k);
% [real_roots,real_roots_num,s2] = Krawczyk_Moore(F,X2,n,real_roots,real_roots_num,k);
s1 = Krawczyk_Moore(F,X1,n,k);
s2 = Krawczyk_Moore(F,X2,n,k);

if s1 || s2
    flag = true;
else
    flag = false;
end
