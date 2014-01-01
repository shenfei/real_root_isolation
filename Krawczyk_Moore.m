%function [real_roots,real_roots_num,flag] = Krawczyk_Moore(F,X0,n,real_roots,real_roots_num,k)
function flag = Krawczyk_Moore(F,X0,n,k)

global real_roots real_roots_num iter_num

gx = gradientinit(X0);
gy = F(gx);
Y = mid(gy.dx);
% if det(Y) < 1e-10
%     Y = rand(n);
%     while det(Y) < 1e-10
%         Y = rand(n);
%     end
% end

Xs = Krawczyk(F,X0,Y,gy,n);
%Xs = Krawczyk(F,X0,Y);
iter_num = iter_num + 1;

if all(in0(Xs,X0))
    flag = true;
    real_roots_num = real_roots_num + 1;
    real_roots{real_roots_num,1} = Xs;
    real_roots{real_roots_num,2} = k;
    return
end

Xt = intersect(Xs,X0);
%[empty,Xt] = emptyintersect(Xs,X0);
if any(isnan(Xt))
%if any(empty)
    flag = false;
    return
end

%if max(rad(Xt)) > 0.5 * 1e-14
    if int_singular(gy.dx)
        %[real_roots,real_roots_num,flag] = special_judge(F,Xt,n,real_roots,real_roots_num,k);
        flag = special_judge(F,Xt,n,k);
        return
    end
%end

Xt = Xs;
Xs = X0;
while ~all(in0(Xt,Xs))
    iter_num = iter_num + 1;
    X = intersect(Xt,Xs);
    %[empty,X] = emptyintersect(Xt,Xs);
%     if X == Xs
%         if max(rad(X)) > 0.5 * 1e-14
%             [real_roots,real_roots_num,flag] = special_judge(F,X,n,real_roots,real_roots_num,k);
%             return
%         else
%             Xt = X;
%             break;
%         end
%     end
    if X == Xs
        %[real_roots,real_roots_num,flag] = special_judge(F,X,n,real_roots,real_roots_num,k);
        flag = special_judge(F,Xt,n,k);
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
    
%     if norm(eye(n) - Y \ gy.dx,inf) < 1
%         break
%     end
    
    Xs = X;
    Xt = Krawczyk(F,Xs,Y,gy,n);
    %[empty,interY] = emptyintersect(Xt,Xs);
    %if any(empty)
    if any(isnan(intersect(Xt,Xs)))
        flag = false;
        return
    end
end

flag = true;
real_roots_num = real_roots_num + 1;
real_roots{real_roots_num,1} = Xt;
real_roots{real_roots_num,2} = k;
