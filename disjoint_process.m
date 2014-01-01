function [real_roots,real_roots_num] = disjoint_process(real_roots,real_roots_num,F,n)

k = 0;      %用来标记已验证的实根的index
for i = 1:real_roots_num
    X1 = real_roots{i,1};
    if any(isnan(X1))
        continue;
    end
    new_root = true;       %用于标记是否为新的根
    %初始假设为新的零点，后面来判断是否有重复
    for j = 1:k
        X2 = real_roots{j,1};
        intersect_flag = false;     %只用于标记部分相交的情况
        disjoint_flag = false;
        Y = intval(zeros(n,1));
        for d = 1:n
            if in0(X1(d),X2(d))
                X2(d) = X1(d);
                Y(d) = X1(d);
            elseif in0(X2(d),X1(d))
                X1(d) = X2(d);
                Y(d) = X2(d);
            else
                Y(d) = intersect(X1(d),X2(d));
                if isnan(Y(d))
                    disjoint_flag = true;
                else
                    intersect_flag = true;
                end
            end
        end
        if disjoint_flag
            continue;       %如果有的维度区间不交，则直接判断下面的
        end
        if intersect_flag
            [Yt,root_flag] = zero_check(F,Y,n);
            if root_flag
                real_roots{j,1} = Yt;
                new_root = false;
                break;
            else
                real_roots{j,1} = sub_intval(X2,Yt,n);
                X1 = sub_intval(X1,Yt,n);
                % 待优化
            end
        else
            %这里对应每个维度的区间都是相互包含的情况
            real_roots{j,1} = X2;   %narrowing
            new_root = false;
            break;          %说明重复，跳出内层for j循环，判断下一个点
        end
    end
    if new_root
        k = k + 1;
        real_roots{k,1} = X1;
    end
end
real_roots_num = k;

% for i = 1:real_roots_num
%     X1 = real_roots{i,1};
%     for j = i+1:real_roots_num
%         X2 = real_roots{j,1};
%         if ~any(isnan(intersect(X1,X2)))
%             disp('intersection')
%         end
%     end
% end

    