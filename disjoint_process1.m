%function [real_zeros,real_roots_num,process_num] = disjoint_process1(real_roots,real_roots_num,F,n,process_num)
function process_num = disjoint_process1(F,n)
global real_roots real_roots_num

process_num = 0;
k = 0;
real_zeros = {};
for i = 1:real_roots_num
    X1 = real_roots{i,1};
    
%     [X1,zero_flag] = zero_check(F,X1,n);
%     if ~zero_flag
%         disp('WTF')
%     end
    
    new_root = true;
    for j = 1:k
        X2 = real_roots{j,1};
        Y = intersect(X1,X2);
        %[empty,Y] = emptyintersect(X1,X2);
        if any(isnan(Y))
            continue
        end
%         if any(empty)
%             continue
%         end
        [Yt,flag] = zero_check(F,Y,n);
        if flag
            %real_roots{j,1} = Yt;
            new_root = false;
            process_num = process_num + 1;
            break;
        end
    end
    if new_root
        k = k + 1;
        real_zeros{k,1} = X1;
        real_zeros{k,2} = real_roots{i,2};
    end
end
real_roots_num = k;

%以上先去掉重复的根
%下面再保证各个区间不交

for i = 1:real_roots_num
    X1 = real_zeros{i,1};
    
%     [X1,zero_flag] = zero_check(F,X1,n);
%     if ~zero_flag
%         disp('WTF')
%     end
    
    for j = i+1:real_roots_num
        X2 = real_zeros{j,1};
        Y = intersect(X1,X2);
        if any(isnan(Y))
            continue
        end
%         [empty,Y] = emptyintersect(X1,X2);
%         if any(empty)
%             continue
%         end
        X1 = sub_intval(X1,Y,n);
        real_zeros{j,1} = sub_intval(X2,Y,n);
        process_num = process_num + 1;
    end
    real_zeros{i,1} = X1;
end

% for i = 1:real_roots_num
%     X1 = real_zeros{i,1};
%     for j = i+1:real_roots_num
%         X2 = real_zeros{j,1};
%         if ~any(isnan(intersect(X1,X2)))
%             disp('intersection')
%         end
%     end
% end
real_roots = real_zeros;
