function [real_roots,real_roots_num] = disjoint_process(real_roots,real_roots_num,F,n)

k = 0;      %�����������֤��ʵ����index
for i = 1:real_roots_num
    X1 = real_roots{i,1};
    if any(isnan(X1))
        continue;
    end
    new_root = true;       %���ڱ���Ƿ�Ϊ�µĸ�
    %��ʼ����Ϊ�µ���㣬�������ж��Ƿ����ظ�
    for j = 1:k
        X2 = real_roots{j,1};
        intersect_flag = false;     %ֻ���ڱ�ǲ����ཻ�����
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
            continue;       %����е�ά�����䲻������ֱ���ж������
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
                % ���Ż�
            end
        else
            %�����Ӧÿ��ά�ȵ����䶼���໥���������
            real_roots{j,1} = X2;   %narrowing
            new_root = false;
            break;          %˵���ظ��������ڲ�for jѭ�����ж���һ����
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

    