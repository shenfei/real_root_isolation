function [rootCell,complex_roots_num,real_roots,real_roots_num] = real_root_isolate(equ_path)

[P,file_name] = trans_equation(equ_path);
% P = {
%     'x1*x5+x5*x1*x2+x5*x2*x3+x3*x4*x5-1';
%     'x2*x5+x5*x1*x3+x5*x2*x4-2';
%     'x3*x5+x4*x5*x1-3';
%     'x4*x5-4';
%     'x1+x2+x3+x4+1';
%     };

format long

total_time = cputime;   %程序总时间
%preprocess begin
n = length(P);
H = '{';
for k = 1:n
    H = strcat(H,[P{k},';']);
end
H = strcat(H,'}');
Q = H;
Q(1) = '[';
Q(length(Q)) = ']';
Fin = inline(Q);

Fsym = sym(Fin);
%Jac = jacobian(Fsym);
for k = 1:n
    Fsym(k) = horner(Fsym(k));
end
Fin = inline(Fsym);

var_num = nargin(Fin);
sfun = 'F =@(X)Fin(X(1)';
for k = 2:var_num
    s = sprintf(',X(%d)',k);
    sfun = strcat(sfun,s);
end
sfun = strcat(sfun,');');
eval(sfun);

hom_time = cputime;
[rootCell,complex_roots_num] = HOM4PS2(n,H);
hom_time = cputime - hom_time;      %记录同伦的时间

% Resort the variables' order:
extra_time = cputime;

cd Hom4PS2\
command = sprintf('hom4ps2.bat ..\\%s 1',equ_path);
system(command);
cd ..
command = sprintf('var_orders.exe < hom4ps2/data.roots > data.var');
system(command);

idx_array = resort_var(Fin,'data.var',n);
for i = 1:complex_roots_num
    z = rootCell{i,1};
    x = zeros(n,1);
    for j = 1:n
        x(j) = z(idx_array(j));
    end
    rootCell{i,1} = x;
end

extra_time = cputime - extra_time;
total_time = total_time + extra_time;
% Resort ends
% Preprocess ends

global real_roots_num real_roots
real_roots_num = 0;
real_roots = {};

width_array = zeros(complex_roots_num,1);   %记录初始区间半径
width_time_array = zeros(complex_roots_num,1);  %记录计算半径的时间
Krawcyzk_time_array = zeros(complex_roots_num,1);   %记录K迭代时间
iter_count = zeros(complex_roots_num,1);    %记录迭代次数
priori_judge = 0;       %提前判断的复根个数
intval_time = cputime;  %区间算术总时间

for k = 1:complex_roots_num
    %width = init_width(F,Fin,rootCell{k,1},norm(F(rootCell{k,1}),1),Jac,n);
    %width = init_width(F,Fin,rootCell{k,1},rootCell{k,2},Jac,n);
    %width = init_width4(F,rootCell{k,1});
    temp = cputime;
    width = init_width5(F,rootCell{k,1},n);
    %[width,r2] = init_width7(F,rootCell{k,1},n);
    %r3 = init_width5(F,rootCell{k,1},n)
    if width < eps
        width = max(rootCell{k,3} * 1e-15,1e-5);
        %disp('estimation fail')
    end
    %k
    width_time_array(k) = cputime - temp;
    width_array(k) = width;
    
    if any(abs(imag(rootCell{k,1})) > width)
        priori_judge = priori_judge + 1;
        continue
    end
    X0 = midrad( real(rootCell{k,1}), width);
    %infsup(X0)
    
    global iter_num
    iter_num = 0;
    temp = cputime;
    %[real_roots,real_roots_num,flag] = Krawczyk_Moore(F,X0,n,real_roots,real_roots_num,k);
    flag = Krawczyk_Moore(F,X0,n,k);
%     if ~flag
%         infsup(X0)
%         k
%     end
    Krawcyzk_time_array(k) = cputime - temp;
    iter_count(k) = iter_num;
end
intval_time = cputime - intval_time;

post_process_time = cputime;
%process_num: 真实处理的相交情况次数
%[real_roots,real_roots_num,process_num] = disjoint_process1(real_roots,real_roots_num,F,n);
process_num = disjoint_process1(F,n);

tau = 1e-6;
%[real_roots,real_roots_num] = narrowing(real_roots,real_roots_num,F,n,tau);
narrowing(F,n,tau);
post_process_time = cputime - post_process_time;

total_time = cputime - total_time;

for k = 1:real_roots_num
    infsup(real_roots{k,1})
end

fprintf('The order of variables:\n')
disp(symvar(Fin));

fprintf('The number of real roots: %d\n',real_roots_num)

%debug1(rootCell,real_roots,real_roots_num,F);

%以下将结果输出至文件
%fout = fopen('RealRootsData.roots','w');
fout = fopen(strcat(file_name,'.roots'),'w');
fprintf(fout,'SUMMURY:\n');
div_line = '--------------------\n';
fprintf(fout,div_line);
fprintf(fout,'The # of total roots of Homotopy: %d\n',complex_roots_num);
fprintf(fout,'The # of real roots of Homotopy: \n');
fprintf(fout,'The # of real roots verified: %d\n',real_roots_num);
fprintf(fout,'The # of real roots of Discoverer: \n');
fprintf(fout,'The complex roots prior detected: %d / %d\n',priori_judge,complex_roots_num);
fprintf(fout,'The # of repeat roots processed: %d\n\n',process_num);
fprintf(fout,'The average iteration times: %.3f\n',sum(iter_count)/real_roots_num);

fprintf(fout,'The total time: %f\n',total_time);
fprintf(fout,'The time of Homotopy: %f\n',hom_time);
fprintf(fout,'The time of interval arithmetic: %f\n',intval_time);
fprintf(fout,'The time of disjoint and narrowing process: %f\n',post_process_time);
fprintf(fout,'The time of Discoverer: \n\n');

%idx_array = [real_roots{:,2}];
fprintf(fout,'The average radius of result interval: %e\n',mean(max(rad([real_roots{:,1}]))));

fprintf(fout,'The max initial width: %e\n',max(width_array));
%fprintf(fout,'The min initial width: %f\n',min(width_array));
%fprintf(fout,'The average of initial width: %e\n\n',mean(width_array(idx_array)));
fprintf(fout,'The average of initial width: %e\n\n',mean(width_array));

fprintf(fout,'The max width calc time: %f\n',max(width_time_array));
%fprintf(fout,'The min width calc time: %f\n',min(width_time_array));
%fprintf(fout,'The average of width calc time: %f\n\n',mean(width_time_array(idx_array)));
fprintf(fout,'The average of width calc time: %f\n\n',mean(width_time_array));

fprintf(fout,'The max Krawczyk iteration time: %f\n',max(Krawcyzk_time_array));
%fprintf(fout,'The min Krawcyzk iteration time: %f\n',min(Krawcyzk_time_array));
temp_array = Krawcyzk_time_array(Krawcyzk_time_array > eps);
fprintf(fout,'The average of Krawczyk iteration time: %f\n',mean(temp_array));
%fprintf(fout,'The average of Krawcyzk iteration time: %f\n',mean(Krawcyzk_time_array(idx_array)));

fprintf(fout,div_line);

fprintf(fout,'\nISOLATED REAL ROOTS\n');
fprintf(fout,div_line);

for i = 1:real_roots_num
    fprintf(fout,'%d:\n',i);
    k = real_roots{i,2};
    X = sprintf(infsup(real_roots{i,1}));
    len = length(X);
    step = len / n;
    j = 1;
    strX = '';
    while j <= len
        strX = strcat(strX,sprintf(X(j:j + step - 1)));
        strX = strcat(strX,'\n');
        j = j + step;
    end
    strX = sprintf(strX);
    fprintf(fout,strX);
    fprintf(fout,'Initial width: %e\n',width_array(k));
    fprintf(fout,'Time of width calc: %f\n',width_time_array(k));
    fprintf(fout,'Time of Krawczyk iteration: %f\n',Krawcyzk_time_array(k));
    fprintf(fout,div_line);
    fprintf(fout,'\n');
end

fprintf(fout,'The order of variables:\n');
F_vars = symvar(Fin);
for i = 1:n
    fprintf(fout,'%s\n',F_vars{i});
end
fprintf(fout,div_line);

fclose(fout);
