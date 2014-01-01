function idx_array = resort_var(Fin,vardata,n)

idx_array = zeros(n,1);
fid = fopen(vardata,'r');
F_var = symvar(Fin);
hom_var = cell(n,1);
for i = 1:n
    hom_var{i} = fgetl(fid);
end
for i = 1:n
    for j = 1:n
        if F_var{i} == hom_var{j}
            idx_array(i) = j;
            break;
        end
    end
end
