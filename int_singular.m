function [flag] = int_singular(A)

n = dim(A);
for i = 1:n-1
    j = i;
    while j <= n && in0(0,A(j,i))
        j = j + 1;
    end
    if j > n
        flag = true;
        return
    end
    if j ~= i
        temp = A(i,i:n);
        A(i,i:n) = A(j,i:n);
        A(j,i:n) = temp;
    end
    for j = i+1:n
        A(j,i:n) = A(j,i:n) - A(i,i:n) * A(j,i) / A(i,i);
    end
end
intdet = 1;
for i = 1:n
    intdet = intdet * A(i,i);
end
if in0(0,intdet)
    flag = true;
else
    flag = false;
end

    