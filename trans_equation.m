function [P,file_name] = trans_equation(input_file)

flen = length(input_file);
j = flen;
k = flen;
while j>0 && input_file(j) ~= '/'
    j = j - 1;
end
while k>0 && input_file(k) ~= '.'
    k = k - 1;
end
file_name = input_file(j+1:k-1);

fin = fopen(input_file);
P = {};
i = 1;
line = fgetl(fin);     % read the '{'
equ = '';
while true
    line = fgetl(fin);
    if line == '}'
        break;
    end
    len = length(line);
    if len == 0
        continue;
    end
    if line(len) == ';'
        equ = strcat(equ,line(1:len-1));
        P{i} = equ;
        i = i+1;
        equ = '';
    else
        equ = strcat(equ,line);
    end
end
fclose(fin);
P = P';
