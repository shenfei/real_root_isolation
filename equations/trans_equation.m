function P = trans_equation(input_file)

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
