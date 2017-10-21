function output = sub (input,max)
if input>1 && input<max
    output=input;
elseif input<=1
    output=1;
else output = max;
end