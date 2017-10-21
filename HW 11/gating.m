function outputIndex = gating( input, max )
if input>1 && input<max
    outputIndex=input;
elseif input<=1
    outputIndex=1;
else
    outputIndex = max;
end

