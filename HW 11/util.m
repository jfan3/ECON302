function ut = util( c )


global gamma
if c<0
    c = 10^-8;
end
ut = c^(1-gamma)/(1-gamma);
end

