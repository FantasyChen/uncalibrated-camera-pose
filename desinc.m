function y = desinc(x)

if(x == 0)
    y = 0;
else
    y = (cos(x)/x) - (sin(x)/(x^2));
end
end

