function v = vector_parameterize(v_bar)    
a = v_bar(1);
b = v_bar(2:end);
v = (2/mysinc(acos(a)))*b;

if(norm(v) > pi)
    v = (1 - (2*pi/norm(v))*ceil((norm(v)-pi)/(2*pi)))*v;
end
end