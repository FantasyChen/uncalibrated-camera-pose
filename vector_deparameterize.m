function v_bar = vector_deparameterize(v)
    a = cos(norm(v)/2);
    b = mysinc(norm(v)/2)/2 * v;
    v_bar = [a; b];
end