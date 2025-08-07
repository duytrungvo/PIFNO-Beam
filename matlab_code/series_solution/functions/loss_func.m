function res = loss_func(BC, x, b, c, d, E, L)

b = series_coeff(BC, x, b, c, d, E);

if BC == "CC"
    con1 = displacement(b, L);
    con2 = slope(b, L);    
end

if BC == "CS"
    con1 = displacement(b, L);
    con2 = moment(b, L);
end

if BC == "CF"
    con1 = moment(b, L);
    con2 = shearforce(b, L);
end

if BC == "SS"
    con1 = displacement(b, L);
    con2 = moment(b, L);
end

res = con1^2 + con2^2;

% obj1 = sum(b)^2;
% obj2 = sum(b.* (0:length(b)-1)')^2;
% res = obj1 + obj2;
end