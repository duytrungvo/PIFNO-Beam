function res = series_coeff(BC, x, b, c, d, E)

N = length(b);
if BC == "CC" || BC == "CH" || BC == "CF"
    b(3) = x(1);    % a2
    b(4) = x(2);    % a3
end
if BC == "HH"
    b(2) = x(1);    % a1
    b(4) = x(2);    % a3
end
for n = 1:N-4
    sum_np = 0;
    for p = 1:n+1
        sum_np = sum_np + coeff(n-1, p) * b(n - p + 4) * d(p + 1);
    end

    b(n + 4) = 1 / N4(n-1) / d(1) * (c(n) / E - sum_np);
end
res = b;
end