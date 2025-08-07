function res = coeff(n, p)

res = N2(n - p + 2) * I2(p -2) + 2 * N3(n - p + 1) * I1(p - 1) + N4(n - p);
end