function mL = moment(b, L)

d2xL = zeros(length(b), 1);
for i = 3:length(b)
    d2xL(i) = (i - 1) * (i - 2) * L^(i - 3);
end

mL = sum(b .* d2xL);

end