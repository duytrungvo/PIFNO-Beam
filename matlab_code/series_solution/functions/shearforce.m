function mL = shearforce(b, L)

d3xL = zeros(length(b), 1);
for i = 4:length(b)
    d3xL(i) = (i - 1) * (i - 2) * (i - 3) * L^(i - 4);
end

mL = sum(b .* d3xL);

end