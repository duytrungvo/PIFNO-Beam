function dwL = slope(b, L)

% dwL = 0;
% for n = 2: length(b)
%     dwL = dwL + (n - 1) * b(n) * L^(n - 2);
% end

dxL = zeros(length(b),1);
for i = 2:length(b)
    dxL(i) = (i - 1) * L^(i-2);
end

dwL = sum(b .* dxL);

end