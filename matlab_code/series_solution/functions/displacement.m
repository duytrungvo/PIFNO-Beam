function wL = displacement(b, L)

% wL = 0;
% for n = 1:length(b)
%     wL = wL + b(n) * L^(n - 1);
% end

xL = zeros(length(b),1);
for i = 1:length(b)
    xL(i) = L^(i-1);
end

wL = sum(b .* xL);

end