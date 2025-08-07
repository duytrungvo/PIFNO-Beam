function [w, m] = series_solution(BC, x, E, L, N, alpha, beta, form_shape)

nx = length(x);
w = zeros(1, nx);
I = zeros(1, nx);
m = zeros(1, nx);

sol_coeff = zeros(N,1); % solution coefficients
load_coeff = zeros(N,1); % load coeffs
moment_coeff = zeros(N,1); % moment of inertia coeffs

% alpha = -0.5;
load_coeff(1) = -1;  % q0
load_coeff(2) = -alpha; % q1


% cubic : rectangular section, linear 
if form_shape == "linear_rectang"
    moment_coeff(1) = 1; moment_coeff(2) = 3 * beta;
    moment_coeff(3) = 3 * beta^2; moment_coeff(4) = beta^3;
end

% rectangular section, parabolic
if form_shape == "parabolic_rectang"
    moment_coeff(1) = 1; moment_coeff(2) = 6 * beta;
    moment_coeff(3) = 15 * beta^2; moment_coeff(4) = 20 * beta^3;
    moment_coeff(5) = 15 * beta^4; moment_coeff(6) = 6 * beta^5;
    moment_coeff(7) = beta^6;
end

% rectangular section, parabolic symmetric
if form_shape == "parabolic_sym_rectang"
    h3 = beta; h2 = -h3 * L;
    moment_coeff(1) = 1; moment_coeff(2) = 3 * h2;
    moment_coeff(3) = 3 * h2^2 + 3 * h3; moment_coeff(4) = h2^3 + 6 * h2 * h3;
    moment_coeff(5) = 3 * h2^2 * h3 + 3 * h3^2; moment_coeff(6) = 3 * h2 * h3^2;
    moment_coeff(7) = h3^3;
    % I0 = 1;
    % moment_coeff = I0 * moment_coeff;
    % xx = L/3;
    % mm = dot(moment_coeff(1:7), [1, xx, xx^2, xx^3, xx^4, xx^5, xx^6])
end

% fourth order: circular section, linear
if form_shape == "linear_cir"
    moment_coeff(1) = 1; moment_coeff(2) = 4 * beta;
    moment_coeff(3) = 6 * beta^2; moment_coeff(4) = 4 * beta^3;
    moment_coeff(5) = beta^4;
end

% third order: circular hollow section, linear 
if form_shape == "linear_cir_hollow"
    % t = c*h0; c = 0.1 
    c = 0.1;
    g = 1 - c;
    moment_coeff(1) = 1 - g^4; moment_coeff(2) = 4 * (1 - g^3) * beta;
    moment_coeff(3) = 6 * (1 - g^2) * beta^2; 
    moment_coeff(4) = 4 * (1 - g) * beta^3;
end

objfun = @(x)loss_func(BC, x, sol_coeff, load_coeff, moment_coeff, E, L);

if BC == "CC"
    init_x = [-0.105, -0.02238]; % CC [1, 1] [-0.105, -0.02238]
elseif BC == "SS"
    init_x = [-0.774, 0.069]; % HH [-0.77, 1]
elseif BC == "CS"
    init_x = [-0.1162, -0.031]; % CH [-0.1519, 0.01]
else
    init_x = [-0.167, -0.075]; % CF [-0.1671, -0.074], [-0.1667,    0.0139]
end

% fin = objfun(init_x)

options = optimoptions('fminunc','Display','none','Algorithm','quasi-newton');
% options = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton');
sol_x = fminunc(objfun, init_x, options);

coeff = series_coeff(BC, sol_x, sol_coeff, load_coeff, moment_coeff, E);

for i = 1:N
    w = w + coeff(i) * x.^(i - 1);
    % q = q + load_coeff(i) * x.^(i - 1);
    I = I + moment_coeff(i) * x.^(i - 1);
    if i >= 3
        m = m + (i - 1) * (i - 2) * coeff(i) * x.^(i - 3);
    end
end
m = E * I .* m;

end

