% generate data for section 5.1 in the paper.
% solution of euler-bernoulli beam with linear moment of inertia and linear
% distributed load

clear
close all

q0 = 1;% coefficient of load q(x) = q0 + q1x
E = 1;
I0 = 1; % coefficient of linear moment of inertia I(x) = I0(h1 + h2x)^3
b = 12;  % the width of cross section
h0 = (12 * I0 / b)^(1/3); % coefficient of the depth h(x) = h0(h1 + h2x) of cross section
L = 1;

addpath ultils\

non_dim_w = abs(E*I0/q0/L^4);
non_dim_m = abs(1/q0/L^2);
num_data = 1000;
input_dim = 2;
output_dim = 2;
nx = 2^10+1;
input = zeros(num_data,nx,input_dim);
output = zeros(num_data,nx,output_dim);

alpha = zeros(num_data,1);
beta = zeros(num_data,1);
x = linspace(0,L,nx);

BC = 'CS';      % boundary condition
savefil = 'no'; % save data option

seed_number = 0;
rng(seed_number)
a = -0.9; b = -0.02;

for i = 1:num_data
    q1 = 1;
    q3 = 0;
    q2 = -rand / L;
    alpha(i) = q2;                   % parameter of load

    beta_i = (a + (b - a) * rand);   % parameter of moment of inertia
    h1 = 1;
    h2 = beta_i / L;
    beta(i) = h2;

    [C1, C2, C3, C4, U, V, W] = coefficients_linear(BC, h1, h2, -q0, q1, q2, q3, E, I0, 0, L);

    q = q0 * (q1  + q2 * x + q3 * x.^2);

    seg = segment(h1, h2, x);

    h = h0 * seg;

    I = I0 * seg.^3;

    w = C1 + C2 * log(seg) + C3 * seg ...
        + C4 ./ seg + U * seg .* log(seg) ...
        + V * seg.^2 + W * seg.^3;

    m = E * I0 * h2^2 * (-C2 * seg + 2 * C4 ...
        + U * seg.^2 + 2 * V * seg.^3 + 6 * W * seg.^4);

    input(i,:,1) = I;
    input(i,:,2) = q;
    output(i,:,1) = w;
    output(i,:,2) = m;

end

if strcmp(savefil, 'yes')
    save_data('', 'eb_beam_linear', BC, nx, seed_number, input, ...
        output, x, alpha, beta, q0, E, I0, h0, L)
end

i = 1
% plot a sample solution
plot_solution(x, input, output, "rectangular", i, alpha, beta, ...
    I0, q0, L, non_dim_w, non_dim_m)

rang_indx = 1:num_data;
% plot domain of moment of inertia
[valI, indxI] = sort(beta(rang_indx));
moment_inertia_domain = [(input(rang_indx(indxI(1)),:,1)/I0);...
    (input(rang_indx(indxI(end)),:,1)/I0) - (input(rang_indx(indxI(1)),:,1)/I0)]';

plot_domain_v2(x, input(rang_indx(indxI(1:round(num_data/9):end)),:,1),...
    beta(rang_indx(indxI(1:round(num_data/9):end))), moment_inertia_domain, '$\phi$')

% plot domain of load
[valq, indxq] = sort(alpha(rang_indx));
load_domain = [(input(rang_indx(indxq(1)),:,2)/q0);...
    (input(rang_indx(indxq(end)),:,2)/q0) - (input(rang_indx(indxq(1)),:,2)/q0)]';

plot_domain_v2(x, input(rang_indx(indxq(1:round(num_data/9):end)),:,2), ...
    alpha(rang_indx(indxq(1:round(num_data/9):end))), load_domain, '$\rho$')
