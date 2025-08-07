% generate data for section 5.4 in the paper.

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

BC = 'CC';      % boundary condition
savefil = 'no'; % save data option

seed_number = 0;
rng(seed_number)
a = -0.5; b = 0.5;

for i = 1:num_data

    q1 = 1.0;
    q2 = a + (b - a) * rand;
    omega = 3 * pi / L;
    alpha(i) = q2;                      % parameter of load
    
    h1 = 1;
    h2 = 0;
    beta(i) = h2;

    [C1, C2, C3, C4, U, V] = coefficients_sinusoidal_load(BC, ...
        -q0, q1, q2,omega, E, I0, 0, L);

    q = q0 * (q1 + q2 * sin(omega * x));
    seg = segment(h1, h2, x);

    I = I0 * seg;

    w = C1 + C2 * x + C3 * x.^2 + C4 * x.^3 ...
        + U * x.^4 + V * sin(omega * x);

    m = E * I0 * (2 * C3 + 6 * C4 * x ...
        + 12 * U * x.^2 - omega^2 * V * sin(omega * x));

    input(i,:,1) = I;
    input(i,:,2) = q;
    output(i,:,1) = w;
    output(i,:,2) = m;

end

if strcmp(savefil, 'yes')
    save_data('const', 'eb_beam_sin_load', BC, nx, seed_number, ...
        input, output, x, alpha, beta, q0, E, I0, h0, L)
end

i = 1
% plot a sample solution
plot_solution(x, input, output, "rectangular", i, alpha, beta, ...
    I0, q0, L, non_dim_w, non_dim_m)

rang_indx = 1:num_data;
% plot domain of load
[valq, indxq] = sort(alpha);
load_domain = [(input(indxq(1),:,2)/q0);...
    (input(indxq(end),:,2)/q0) - (input(indxq(1),:,2)/q0)]';
plot_domain_v2(x, input(rang_indx(indxq(1:round(num_data/9):end)),:,2), ...
    alpha(rang_indx(indxq(1:round(num_data/9):end))), load_domain, '$\rho$')
