% generate data for section 5.2 in the paper.
% solution of euler-bernoulli circular hollow section beam with
% linear moment of inertia and linear distributed load

clear
close all

q0 = 1;    % coefficient of load q(x) = q0 + q1x
E = 1;
I0 = 1;     % coefficient of linear moment of inertia I(x) = I0[(h1 + h2x)^4 - (h1 - c + h2x)^4]
h0 = (64 * I0 / pi)^(1/4); % coefficient of the depth h(x) = h0[(h1 + h2x) - (h1 - c + h2x)] of cross section
L = 1;
c = 0.1;    % t = c * h0 is thickness of cross section

addpath ultils\
addpath series_solution\functions\

non_dim_w = abs(E*I0/q0/L^4);
non_dim_m = abs(1/q0/L^2);
num_data = 100;
input_dim = 2;
output_dim = 2;
nx = 2^10+1;
input = zeros(num_data,nx,input_dim);
% output = 1e-3*ones(no_data,nx,output_dim);
output = zeros(num_data,nx,output_dim);

alpha = zeros(num_data,1);
beta = zeros(num_data,1);
x = linspace(0,L,nx);

BC = 'CF';      % boundary condition
savefil = 'no'; % save data option

seed_number = 0;
rng(seed_number)
a = -0.5;  b = -0.02;

for i = 1:num_data
    if mod(i, 10) == 0
        fprintf('iteration = %d\n', i)
    end

    q1 = 1;
    q3 = 0; 
    q2 = -0.5; 
    alpha(i) = q2;

    beta_i = (a + (b - a) * rand);      % parameter of moment of inertia
    h1 = 1;
    h2 = beta_i / L;
    beta(i) = h2;

    q = q0 * (q1 + q2 * x + q3 * x.^2);  % q0 * rho

    seg_in = segment(h1-c, h2, x);
    seg_out = segment(h1, h2, x);
    h = h0 * seg_out;
    I = I0 * (seg_out.^4 - seg_in.^4);

    [w, m] = series_solution(BC, x, E, L, 1000, q2, h2, "linear_cir_hollow");

    input(i,:,1) = I;
    input(i,:,2) = q;
    output(i,:,1) = w;
    output(i,:,2) = m;

end

if strcmp(savefil, 'yes')
    save_data('', 'eb_beam_linear_cir_hollow', BC, nx, seed_number, ...
        input, output, x, alpha, beta, q0, E, I0, h0, L)
end

i = 1
% plot a sample solution
plot_solution(x, input, output, "circular", i, alpha, beta, ...
    I0, q0, L, non_dim_w, non_dim_m)


rang_indx = 1:num_data;
% plot domain of moment of inertia
[valI, indxI] = sort(beta(rang_indx));
moment_inertia_domain = [(input(rang_indx(indxI(1)),:,1)/I0);...
    (input(rang_indx(indxI(end)),:,1)/I0) - (input(rang_indx(indxI(1)),:,1)/I0)]';

plot_domain_v2(x, input(rang_indx(indxI(1:round(num_data/9):end)),:,1), ...
    beta(rang_indx(indxI(1:round(num_data/9):end))), moment_inertia_domain, '$\phi$')
