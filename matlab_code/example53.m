% generate data for section 5.3 in the paper.

clear
close all

q0 = 1;    % coefficient of load q(x) = q0 + q1x
E = 1;
I0 = 1;
b = 12;  % the width of cross section
h0 = (12 * I0 / b)^(1/3); % coefficient of the depth h(x) = h0(h1 + h2x) of cross section
L = 1;

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

BC = 'CC';      % boundary condition
savefil = 'no'; % save data option

seed_number = 0;
rng(seed_number)
a = -0.1;  b = 0.1; % -0.12 a < b < 0.22

for i = 1:num_data
    if mod(i, 10) == 0
        fprintf('iteration = %d\n', i)
    end

    q1 = 1;
    q3 = 0; 
    q2 = -0.5; 
    alpha(i) = q2; 

    h1 = 1;
    h3 = (a + (b - a) * rand) * 4 / L^2;    % parameter of moment of inertia
    h2 = -h3 * L;
    beta(i) = h3;

    q = q0 * (q1 + q2 * x + q3 * x.^2);     % q0 * rho
    h = h0 * (h1 + h2 * x + h3 * x.^2);     % h0 * psi
    I = I0 * (h1 + h2 * x + h3 * x.^2).^3;  % I0 * phi

    [w, m] = series_solution(BC, x, E, L, 1000, q2, h3, "parabolic_sym_rectang");

    input(i,:,1) = I;
    input(i,:,2) = q;
    output(i,:,1) = w;
    output(i,:,2) = m;

end

if strcmp(savefil,'yes')
    save_data('', 'eb_beam_parabolic_sym', BC, nx, seed_number, ...
        input, output, x, alpha, beta, q0, E, I0, h0, L)
end

i = 1
% plot a sample solution
plot_solution(x, input, output, "rectangular", i, alpha, beta, ...
    I0, q0, L, non_dim_w, non_dim_m)

rang_indx = 1:num_data;
% plot domain of moment of inertia
[valI, indxI] = sort(beta);
moment_inertia_domain = [(input(indxI(end),:,1)/I0);...
    (input(indxI(1),:,1)/I0) - (input(indxI(end),:,1)/I0)]';

plot_domain_v2(x, input(rang_indx(indxI(1:round(num_data/9):end)),:,1), ...
    beta(rang_indx(indxI(1:round(num_data/9):end))), moment_inertia_domain, '$\phi$')
