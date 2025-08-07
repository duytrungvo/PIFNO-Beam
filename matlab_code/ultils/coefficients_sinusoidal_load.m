function [C1, C2, C3, C4, U, V] = coefficients_sinusoidal_load(BC, q0, q1, q2, omega, E, I0, a, b)
T = zeros(4);
r = zeros(4,1);

Q0 = q0 * q1 / E / I0;
Q1 = q0 * q2 / E / I0;
U = Q0 / 24;
V = Q1 / omega^4;

if BC == "CC"
% w(0) = 0
T(1,1) = 1; T(1,2) = a; T(1,3) = a^2; T(1,4) = a^3;
r(1,1) = -U * a^4 - V * sin(omega * a);

% w'(0) = 0
T(2,2) = 1; T(2,3) = 2 * a; T(2,4) = 3 * a^2;
r(2,1) = -4 * U * a^3 - omega * V * cos(omega * a);

% w(L) = 0
T(3,1) = 1; T(3,2) = b; T(3,3) = b^2; T(3,4) = b^3;
r(3,1) = -U * b^4 - V * sin(omega * b);

% w'(L) = 0
T(4,2) = 1; T(4,3) = 2 * b; T(4,4) = 3 * b^2;
r(4,1) = -4 * U * b^3 - omega * V * cos(omega * b);
end

if BC == "CS"
% w(0) = 0
T(1,1) = 1; T(1,2) = a; T(1,3) = a^2; T(1,4) = a^3;
r(1,1) = -U * a^4 - V * sin(omega * a);

% w'(0) = 0
T(2,2) = 1; T(2,3) = 2 * a; T(2,4) = 3 * a^2;
r(2,1) = -4 * U * a^3 - omega * V * cos(omega * a);

% w(L) = 0
T(3,1) = 1; T(3,2) = b; T(3,3) = b^2; T(3,4) = b^3;
r(3,1) = -U * b^4 - V * sin(omega * b);

% M(L) = 0
T(4,3) = 2; T(4,4) = 6 * b;
r(4,1) = -12 * U * b^2 + omega^2 * V * sin(omega * b);
end

if BC == "CF"
% w(0) = 0
T(1,1) = 1; T(1,2) = a; T(1,3) = a^2; T(1,4) = a^3;
r(1,1) = -U * a^4 - V * sin(omega * a);

% w'(0) = 0
T(2,2) = 1; T(2,3) = 2 * a; T(2,4) = 3 * a^2;
r(2,1) = -4 * U * a^3 - omega * V * cos(omega * a);

% M(L) = 0
T(3,3) = 2; T(3,4) = 6 * b;
r(3,1) = -12 * U * b^2 + omega^2 * V * sin(omega * b);

% V(L) = 0
T(4,4) = 6;
r(4,1) = -24 * U * b + omega^3 * V * cos(omega * b);
end

if BC == "SS"
% w(0) = 0
T(1,1) = 1; T(1,2) = a; T(1,3) = a^2; T(1,4) = a^3;
r(1,1) = -U * a^4 - V * sin(omega * a);

% M(0) = 0
T(2,3) = 2; T(2,4) = 6 * a;
r(2,1) = -12 * U * a^2 + omega^2 * V * sin(omega * a);

% w(L) = 0
T(3,1) = 1; T(3,2) = b; T(3,3) = b^2; T(3,4) = b^3;
r(3,1) = -U * b^4 - V * sin(omega * b);

% M(L) = 0
T(4,3) = 2; T(4,4) = 6 * b;
r(4,1) = -12 * U * b^2 + omega^2 * V * sin(omega * b);
end

x = T\r;
C1 = x(1); C2 = x(2); C3 = x(3); C4 = x(4);
end