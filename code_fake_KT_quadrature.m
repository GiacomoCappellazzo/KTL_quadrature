% Example for KTI & KTL (Kosloff Tal-Ezer least-squares) quadrature formula.
%
% (C) G. Cappellazzo (*), W. Erb (*), F. Marchetti (*), D. Poggiali (**)
% (*)  Dipartimento di Matematica ''Tullio Levi-Civita''
% (**) PNC - Padova Neuroscience Center
% University of Padua, 2021

% Analytic function to integrate
f = @(x) 1./(1+100.*(x.^2));

% Kosloff Tal-Ezer map parameter (0 <= alpha <= 1)
alpha = 0.99;

% Number of intervals
M = 50;
% Polynomial degree of the interpolant (N <= M)
N = ceil(4*sqrt(M));
% M+1 equispaced nodes on the interval [-1,1]
nodes = linspace(-1,1,M+1);

% KTI/KTL quadrature formula
[integral,weights, coefficients] = fake_KT_quadrature(nodes, f(nodes), alpha, N);

integral  