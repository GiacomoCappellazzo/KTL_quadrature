% Example for KTI (Kosloff Tal-Ezer Interpolation) quadrature formula.
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
M = 20;
% M+1 equispaced nodes on the interval [-1,1]
nodes = linspace(-1,1,M+1);

% KTI quadrature formula
w = weights_KTI(nodes, alpha);

% approximation of the integral
int_w = f(nodes)*w