% Relative Errors for KTL quadrature formula (perturbed nodes).
%
% (C) G. Cappellazzo (*), W. Erb (*), F. Marchetti (*), D. Poggiali (**)
% (*)  Dipartimento di Matematica ''Tullio Levi-Civita''
% (**) PNC - Padova Neuroscience Center
% University of Padua, 2021

% Analytic function to integrate
f = @(x) 1./(1+100.*(x.^2));
% Exact value of the integral
int_f = integral(f,-1,1,'AbsTol',1e-15,'RelTol',1e-9);

% tollerance --> alpha dynamic choice
tol = 10^(-12);

int_KTL_1 = [];
int_KTL_2 = [];

% Number of intervals
ppp = 50:25:500;

for M = ppp
    % M+1 equispaced nodes on the interval [-1,1]
    xx = linspace(-1,1,M+1);
    
    % M+1 perturbed nodes on the interval [-1,1]
    % Starting nodes: M+1 equispaced nodes on the interval [-1,1] 
    xx_pert = [];
    xx_pert(1) = 0;
    xx_pert(M+1) = 0;
    % Uniform random variable on [-1/M,1/M]
    xx_pert(2:M) = (2.*rand(1,M-1)-1).*(1/M);
    xx_pert = xx+xx_pert;
    
    % KTL quadrature formula (alpha = dynamic, degree=M/2)
    alp = 1+((2*log(tol))/(pi*M/2));
    
    [weights, coefficients] = weights_KTL(xx, f(xx), alp, ceil(M/2));
    % approximation of the integral (equispaced nodes)
    int_KTL_1 = [int_KTL_1,weights*coefficients];
    
    [weights, coefficients] = weights_KTL(xx_pert, f(xx_pert), alp, ceil(M/2));
    % approximation of the integral (perturbed nodes)
    int_KTL_2 = [int_KTL_2,weights*coefficients];
    
end

% Relative Errors
% alpha = dynamic, degree=M/2, equispaced nodes
err_rel_KTL_1 = abs(int_KTL_1 - int_f )./abs(int_f);
% alpha = dynamic, degree=M/2, perturbed nodes
err_rel_KTL_2 = abs(int_KTL_2 - int_f )./abs(int_f);

% Begin Plot --------------------------------------------------------------
P = semilogy(ppp,err_rel_KTL_1,'red',...
    ppp,err_rel_KTL_2,'blue');

P(1).LineWidth = 3;
P(2).LineWidth = 3;
grid on;
legend({'Equispaced Nodes','Perturbed Nodes'},'Location','southwest');

axis([50 500 1e-14 1]);
axis square;
% End Plot ----------------------------------------------------------------