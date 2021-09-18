% Relative Errors for KTL quadrature formula.
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

int_trap = [];
int_KTL_1 = [];
int_KTL_2 = [];
int_KTL_var = [];

% Number of intervals
ppp = 50:25:500;

for M = ppp
    % M+1 equispaced nodes on the interval [-1,1]
    xx = linspace(-1,1,M+1);
    
    % KTI quadrature formula (alpha = 1, degree=M)
    w = weights_KTI(xx,1);
    % approximation of the integral
    int_trap = [int_trap,f(xx)*w];
    
    % KTL quadrature formula (alpha = 0.9, degree=4*sqrt(M))
    [weights, coefficients] = weights_KTL(xx, f(xx), 0.9, ceil(4*sqrt(M)));
    % approximation of the integral
    int_KTL_1 = [int_KTL_1,weights*coefficients];
    
    % KTL quadrature formula (alpha = 0.7, degree=4*sqrt(M))
    [weights, coefficients] = weights_KTL(xx, f(xx), 0.7, ceil(4*sqrt(M)));
    % approximation of the integral
    int_KTL_2 = [int_KTL_2,weights*coefficients];
    
    % KTL quadrature formula (alpha = dynamic, degree=M/2)
    alp = 1+((2*log(tol))/(pi*M/2));
    [weights, coefficients] = weights_KTL(xx, f(xx), alp, ceil(M/2));
    % approximation of the integral
    int_KTL_var = [int_KTL_var,weights*coefficients];
end

% Relative Errors
% alpha = 1, degree=M
err_rel_trap = abs(int_trap - int_f )./abs(int_f);
% alpha = 0.9, degree=4*sqrt(M)
err_rel_KTL_1 = abs(int_KTL_1 - int_f )./abs(int_f);
% alpha = 0.7, degree=4*sqrt(M)
err_rel_KTL_2 = abs(int_KTL_2 - int_f )./abs(int_f);
% alpha = dynamic, degree=M/2
err_rel_KTL_var = abs(int_KTL_var - int_f )./abs(int_f);

% Begin Plot --------------------------------------------------------------
P = semilogy(ppp,err_rel_trap,'blue',...
    ppp,err_rel_KTL_1,'black',...
    ppp,err_rel_KTL_2,'magenta',...
    ppp,err_rel_KTL_var,'red');

P(1).LineWidth = 3;
P(2).LineWidth = 3;
P(3).LineWidth = 3;
P(4).LineWidth = 3;
grid on;

legend({'Trapezoidal Rule','\alpha = 0.9','\alpha = 0.7',...
    '\alpha = dynamic'},'Location','southwest');

axis([50 500 1e-14 1]);
axis square;
% End Plot ----------------------------------------------------------------