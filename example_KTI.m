% Relative Errors for KTI quadrature formula.
%
% (C) G. Cappellazzo (*), W. Erb (*), F. Marchetti (*), D. Poggiali (**)
% (*)  Dipartimento di Matematica ''Tullio Levi-Civita''
% (**) PNC - Padova Neuroscience Center
% University of Padua, 2021

% Analytic function to integrate
f = @(x) 1./(1+100.*(x.^2));
% Exact value of the integral
int_f = integral(f,-1,1,'AbsTol',1e-15,'RelTol',1e-9);

int_trap = [];
int_KTI_1 = [];
int_KTI_2 = [];

% Number of intervals
ppp = 50:25:500;

for M = ppp
    % M+1 equispaced nodes on the interval [-1,1]
    xx = linspace(-1,1,M+1);
    
    % KTI quadrature formula (alpha = 1)
    w = weights_KTI(xx, 1);
    % approximation of the integral
    int_trap = [int_trap,f(xx)*w];
    
    % KTI quadrature formula (alpha = 0.99)
    w = weights_KTI(xx,0.99);
    % approximation of the integral
    int_KTI_1 = [int_KTI_1,f(xx)*w];
    
    % KTI quadrature formula (alpha = 0.98)
    w = weights_KTI(xx,0.98);
    % approximation of the integral
    int_KTI_2 = [int_KTI_2,f(xx)*w];
end

% Relative Errors
% alpha = 1
err_rel_trap = abs(int_trap - int_f )./abs(int_f);
% alpha = 0.99
err_rel_KTI_1 = abs(int_KTI_1 - int_f )./abs(int_f);
% alpha = 0.98
err_rel_KTI_2 = abs(int_KTI_2 - int_f )./abs(int_f);

% Begin Plot --------------------------------------------------------------
P = semilogy(ppp,err_rel_trap,'blue',...
    ppp,err_rel_KTI_1,'magenta',...
    ppp,err_rel_KTI_2,'black');

P(1).LineWidth = 3;
P(2).LineWidth = 3;
P(3).LineWidth = 3;
grid on;

legend({'Trapezoidal Rule','\alpha = 0.99','\alpha = 0.98'},...
    'Location','southwest');

axis([50 500 1e-11 inf]);
axis square;
% End Plot ----------------------------------------------------------------