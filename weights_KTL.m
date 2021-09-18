function [weights, coefficients] = weights_KTL(nodes, values, alpha, degree)

% KTL (Kosloff Tal-Ezer least-squares) quadrature formula.
%
% (C) G. Cappellazzo (*), W. Erb (*), F. Marchetti (*), D. Poggiali (**)
% (*)  Dipartimento di Matematica ''Tullio Levi-Civita''
% (**) PNC - Padova Neuroscience Center
% University of Padua, 2021
% -------------------------------------------------------------------------
% INPUT:
% nodes = [x_0, ..., x_M]        : quadrature nodes on a compact interval
% values = [f(x_0), ..., f(x_M)] : function samples on quadrature nodes
% alpha                          : KT map parameter (0 <= alpha <= 1)
% degree = N                     : polynomial degree of the interpolant
%
% OUTPUT:
% weights [1 x N+1]         : integrals of T_i(M_alpha(x)) i=0, ..., degree
%                             T_i     --> i-th Chebyshev Polynomial
%                             M_alpha --> KT map 
% coefficients [N+1 x 1]    : degree+1 linear combiners of the approximant

    % Polynomial degree of the interpolant (N <= M)
    M = length(nodes) - 1;
    N = degree;

    % normalizing the nodes to [-1,1]
    scale = (max(nodes) - min(nodes));
    nodes = (nodes - min(nodes)) / scale;
    nodes = 2 * nodes - 1;
    
    if alpha==0
        % KT map
        m_alpha = @(x) x;
        % moments vector : integrals of basis functions
        intcosw = zeros(1,N+1);
        for i = 0:2:N
            intcosw(i+1) = 2/(1-i^2);
        end
    else
        % KT map
        beta = sin(alpha*(pi/2));
        m_alpha = @(x) (sin((alpha*(pi/2)).*x))./(sin(alpha*(pi/2)));
        % auxiliary function for the discrete cosine transform (DCT)
        g_alpha = @(t) (sin(t)./sqrt(1/beta^2 - cos(t).^2)/alpha);
        
        % integrals of DCT moments
        n = 1e06;
        h = pi/n;
        x = linspace(0,pi,n+1) + h/2;
        x = x(1:n);
        intcos = dct(g_alpha(x))*sqrt(2)/sqrt(n);
        intcos(1) = intcos(1)*sqrt(2);
        % moments vector : integrals of basis functions
        intcosw = intcos(1:N+1);
    end
    
    % Fake nodes
    mnodes = m_alpha(nodes);
    
    % Interpolation matrix with Chebyshev polynomials
    NDCT = zeros(M+1,N+1);
    for i = 1:N+1
        NDCT(:,i) = cos((i-1)*acos(mnodes'));
    end
    
    % least-squares weights
    mu = zeros(1,M+1);
    mu(1) = (1/2)*(asin(mnodes(2))-asin(-1));
    mu(M+1) = (1/2)*(asin(1)-asin(mnodes(M)));
    for i = 2:M
        mu(i) = (1/2)*(asin(mnodes(i+1))-asin(mnodes(i-1)));
    end
    mu = sqrt(mu);
    W = diag(mu);
    
    % least-squares weighted problem : find an approximant
    A = W*NDCT;    
    b = W*values';

    % linear combiners of the approximant
    coefficients = A\b;
    % moments vector : integrals of basis functions
    weights = intcosw;