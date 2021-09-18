function [weights] = weights_KTI(nodes, alpha)

% KTI (Kosloff Tal-Ezer Interpolation) quadrature formula.
%
% (C) G. Cappellazzo (*), W. Erb (*), F. Marchetti (*), D. Poggiali (**)
% (*)  Dipartimento di Matematica ''Tullio Levi-Civita''
% (**) PNC - Padova Neuroscience Center
% University of Padua, 2021
% ---------------------------------------------------------------------
% INPUT:
% nodes = [x_0, ..., x_M] : quadrature nodes on a compact interval
% alpha                   : KT map parameter (0 <= alpha <= 1)
%
% OUTPUT:
% weights [M+1 x 1]       : M+1 interpolatory quadrature weights

    % Polynomial degree of the interpolant
    M = length(nodes) - 1;
    
    % normalizing the nodes to [-1,1]
    scale = (max(nodes) - min(nodes));
    nodes = (nodes - min(nodes)) / scale;
    nodes = 2 * nodes - 1;
    
    %Fake Nodes Approach (FNA) with KT map
   
    if alpha==0
        % KT map
        m_alpha = @(x) x;
        % moments vector : integrals of basis functions
        intcosw = zeros(M+1,1);
        for i = 0:2:M
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
        intcosw = intcos(1:M+1)';
    end
    
    % Fake nodes
    mnodes = m_alpha(nodes);
    
    % Quadrature matrix with Chebyshev polynomials
    NDCT = zeros(M+1,M+1);
    for i = 1:M+1
        NDCT(i,:) = cos((i-1)*acos(mnodes));
    end
    
    % KTI quadrature weights
    weights = NDCT\intcosw;

end