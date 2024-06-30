function IMEXRiemannData(a, D, r, M, dt, T)
%% Setup
    format("long");
    A = 0; B = 10;
    dx = (B-A)./M;
    x = (A:dx:B)'; t = 0;
    uprev = initfun(x);
    unew = 0*x; 
    ubar = 0*x;
    I = eye(M+1,M+1);
    %plot(x, initfun(x), 'r-');
    %pause(0.1);

    % Advection step: we solve (ubar(j) - uprev)/dt +
    % a*(uprev(j)-uprev(j-1))/dx = 0 for ubar(j), j=2:M+1. This yields the
    % following matrix:

    AdvectionMatrix = sparse(M+1,M+1);
    for j=2:M+1
        AdvectionMatrix(j, j) = -1;
        AdvectionMatrix(j, j-1) = 1;
    end
    AdvectionMatrix = I + AdvectionMatrix.*(a.*dt./dx);
    % Then we have Dirichlet BC at x(1) = A: ubar(1) = 0.
    AdvectionMatrix(1, 1) = 0; 
    
    % Diffusion step: we have (unew(j) - ubar(j))/dt + D*(2*unew(j) -
    % unew(j-1) - unew(j+1))/(2*dx) + r*unew(j) = 0 for j=2:M. At x(1) = A,
    % we have unew(1) = 0 by Dirichlet BC. But, at x(M+1) = B, we have an
    % outflow boundary. So we impose BC u_x = 0 at x(M+1) and approximate 
    % u_x with a one-sided approximation u(M+2) - u(M+1)/dx = 0 yielding
    % u(M+2) = u(M+1). This changes the diffusion term to be
    % D*(unew(j)-unew(j-1))/(2*dx) for j = M+1.
    
    % Reaction term yields this matrix:
    ReactionMatrix = r.*dt.*I;
    ReactionMatrix(1, 1) = 0;

    % Diffusion term is the same for interior nodes
    DiffusionMatrix = sparse(M+1,M+1);
    for j=2:M
        DiffusionMatrix(j, j) = 2;
        DiffusionMatrix(j, j-1) = -1;
        DiffusionMatrix(j, j+1) = -1;
    end
    % At outflow boundary, we have outflow BC
    DiffusionMatrix(M+1,M+1) = 1;
    DiffusionMatrix(M+1,M) = -1;
    DiffusionMatrix = DiffusionMatrix.*((dt.*D)./(dx.^2));
    % At inflow boundary, we have Dirichlet BC
    DiffusionMatrix(1, 1) = 0;

%% Time stepping loop
    while t < T
        t = t + dt;
        % Advection step: solving for ubar
        ubar = AdvectionMatrix*uprev;
        % plot(x, ubar, 'r-');
        % title(sprintf('IMEX numerical solution to D-A-R equation')); axis([0 10 0 1.2]);
        % pause(0.1);
        % Diffusion step: solving for unew
        unew = (I + DiffusionMatrix + ReactionMatrix) \ ubar;
        unew(1) = 0;
        plot(x, unew, 'r-');
        title(sprintf('IMEX numerical solution to D-A-R equation, t = %g', t)); axis([0 10 0 1.2]);
        pause(0.1);
        uprev = unew;
    end
end

%% Needed functions

function H = MyHeaviside(x)
    H=0*x; H(x>0)=1;
    for j=1:length(x)        
        if x(j)<0
            H(j)=0;
        else
            H(j)=1;
        end
    end
end

function v = initfun(x)
    v = MyHeaviside(x - 4) - MyHeaviside(x - 6);
end