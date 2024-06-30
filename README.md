# AdvectionDiffusionReactionOperatorSplitting

A MATLAB implementation of the advection-diffusion-reaction equation using the operator splitting scheme.

## Description

This repository contains a MATLAB code that demonstrates the solution of the advection-diffusion-reaction equation using the operator splitting scheme. The initial condition is set using a custom Heaviside function, and the numerical solution is visualized over time.

## Code Overview

The main function `IMEXRiemannData` advances the solution of the advection-diffusion-reaction equation using an implicit-explicit (IMEX) operator splitting scheme. The code includes:
- Initialization of the spatial grid, initial condition, and matrices for advection, diffusion, and reaction terms.
- A time-stepping loop that advances the solution, updates the plot, and visualizes the numerical solution at each time step.

### Main Function

```matlab
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
 
    AdvectionMatrix = sparse(M+1,M+1);
    for j=2:M+1
        AdvectionMatrix(j, j) = -1;
        AdvectionMatrix(j, j-1) = 1;
    end
    AdvectionMatrix = I + AdvectionMatrix.*(a.*dt./dx);
    AdvectionMatrix(1, 1) = 0; 
    
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
        unew = (I + DiffusionMatrix + ReactionMatrix) \ ubar;
        unew(1) = 0;
        plot(x, unew, 'r-');
        title(sprintf('IMEX numerical solution to D-A-R equation, t = %g', t)); axis([0 10 0 1.2]);
        pause(0.1);
        uprev = unew;
    end
end
```
## Helper Functions
```matlab
function H = SFHeaviside(x)
    H = 0*x; 
    H(x > 0) = 1;
    for j = 1:length(x)        
        if x(j) < 0
            H(j) = 0;
        else
            H(j) = 1;
        end
    end
end

```
## Initial condition
```matlab
function v = initfun(x)
    v = MyHeaviside(x - 4) - MyHeaviside(x - 6);
end
```
## Usage
To run the code, call the IMEXRiemannData function with the desired parameters for advection speed a, diffusion coefficient D, reaction rate r, number of grid points M, time step dt, and total simulation time T. For example:
```matlab
a = 1;    % advection speed
D = 0.1;  % diffusion coefficient
r = 0.01; % reaction rate
M = 100;  % number of grid points
dt = 0.01; % time step
T = 2;    % total simulation time
IMEXRiemannData(a, D, r, M, dt, T);

```
## License
This project is licensed under the MIT License - see the LICENSE file for details.
```
Feel free to adjust any part of this README to better fit your specific needs or preferences.
