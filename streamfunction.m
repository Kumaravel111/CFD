clc;            % Clear the command window
clear;          % Clear all variables from the workspace
close all;      % Close all figure windows

% Grid dimensions
m = 41;             % Number of grid points along x direction
n = 31;             % Number of grid points along y direction.
dx = 7.0 / (m - 1); % grid size at x direction  
dy = 5.0 / (n - 1); % grid size at y direction

% create 2 matrixs for store the stream function old and new values  
psi_old = zeros(m, n);
psi_new = zeros(m, n);

% create of stream function stream function coefficients S,W,P,E,N
S = 1.0 / dy^2;
W = 1.0 / dx^2;
P = 2.0 * ((1.0 / dx^2) + (1.0 / dy^2));
E = 1.0 / dx^2;
N = 1.0 / dy^2;

% Initialize error value  and iteration count
e = 1.0; % Start with a large error to enter the loop for solving stream function
itr = 0; % For know how much iteration it took 


% Apply boundary conditions
for i = 1:m
    for j = 1:n
        if j == n
            psi_new(i, j) = 0.0; % top boundary
        elseif i == 1
            psi_new(i, j) = 0.0; % left boundary
        elseif i <= (1/dx) 
            psi_new(i, j) = 0.0; % bottom upto x = 1.0 boundary
        elseif i >= (1.3/dx)
            psi_new(i, j) = 100.0; % bottom up to x = 1.3 to 7 boundary
        else
            psi_new(i, j) = 0.0; 
        end
    end
end

% jacobi iteration method for solving stream function equation

while e > 1e-8
    
    psi_old = psi_new; % replace psiold with psi new

    % Compute psi_new for internal grid points 

    for i = 2:(m - 1) 
        for j = 2:(n - 1) 
            psi_new(i, j) = (1.0 / P) * ( ...
                + S * psi_old(i, j - 1) ...
                + W * psi_old(i - 1, j) ...
                + E * psi_old(i + 1, j) ...
                + N * psi_old(i, j + 1));
        end
    end

    % Neumann boundary condition at the right wall (outlet)

    psi_new(m, :) = psi_new(m - 1, :);

    % Calculate the error

    e = 0; % Reset error
    for i = 1:m
        for j = 1:n
            e = e + (psi_new(i, j) - psi_old(i, j))^2;
        end
    end
    e = sqrt(e / (m * n));

    % Increment iteration counter
    itr = itr + 1;
end

%  results

disp(['Converged in ', num2str(itr), ' iterations.']);
disp(['Final Error: ', num2str(e)]);

% Generate grid for plotting
x = linspace(0, 7, m);
y = linspace(0, 5, n);
[X, Y] = meshgrid(x, y);

% Plot psi as streamlines
figure;
contour(X, Y, psi_new', 25, 'LineWidth', 1.2); % Transpose psi_new for correct orientation
colorbar;
caxis([min(psi_new(:)), max(psi_new(:))]);
title('Streamlines of \psi');
xlabel('X-axis');
ylabel('Y-axis');
grid on;