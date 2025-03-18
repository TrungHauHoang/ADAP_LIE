% Define the mesh
mesh = [-1, -0.5, -0.1, 0.1, 0.5, 1]; % Example non-equidistant mesh with 6 points

% Construct the finite difference matrix for uxx
N = length(mesh);
D = zeros(N, N);

% Interior points
for i = 2:N-1
    h_left = mesh(i) - mesh(i-1);
    h_right = mesh(i+1) - mesh(i);
    D(i, i-1) = 1 / (h_left * (h_left + h_right));
    D(i, i) = -2 / (h_left * h_right);
    D(i, i+1) = 1 / (h_right * (h_left + h_right));
end

% Boundary points (assuming Dirichlet boundary conditions)
% You need to specify boundary conditions for the first and last points
% Here, I'm assuming Dirichlet boundary conditions where u' = 0 at the boundaries
D(1, 1) = -2 / (mesh(2) - mesh(1));
D(1, 2) = 1 / (mesh(2) - mesh(1));
D(N, N-1) = 1 / (mesh(N) - mesh(N-1));
D(N, N) = -2 / (mesh(N) - mesh(N-1));

disp(D);
