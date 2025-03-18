% Sample data
x = randn(100, 1);  % x-coordinate
y = randn(100, 1);  % y-coordinate
z = randn(100, 1);  % z-coordinate
v = randn(100, 1);  % value

% Remove NaN values
valid_indices = ~isnan(x) & ~isnan(y) & ~isnan(z) & ~isnan(v);
x = x(valid_indices);
y = y(valid_indices);
z = z(valid_indices);
v = v(valid_indices);

% Create a grid
[X, Y] = meshgrid(linspace(min(x), max(x), 100), ...
                  linspace(min(y), max(y), 100));

% Interpolate the scattered data onto the grid
Z = griddata(x, y, z, v, X, Y);

% Plot 3D surface plot
surf(X, Y, Z);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Surface Plot');
