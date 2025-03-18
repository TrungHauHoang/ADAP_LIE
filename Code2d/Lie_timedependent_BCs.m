%% editing...

clear all
close all
format short

% Solve the diffusion-reaction equations by Strang splitting methods

%% Parameters
ax = -0.5;
bx = 1;
nx= 99;   %Number of interior grid points in space(x)
h=(bx-ax)/(nx+1); %Width of space step(x)
x=ax:h:bx; %Range of x(0,1)and specifying the grid points
y=ax:h:bx; 
x = transpose(x);
y = transpose(y);
%% Setting up for time
t = 0; % starting time
% k = 0.00125/2; % time step size(t)
k = 0.05;
T = 0.1; % final time
k_ori = k;
numk = 4;
M = 1;
%% Initial data
mesh_x = x;
mesh_x(1) = [];
mesh_x(end) = [];
mesh_y = y;
mesh_y(1) = [];
mesh_y(end) = [];
[X,Y] = meshgrid(mesh_x,mesh_y);
X = X.';
Y = Y.';


a = 0.01;
u0 = arrayfun(@function_U0,X,Y,3*ones(length(mesh_x)));


u0 = u0(:);
u0_ori = u0;

%%  Compute the matrices
A = a*matrixA(nx,1);
null_matrix = sparse(size(A,1),size(A,1));
x_repeated = repmat(mesh_x,length(mesh_x),1);
y_repeated = Y(:);

theta = min(1./(abs(x_repeated)),M);

% C = log(abs(x_repeated)).*lap2d_nabla(nx,1);
C = (1./(abs(x_repeated))+1./(abs(y_repeated))).*lap2d_nabla(nx,1);
D = (1./(abs(x_repeated))+1./(abs(y_repeated)) - theta).*lap2d_nabla(nx,1);
E = (theta).*lap2d_nabla(nx,1);
I = eye(size(A));

%% Deefine the solution
v1 = zeros(length(x),1); % solution of diffusion equation
w1 = zeros(length(x),1); % solution of reaction equation
u1 = zeros(length(x),1); % numerical solution

%% Compute the reference solution by ode45
tStart = tic;
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_ref,u_ref_array] = ode45( @(t,y) myode_time(t,y,A,C,6,X(:),Y(:),a,nx,theta) , [0 T], u0_ori, opts);
u_ref = (u_ref_array(end,:)).';
tEnd = toc(tStart)

%%  Run simulation for Lie splitting method with correct BCs

for i = 1:5
    u0 = u0_ori; % turn back the original inital data for the next step
    nt = ceil(T/k); % number time steps
    tt = 0;
    for j= 1:nt
        
        ttend = tt + k;
        %% Compute the solution for the advection equation w_t = w_x
%         tStart = tic;
        [t_ref,u_adv_array] = ode45( @(t,y) myode_time(t,y,null_matrix,C,7,X(:),Y(:),a,nx,theta), [tt ttend], u0, opts);
        u_adv = (u_adv_array(end,:)).';
%         tEnd = toc(tStart)
        %% Compute again the numerical solution of v equation in haft step
        [t_ref,u_dif_array] = ode45( @(t,y) myode_time(t,y,A,null_matrix,8,X(:),Y(:),a,nx,theta), [tt ttend], u_adv, opts);
        u1 = (u_dif_array(end,:)).';
        tt = ttend;
        u0 = u1;
        
    end
    
    temp = abs(u1-u_ref);
    e1 = h*norm(temp,2)
    vec_err2(i) = e1;
    vectork(i) = k;
    k = k/2;
    
end

k = k_ori;
for i = 1:5
    
    u0 = u0_ori; % turn back the original inital data for the next step
    nt = ceil(T/k); % number time steps
    tt = 0;
    for j= 1:nt
        %% Compute the solution for the advection equation w_t = w_x
        ttend = tt + k;
        [t_ref,u_adv_array_mod] = ode45( @(t,y) myode_time(t,y,null_matrix,E,9,X(:),Y(:),a,nx,theta) , [tt ttend], u0, opts);
        u_adv_mod = (u_adv_array_mod(end,:)).';
        
        %% Compute again the numerical solution of v equation in haft step
        [t_ref,u_dif_array_mod] = ode45( @(t,y) myode_time(t,y,A,D,10,X(:),Y(:),a,nx,theta) , [tt ttend], u_adv_mod, opts);
        u1 = (u_dif_array_mod(end,:)).';
        tt = ttend;
        u0 = u1;
        
    end
    
    temp = abs(u1-u_ref);
    e2 = h*norm(temp,2)
    vec_err2_mod(i) = e2;
    vectork(i) = k;
    k = k/2;
    
end




%% Plot the solution
u_min = min(u_ref_array(end,:));
u_max = max(u_ref_array(end,:));
matrix_u0 =  reshape(u0_ori,[nx,nx]);
matrix_u_Strang = reshape(u1,[nx,nx]);
matrix_u_ref = reshape(u_ref,[nx,nx]);
figure
s = surfl(X,Y,matrix_u0);
legend('u\_ori');
s.MeshStyle = 'none';
xlabel ( '<-- X -->' )
ylabel ( '<-- Y -->' )
zlabel ( '<-- Value -->' )

h_plot = figure;
s = pcolor(X,Y,matrix_u_ref);
xlabel ( '<-- X -->' )
ylabel ( '<-- Y -->' )
% str=sprintf('T = %.f  ',T);
% title(str);
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';
% 
% 
%% the L2 error
h_plot = figure;
loglog( vectork, vec_err2,'-ro',...
    'LineWidth',1,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',5);
hold on
loglog( vectork, vec_err2_mod,'-mx',...
    'LineWidth',1,...
    'MarkerEdgeColor','m',...
    'MarkerFaceColor','m',...
    'MarkerSize',5);

loglog(vectork, vectork,'-.k');
xlabel('\tau');
ylabel('L_2 Error');
% title('Errors');
xlim([min(vectork)   max(vectork)])
legend('Classical Lie','Adapted Lie','Slope 1','Location','southeast');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
s.MeshStyle = 'none';
s.FaceColor = 'interp';

% 
% 
% 











