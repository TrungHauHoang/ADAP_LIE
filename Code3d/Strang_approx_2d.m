%% editing...

clear all
close all
format short

% Solve the diffusion-reaction equations by Strang splitting methods

%% Parameters
ax = 0;
bx = 1;

nx= 39;   %Number of interior grid points in space(x)
h=1/(nx+1); %Width of space step(x)
x=0:h:1; %Range of x(0,1)and specifying the grid points
y=0:h:1;
z=0:h:1;
x = transpose(x);
y = transpose(y);
z = transpose(z);

%% Setting up for time
t = 0; % starting time
k = 0.001; % time step size(t)
T = 0.1; % final time
k_ori = k;
numk = 3;


%% Initial data
mesh_x = x;
mesh_x(1) = [];
mesh_x(end) = [];
mesh_y = y;
mesh_y(1) = [];
mesh_y(end) = [];
mesh_z = z;
mesh_z(1) = [];
mesh_z(end) = [];
[X,Y,Z] = meshgrid(mesh_x,mesh_y,mesh_z);
% X = X.';
% Y = Y.';
% Z = Z.';


a = 1;
% u0 = arrayfun(@function_U0,X,Y,Z,1*ones(length(mesh_x)));
% u0 = u0(:);
% u0_ori = u0;

u0 = zeros(nx*nx*nx,1);
u0_ori = u0;
%%  Compute the matrices

A = a*matrixA(nx,1);
x_repeated = repmat(mesh_x,length(mesh_x),1);
x_repeated3d = repmat(x_repeated,length(mesh_x),1);
C = log(x_repeated3d).*lap2d_nabla(nx,1);
% C = (1./(x_repeated)).*lap2d_nabla(nx,1);
vector_d = 1./(x_repeated3d.^2);
I = speye(size(A));

%% Deefine the solution
% v1 = zeros(length(x),1); % solution of diffusion equation
% w1 = zeros(length(x),1); % solution of reaction equation
% u1 = zeros(length(x),1); % numerical solution

%% Compute the reference solution by ode45
% tStart = tic;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_ref,u_ref_array] = ode45( @(t,y) myode(t,y,vector_d,A,C,1) , [0 T], u0_ori, opts);
u_ref = (u_ref_array(end,:)).';
% tEnd = toc(tStart)


%%  Run simulation for Lie splitting method with correct BCs


for i = 1:5
    
    u0 = u0_ori; % turn back the original inital data for the next step
    nt = ceil(T/k); % number time steps
    extreigs = gersh((k/2)*A);
    dif_min = extreigs.SR;
    dif_max = 0;
    
    for j= 1:nt
        
        tau = k/2;
        
        %% Compute the solution for the advection equation w_t = w_x
        
        W_2nd = vector_d;
        C_tilda = [C W_2nd; zeros(1,size(C,1)) 0];
        V0_2nd = [u0; 1];
        [~,ehAv] = real_leja_sparse(k,C_tilda,dif_min*2,dif_max,V0_2nd/norm(V0_2nd),1e-14,0);
        ehAv = norm(V0_2nd)*ehAv;
        %         ehAv = expm(k*C_tilda)*V0_2nd;
        w1 = [ I zeros(size(C,1),1) ]*ehAv;
        
        %% Compute again the numerical solution of v equation in haft step
        [~,u1] = real_leja_sparse(k,A,dif_min*2,dif_max,w1/norm(w1),1e-14,0);
        u1 = norm(w1)*u1;
%         u1 = expm(k*A)*w1;
        u0 = u1;
        
    end
    
    temp = abs(u1-u_ref);
    e2 = h*norm(temp,2)
    vec_err2(i) = e2;
    vectork(i) = k;
    k = k/2;
    
end

% e2 =
% 
%     0.4345
% 
% 
% e2 =
% 
%     0.2602
% 
% 
% e2 =
% 
%     0.1436
% 
% 
% e2 =
% 
%     0.0762
% 
% 
% e2 =
% 
%     0.0394

% %%  Run simulation for Strang splitting method with correct BCs
%
%
% for i = 1:3
%
%     u0 = u0_ori; % turn back the original inital data for the next step
%     nt = ceil(T/k); % number time steps
%     extreigs = gersh((k/2)*A);
%     dif_min = extreigs.SR;
%     dif_max = 0;
%
%     for j= 1:nt
%
%         tau = k/2;
%
%         %% Compute the solution of v-equation v_t = v_xx
%         if norm(u0) > 1e-7
%             [~,v1] = real_leja_sparse(tau,A,dif_min,dif_max,u0/norm(u0),1e-10,0);
%             v1 = norm(u0)*v1;
%         elseif norm(u0) < 1e-7
%             v1 = expm(tau*A)*u0;
%         end
%
%         %         v1 = expm(tau*A)*u0;
%         %% Compute the solution for the advection equation w_t = w_x
%
%         %         [~,w1] = real_leja_sparse(k,C,dif_min*2,dif_max,v1/norm(v1),1e-10,0);
%         %         w1 = norm(v1)*w1;
%
%         W_2nd = vector_d;
%         C_tilda = [C W_2nd; zeros(1,size(C,1)) 0];
%         V0_2nd = [v1; 1];
%         %% exp(kA) computed by exact formulation
%         ehAv = expm(k*C_tilda)*V0_2nd;
%         %% Compute the solution
%         w1 = [ I zeros(size(C,1),1) ]*ehAv;
%
%
%         %% Compute again the numerical solution of v equation in haft step
%
%         if norm(w1) > 1e-7
%             [~,u1] = real_leja_sparse(tau,A,dif_min,dif_max,w1/norm(w1),1e-10,0);
%             u1 = norm(w1)*u1;
%             u0 = u1;
%         elseif norm(w1) < 1e-7
%             u1 = expm(tau*A)*w1;
%             u0 = u1;
%         end
%
%     end
%
%     temp = abs(u1-u_ref);
%     e2 = h*norm(temp,2)
%     vec_err2(i) = e2;
%     vectork(i) = k;
%     k = k/2;
%
% end
%% Plot the solution
% u_min = min(u_ref_array(end,:));
% u_max = max(u_ref_array(end,:));
% matrix_u0 =  reshape(u0_ori,[nx,nx]);
% matrix_u_Strang = reshape(u1,[nx,nx]);
% matrix_u_ref = reshape(u_ref,[nx,nx]);
% 
% % u0 = u0_ori;
% figure
% s = surfl(X,Y,matrix_u0);
% legend('u\_ori');
% s.MeshStyle = 'none';
% xlabel ( '<-- X -->' )
% ylabel ( '<-- Y -->' )
% zlabel ( '<-- Value -->' )
% 
% 
% figure
% s = surfl(X,Y,matrix_u_Strang);
% legend('u\_num');
% zlim([u_min-1.0 u_max+1.0])
% s.MeshStyle = 'none';
% xlabel ( '<-- X -->' )
% ylabel ( '<-- Y -->' )
% zlabel ( '<-- Value -->' )

%% the L2 error
figure
loglog( vectork, vec_err2,'-ro',...
    'LineWidth',1,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',5);
hold on
loglog(vectork, vectork.^2,'-.k');
xlabel('\tau');
ylabel('L_2 Error');
title('Errors');
xlim([min(vectork)   max(vectork)])
legend('L_2 norm','slope 2','Location','southeast');
% %% Create a video
% figure
% myVideo = VideoWriter('myfile.avi');
% uncompressedVideo = VideoWriter('myfile.avi', 'Uncompressed AVI');
% myVideo.FrameRate = 15;
% open(myVideo);
%
% time_step_video = floor(length(t_ref)/40);
%
% for j = 1:time_step_video:length(t_ref)
%     time = t_ref(j);
%
%     u_ref_video = (u_ref_array(j,:)).';
%     matrix_u_ref = reshape(u_ref_video,[nx,nx]);
%
%     s = surfl(X,Y,matrix_u_ref);
%     legend('u\_ref');
%     %     zlim([u_min-1.0 u_max+1.0])
%     s.MeshStyle = 'none';
%
%     axis ( [ ax, bx,ax ,bx, u_min-1.0, u_max+1.0 ] );
%     title ( sprintf ( ' Time %f\n', time ) )
%     xlabel ( '<-- X -->' )
%     ylabel ( '<-- Y -->' )
%     zlabel ( '<-- Value -->' )
%     legend('u\_ref');
%     frame=getframe(gcf);
%     writeVideo(myVideo,frame);
% end
%
% s = surfl(X,Y,matrix_u_ref);
% s.MeshStyle = 'none';
% axis ( [ ax, bx,ax ,bx, u_min-1.0, u_max+1.0 ] );
% title ( sprintf ( ' Time %f\n', T ) )
% xlabel ( '<-- X -->' )
% ylabel ( '<-- Y -->' )
% zlabel ( '<-- Value -->' )
% legend('u\_ref');
% frame=getframe(gcf);
% writeVideo(myVideo,frame);
%
% close(myVideo);




