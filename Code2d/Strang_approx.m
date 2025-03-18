
%% editing...

clear all
close all
format short

% Solve the diffusion-reaction equations by Strang splitting methods

%% Parameters
ax = 0;
bx = 1;

nx= 99;   %Number of interior grid points in space(x)
h=1/(nx+1); %Width of space step(x)
x=0:h:1; %Range of x(0,1)and specifying the grid points
x = transpose(x);
dim = 1;
%% Setting up for time
t = 0; % starting time
k = 0.1; % time step size(t)
T = 2; % final time
k_ori = k;
numk = 5;


%% Initial data
% a = 1;
% u0 = 1 + cos(pi*x); % converge to 1, k = 0.1, T = 2, converge to 1 faster
% u0= -3 +  5*cos(pi*x); % converge to -1, k = 0.025, T = 2,
% u0=  cos(pi*x); % converge to 0, k = 0.1, T = 2,

% a = 0.1;
% u0 = 1 + cos(pi*x); % converge to 1, k = 0.2, T = 2, converge to 1 slower
% u0= -3 +  5*cos(pi*x); % converge to -1, k = 0.04, T = 2,
% u0 = cos(pi*x); % converge to 1, k = 0.2, T = 100,

a = 0.01;
% u0 = 1 + cos(pi*x); % converge to 1, k = 0.4, T = 2, converge to 1 even slower
u0= -3 +  5*cos(pi*x); % break, k = 0.1, T = 2,
% u0 = cos(pi*x); % break, k = 0.4, T = 20,

u0_ori = u0;

%%  Compute the matrices

A = a*matrixA(nx,1);


%% Deefine the solution
v1 = zeros(length(x),1); % solution of diffusion equation
w1 = zeros(length(x),1); % solution of reaction equation
u1 = zeros(length(x),1); % numerical solution

%% Compute the reference solution by ode45

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_ref,u_ref_array] = ode45( @(t,y) myode(t,y,A,1) , [0 T], u0_ori(2:end-1), opts);

u_ref = u_ref_array(end,:);
u_ref = [ (4/3)*u_ref(1)-(1/3)*u_ref(1)  u_ref  (4/3)*u_ref(end)-(1/3)*u_ref(end-1)];
u_ref = u_ref.';

%%  Run simulation for Strang splitting method with correct BCs

for i = 1:4
    
    u0 = u0_ori(2:end-1); % turn back the original inital data for the next step
    nt = ceil(T/k); % number time steps
    
    for j= 1:nt
        
        tau = k/2;
        
        %% Compute the solution of v-equation v_t = v_xx
        
        v1 = expm(tau*A)*u0;
                
        %% Compute the solution for the advection equation w_t = w_x
        
       
%         opts1 = odeset('RelTol',1e-6,'AbsTol',1e-6);
%         [t_reaction,w_reaction] = ode23s( @(t,y) myode(t,y,A,2) , [0 k], v1, opts1);
%         w1 = (w_reaction(end,:)).';
        
%         w1 = sqrt(1./(1+ (1./(v1.^2)-1)*exp(-2*k)) );
%         w1 = -sqrt(1./(1+ (1./(v1.^2)-1)*exp(-2*k)) );
        
        step_temp = k/100;
        N_tilda = ceil(k/step_temp);
        step = k/N_tilda;
        t_rk4 = 0:step:k;
        w1= RK4_Strang(t_rk4,v1);
        
        %% Compute again the numerical solution of v equation in haft step
        
        u1 = expm(tau*A)*w1;
        u0 = u1;
        
    end
    
    u1 = [ (4/3)*u1(1)-(1/3)*u1(1); u1 ;(4/3)*u1(end)-(1/3)*u1(end-1)];
    temp = abs(u1-u_ref);
    
    
    e2 = (h)^0.5* norm(temp, 2)
    
    vec_err2(i) = e2;
    
    vectork(i) = k;
    k = k/2;
    
end
%% Plot

figure
plot(x, u0_ori)
hold on
plot(x, u_ref)
plot(x,u1)
legend('u\_ori','u\_ref','u\_num');
xlabel('x');
ylabel('Value');
% title('');

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
xlim([0.9*min(vectork)   1.1*max(vectork)])
legend('L2 norm','slope 2','Location','southeast');
%% Create a video
figure
myVideo = VideoWriter('myfile.avi');
uncompressedVideo = VideoWriter('myfile.avi', 'Uncompressed AVI');
myVideo.FrameRate = 15;
open(myVideo);
u_min = min(u_ref_array(:));
u_max = max(u_ref_array(:));

time_step_video = floor(length(t_ref)/200);

for j = 1:time_step_video:length(t_ref)
    time = t_ref(j);
    u_ref_video = u_ref_array(j,:);
    u_ref_video = [ (4/3)*u_ref_video(1)-(1/3)*u_ref_video(1)  u_ref_video  (4/3)*u_ref_video(end)-(1/3)*u_ref_video(end-1)];
    u_ref_video = transpose(u_ref_video);
    
    plot(x, u_ref_video, 'r*- ' )
    grid on
    axis ( [ ax, bx, u_min-1.0, u_max+1.0 ] );
    title ( sprintf ( ' Time %f\n', time ) )
    xlabel ( '<-- X -->' )
    ylabel ( 'Value' )
    frame=getframe(gcf);
    writeVideo(myVideo,frame);
end

plot(x, u_ref, 'r*- ' )
grid on
axis ( [ ax, bx, u_min-1.0, u_max+1.0 ] );
title ( sprintf ( ' Time %f\n', T ) )
xlabel ( '<-- X -->' )
ylabel ( 'Value' )
frame=getframe(gcf);
writeVideo(myVideo,frame);

close(myVideo);




