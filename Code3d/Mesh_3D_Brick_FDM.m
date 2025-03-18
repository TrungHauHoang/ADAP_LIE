
% This program initialize the creation of the Cubic element and faces to show info in MATLAB
% directly,
% It is useful for finite difference problems (FDM) to show the
% results directly in the 3D matlab
%% By Amin Zaami, Ph.D.  Hengelo, The Netherlands, 2022
% https://www.linkedin.com/in/aminzaami/
%**** Share and redistribute WITHOUT PERMISSION !!   ***
function  FDM_mesh=Mesh_3D_Brick_FDM (xnode,ynode,znode,Lx,Ly,Lz)
% % Brick elements properties, random discritization
% xnode=50; > discritization x
% ynode=25;  > discritization y
% znode=16;   > discritization z
% Lx=10;  >> length along x-axis
% Ly=2;  >> along y-axis
% Lz=4;   >> along z-axis
%% make elements
counter=0;
N_plane=xnode*ynode;
N=xnode*ynode*znode;
eleNx=xnode;
eleNy=ynode;
eleNz=znode;
% number of element
ele_NUM=(eleNx-1)*(eleNy-1)*(eleNz-1) ;
% 8 node brick element
Element=zeros(ele_NUM,8,'int32');
for kk=1:znode-1
    for ii=1: xnode-1
        for jj=1: ynode-1
            counter=counter+1;
            Element(counter,1:8)=[jj+((kk-1)*N_plane)+(ii-1)*ynode  (jj+1)+((kk-1)*N_plane)+(ii-1)*ynode,  (ii*ynode)+(jj)+((kk-1)*N_plane) (ii*ynode)+(jj+1)+((kk-1)*N_plane),...
                jj+((kk)*N_plane)+(ii-1)*ynode  (jj+1)+((kk)*N_plane)+(ii-1)*ynode,  (ii*ynode)+(jj)+((kk)*N_plane) (ii*ynode)+(jj+1)+((kk)*N_plane) ];
        end
    end
end
%% Face making
face_all=zeros(ele_NUM*6,5,'int32');
for ii=1:ele_NUM
    cur_ele= Element(ii,:);
    face= [cur_ele(1) cur_ele(2) cur_ele(4) cur_ele(3);...
        cur_ele(1) cur_ele(2) cur_ele(6) cur_ele(5);...
        cur_ele(1) cur_ele(3) cur_ele(7) cur_ele(5);...
        cur_ele(8) cur_ele(6) cur_ele(5) cur_ele(7);...
        cur_ele(8) cur_ele(6) cur_ele(2) cur_ele(4);...
        cur_ele(8) cur_ele(4) cur_ele(3) cur_ele(7)] ;
    % to put index of element near faces, we know which face belong to
    % which element
    face_all(6*(ii-1) + (1:6),:)=[ones(6,1,'int32')*ii, face];
end
%% Node numbering 
xnode_P=linspace(0,Lx,eleNx);
ynode_P=linspace(0,Ly,eleNy);
znode_P=linspace(Lz,0,eleNz);
counter=0;
Nodes=zeros(N,3,'single');
for kk=1:znode
    for ii=1: xnode
        for jj=1: ynode
            counter=counter+1;
%             fprintf(fileID6,' %f %f  %f %f \r\n', xnode_P(ii), ynode_P(jj), znode_P(kk), T(counter));
       Nodes(counter,:)=[xnode_P(ii), ynode_P(jj), znode_P(kk)];
        
        
        end
    end
end
% only for initialiaztion purposes
% random data to be filled in
Temp_GLob=zeros(N,1,'single')+log2(1:(N))';
%% plotting
figure;
FDM_mesh=patch('Faces', face_all(:,2:end), 'Vertices', Nodes, 'FaceVertexCData', Temp_GLob, 'FaceColor', 'interp',...
    'EdgeColor',[.85 0.85, 0.85],'LineWidth',0.05,'FaceAlpha',0.5,'EdgeAlpha',0.5 );
view([50 30]);
colorbar;
colormap('jet')
axis equal;
