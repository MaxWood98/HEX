%setup 
clearvars
addpath(genpath('matlab_utilities'));

%load geometry
[geom] = import_fv_mesh('test_geometries/NACA_0012.fv');
[geom] = build_full_faces(geom);

%load mesh 
[mesh,cells_full] = import_hex_cell_mesh('grid_cell');
% [mesh.ncell,mesh.nedge,mesh.nvtx,mesh.edge,mesh.vertex,mesh.cell_lr] = import_mesh_cm2d('grid');
% [mesh.ncell,mesh.nedge,mesh.nvtx,mesh.edge,mesh.vertex,mesh.cell_lr] = import_mesh_cm2d('naca0012');


%initialise
cla reset
hold on

%plot geometry
patch('Faces',geom.edges,'Vertices',geom.vertices,'edgecolor','r','marker','none','edgealpha',0.25)
% normals = zeros(length(geom.edges),2);
% midpoints = zeros(length(geom.edges),2);
% for ii=1:length(geom.edges)
%     v1 = geom.edges(ii,1);
%     v2 = geom.edges(ii,2);
%     dx = geom.vertices(v2,1) - geom.vertices(v1,1);
%     dy = geom.vertices(v2,2) - geom.vertices(v1,2);
%     normals(ii,1) = dy;
%     normals(ii,2) = -dx;
%     midpoints(ii,1) = 0.5*(geom.vertices(v2,1) + geom.vertices(v1,1));
%     midpoints(ii,2) = 0.5*(geom.vertices(v2,2) + geom.vertices(v1,2));
% end
% quiver(midpoints(:,1),midpoints(:,2),normals(:,1),normals(:,2))


%plot gradients
% gradient_vol = load('CFD_testing/gradient');
gradient_surf = load('surface_gradient');
hold on 
% quiver(mesh.vertices(:,1),mesh.vertices(:,2),gradient_vol(:,1),gradient_vol(:,2),0,'b','linewidth',0.5,'AutoScale','off');
quiver(geom.vertices(:,1),geom.vertices(:,2),gradient_surf(:,1),gradient_surf(:,2),0,'r','linewidth',0.5,'AutoScale','off');
hold off

%format 
axis equal 
box on 

