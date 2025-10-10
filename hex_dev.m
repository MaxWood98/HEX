%setup 
clearvars
addpath(genpath('matlab_utilities'));

%check cell LR for objects that dont intersect the mesh 
%assign edge external status from surface intersections not from far field flood


%run
% system('hex');

%load geometry
% filename = 'test_geometries/shopt_0p1.fv';
% filename = 'test_geometries/klunkerharder_0p100.fv';
% filename = 'test_geometries/wedge.fv';
% filename = 'test_geometries/NACA_0012.fv';
% filename = 'test_geometries/RAE_2822.fv';
filename = 'test_geometries\cavos_3x3.fv';
[geom] = import_fv_mesh(filename);
% [geom] = build_full_faces(geom);

%load mesh 
[mesh.ncell,mesh.nedge,mesh.nvtx,mesh.edge,mesh.vertex,mesh.cell_lr] = import_mesh_cm2d('CFD_testing/grid');
% [mesh.ncell,mesh.nedge,mesh.nvtx,mesh.edge,mesh.vertex,mesh.cell_lr] = import_mesh_cm2d('naca0012');



%initialise
cla reset
hold on

%plot tree
% vertices = load('vertices');
% faces = load('faces');
% patch('Faces',faces,'Vertices',vertices,'facecolor',[0.5,0.5,0.6],'edgealpha',0.1)

%plot mesh 
patch('Faces',mesh.edge,'Vertices',mesh.vertex,'edgecolor','k')


% %plot geometry
% patch('Faces',geom.edges,'Vertices',geom.vertices,'edgecolor','r','marker','none','edgealpha',0.25)
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



% plot(mesh.vertex(:,1),mesh.vertex(:,2),'k.')

surface_edges = [];
for ii=1:mesh.nedge
    if mesh.cell_lr(ii,1) == -1
        surface_edges = [surface_edges;mesh.edge(ii,:)];
    end
end
patch('Faces',surface_edges,'Vertices',mesh.vertex,'edgecolor','r','linewidth',0.1,'marker','.','markeredgecolor','r')

farfield_edges = [];
for ii=1:mesh.nedge
    if mesh.cell_lr(ii,1) == -2 %|| mesh.cell_lr(ii,2) == -2
        farfield_edges = [farfield_edges;mesh.edge(ii,:)];
    end
end
patch('Faces',farfield_edges,'Vertices',mesh.vertex,'edgecolor','b','linewidth',0.1,'marker','.','markeredgecolor','b')


% farfield_edges = [];
% for ii=1:mesh.nedge
%     if mesh.cell_lr(ii,1) == -3 %|| mesh.cell_lr(ii,2) == -2
%         farfield_edges = [farfield_edges;mesh.edge(ii,:)];
%     end
% end
% patch('Faces',farfield_edges,'Vertices',mesh.vertex,'edgecolor','g','linewidth',0.1,'marker','.','markeredgecolor','g')
% 
% farfield_edges = [];
% for ii=1:mesh.nedge
%     if mesh.cell_lr(ii,1) == -4 %|| mesh.cell_lr(ii,2) == -2
%         farfield_edges = [farfield_edges;mesh.edge(ii,:)];
%     end
% end
% patch('Faces',farfield_edges,'Vertices',mesh.vertex,'edgecolor','m','linewidth',0.1,'marker','.','markeredgecolor','m')



% %edge test
% patch('Faces',mesh.edge(6452,:),'Vertices',mesh.vertex,'edgecolor','r','marker','.','markeredgecolor','r')
% 
% patch('Faces',mesh.edge(5769,:),'Vertices',mesh.vertex,'edgecolor','g','marker','.','markeredgecolor','g')



% patch('Faces',mesh.edge(10029,:),'Vertices',mesh.vertex,'edgecolor','g','marker','.','markeredgecolor','g')
% 
% patch('Faces',mesh.edge(9963,:),'Vertices',mesh.vertex,'edgecolor','m','marker','.','markeredgecolor','m')


    


% mesh.edge(3436,:)
% mesh.cell_lr(3436,:)
% disp('----------')
% mesh.edge(3694,:)
% mesh.cell_lr(3694,:)

% vm = 22189;
% plot(mesh.vertex(vm,1),mesh.vertex(vm,2),'ro','MarkerSize',20)
% 
% vg1 = 88;
% vg2 = 87;
% plot(geom.vertices(vg1,1),geom.vertices(vg1,2),'g.','MarkerSize',20)
% plot(geom.vertices(vg2,1),geom.vertices(vg2,2),'r.','MarkerSize',20)




% mesh.vertex(mesh.edge(3694,:),:)

%format
hold off
axis equal
% axis tight
box on 
xlabel('x')
ylabel('y')

% axis([ -0.2 1.2 -0.6001 0.6001])
% axis([ -1 1 -1 1])
% axis([ -0.6 0.6 -0.6 0.6])

% axis([0.9768    1.0221   -0.0299    0.0089])

%plot edge normals 
% evins = 0;
% enface = zeros(mesh.nedge,2);
% envtx = zeros(2*mesh.nedge,2);
% for ii=1:mesh.nedge
%     evins = evins + 1;
%     enface(ii,1) = evins;
%     emidv =  0.5*(mesh.vertex(mesh.edge(ii,1),:) + mesh.vertex(mesh.edge(ii,2),:));
%     envtx(evins,:) = emidv;
%     dx = mesh.vertex(mesh.edge(ii,2),1) - mesh.vertex(mesh.edge(ii,1),1);
%     dy = mesh.vertex(mesh.edge(ii,2),2) - mesh.vertex(mesh.edge(ii,1),2);
%     evins = evins + 1;
%     enface(ii,2) = evins;
%     envtx(evins,:) = emidv;
%     envtx(evins,1) = envtx(evins,1) - 0.5*dy;
%     envtx(evins,2) = envtx(evins,2) + 0.5*dx;
% end 
% patch('faces',enface,'vertices',envtx,'edgecolor','r')

% hold on
% cell_midpoint = load('cell_midpoint');
% plot(cell_midpoint(:,1),cell_midpoint(:,2),'r.')
% hold off


% %cell associations 
% hold on
% Nedge = mesh.nedge;
% Ncell = mesh.ncell;
% cell_lr = mesh.cell_lr;
% vtx = mesh.vertex;
% edge = mesh.edge;
% cmid = zeros(Ncell,2);
% cmidls = zeros(Ncell,1);
% for ii=1:Nedge
%     cl = cell_lr(ii,1);
%     cr = cell_lr(ii,2);
%     v1 = edge(ii,1);
%     v2 = edge(ii,2);
%     emidx = 0.5*(vtx(v1,1) + vtx(v2,1));
%     emidy = 0.5*(vtx(v1,2) + vtx(v2,2));
%     dx = vtx(v2,1) - vtx(v1,1);
%     dy = vtx(v2,2) - vtx(v1,2);
%     ledge = sqrt(dx^2 + dy^2);
%     if cl > 0
%         cmid(cl,1) = cmid(cl,1) + emidx*ledge;
%         cmid(cl,2) = cmid(cl,2) + emidy*ledge;
%         cmidls(cl) = cmidls(cl) + ledge;
%     end
%     if cr > 0
%         cmid(cr,1) = cmid(cr,1) + emidx*ledge;
%         cmid(cr,2) = cmid(cr,2) + emidy*ledge;
%         cmidls(cr) = cmidls(cr) + ledge;
%     end
% end
% cmid(:,:) = cmid(:,:)./cmidls(:);
% for ii=1:Nedge
%     etgt = ii;
%     v1 = edge(etgt,1);
%     v2 = edge(etgt,2);
%     emidx = 0.5*(vtx(v1,1) + vtx(v2,1));
%     emidy = 0.5*(vtx(v1,2) + vtx(v2,2));
%     cadj = cell_lr(etgt,2);
% 
%     ledge = norm(vtx(v1,:) - vtx(v2,:));
% 
%     % if cadj > 0 && cell_lr(etgt,1) == -1
%     %     plot([emidx cmid(cadj,1)],[emidy cmid(cadj,2)],'b','linewidth',2)
%     % else
%     %     % % plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
%     %     % 
%     %     % if cadj == -1
%     %     % % if cell_lr(etgt,2) > 0
%     %     %     plot([emidx cmid(cell_lr(etgt,1),1)],[emidy cmid(cell_lr(etgt,1),2)],'b','linewidth',2)
%     %     % end 
%     % end
% 
%     % if cell_lr(etgt,1) == 2764 || cell_lr(etgt,2) == 2764
%     %     plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
%     %     plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r.','markersize',20)
%     % end
% 
%     if cell_lr(etgt,1) == 4251 || cell_lr(etgt,2) == 4251
%         plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
%         plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r.','markersize',20)
%     end
% end
% hold off 


% axis([-0.0819    0.0894    0.4314    0.5841])

% axis([0.1459    0.6447    0.1800    0.6060])

% axis([ -0.0115    0.6524   -0.2887    0.2783])
% axis([-0.1435    0.1662   -0.1435    0.1210])
% axis([ -0.0276    0.0466   -0.0404    0.0229])
% axis([-0.0110    0.0271   -0.5257   -0.4932])

% axis([0.0795    0.1248    0.0251    0.0639])

% %plot gradients
% % gradient_vol = load('CFD_testing/gradient.dat');
% gradient_surft = load('CFD_testing/surface_gradient.dat');
% gradient_surf = load('CFD_testing/surface_gradient_mesh.dat');
% hold on 
% % quiver(mesh.vertex(:,1),mesh.vertex(:,2),gradient_vol(:,1),gradient_vol(:,2),0,'b','linewidth',0.5,'AutoScale','off');
% quiver(geom.vertices(:,1),geom.vertices(:,2),gradient_surft(:,1),gradient_surft(:,2),0,'b','linewidth',0.5,'AutoScale','off');
% quiver(geom.vertices(:,1),geom.vertices(:,2),gradient_surf(:,1),gradient_surf(:,2),0,'r','linewidth',0.5,'AutoScale','off');
% hold off


% etgts = [18408,33950,34699,34665];
% etgts = [18380
% 33941
% 18382
% 18381
% 34683
% 34673
% 34674
% 34675
% 34676
% 34677];
% 
% etgts = [18380];
% 
% patch('Faces',mesh.edge(etgts,:),'Vertices',mesh.vertex,'edgecolor','g','linewidth',3)
% 
% axis([0.9952    1.0042   -0.0039    0.0038])