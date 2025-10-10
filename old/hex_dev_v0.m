%setup 
clearvars
addpath(genpath('matlab_utilities'))

%load
vertices = load('vertices');
faces = load('faces');
verticesf = load('verticesf');
facesf = load('facesf');

filename = 'test_geometries/NACA_0012.fv';
[mesh] = import_fv_mesh(filename);
[mesh] = build_full_faces(mesh);


%plot
cla reset
% figure
hold on
% patch('Faces',faces,'Vertices',vertices,'facecolor',[0.5,0.5,0.6],'edgealpha',0.1)
patch('Faces',facesf(:,1:2),'Vertices',verticesf,'facecolor','k')
% plot(vertices(:,1),vertices(:,2),'r.','MarkerSize',2)
% patch('Faces',mesh.faces_full,'Vertices',mesh.vertices,'edgecolor','r','marker','.')

% plot(mesh.vertices(:,1),mesh.vertices(:,2),'r.','MarkerSize',10)
% plot(verticesf(:,1),verticesf(:,2),'g.','MarkerSize',10)

hold off
axis equal
axis tight
box on 

axis([ -0.1311    1.2923   -0.6196    0.6196])

% axis equal 
% axis tight

% vdist = zeros(length(vertices),length(vertices));
% for ii=1:length(vertices)
%     for jj=1:length(vertices)
%         vdist(ii,jj) = norm(vertices(ii,:) - vertices(jj,:));
%         if ii ~= jj
%             if vdist(ii,jj) < 0.00000000001
%                 disp('double vtx')
%             end
%         end
%     end
% end 


% Ncell = max(max(cell_lr(:,1)),max(cell_lr(:,2)));
% Nedge = length(edge);

% hold on
% edge = facesf(:,1:2);
% cell_lr = facesf(:,3:4);
% vtx = verticesf;

 
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
%     if cadj > 0
%         plot([emidx cmid(cadj,1)],[emidy cmid(cadj,2)],'b','linewidth',2)
%     else
%         plot([vtx(v1,1) vtx(v2,1)],[vtx(v1,2) vtx(v2,2)],'r','linewidth',2)
% 
%         % if cadj == -2
%         if cell_lr(etgt,1) > 0
%             plot([emidx cmid(cell_lr(etgt,1),1)],[emidy cmid(cell_lr(etgt,1),2)],'b','linewidth',2)
%         end 
%     end
% end
% hold off 


%plot edge normals 
evins = 0;
enface = zeros(2*nedge,2);
envtx = zeros(2*nedge,2);
for ii=1:Nedge
    evins = evins + 1;
    enface(ii,1) = evins;
    emidv =  0.5*(vtx(edge(ii,1),:) + vtx(edge(ii,2),:));
    envtx(evins,:) = emidv;
    dx = vtx(edge(ii,2),1) - vtx(edge(ii,1),1);
    dy = vtx(edge(ii,2),2) - vtx(edge(ii,1),2);
    evins = evins + 1;
    enface(ii,2) = evins;
    envtx(evins,:) = emidv;
    envtx(evins,1) = envtx(evins,1) - 0.5*dy;
    envtx(evins,2) = envtx(evins,2) + 0.5*dx;
end 
patch('faces',enface,'vertices',envtx,'edgecolor','r')