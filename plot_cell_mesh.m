%plot cell based mesh 
clc
clearvars


filename = 'grid_cell';

%Load mesh
fid = fopen(filename);
meshi = textscan(fid,'%s %s %s %s');
mesh1 = meshi{1};
mesh2 = meshi{2};
mesh3 = meshi{3};
mesh4 = meshi{4};
fclose(fid);

%read cells 
for ii=1:length(mesh1)
    if strcmp(mesh1(ii),'ncell')
        ncell = str2double(mesh3(ii));
        cells = cell(ncell,1);
        rowc = ii;
        for cc=1:ncell
            rowc = rowc + 1;
            cells{cc}.nedge = str2double(mesh2(rowc));
            cells{cc}.edges = zeros(cells{cc}.nedge,4);
            for ee=1:cells{cc}.nedge
                rowc = rowc + 1;
                cells{cc}.edges(ee,1) = str2double(mesh1(rowc));
                cells{cc}.edges(ee,2) = str2double(mesh2(rowc));
                cells{cc}.edges(ee,3) = str2double(mesh3(rowc));
                cells{cc}.edges(ee,4) = str2double(mesh4(rowc));
            end
        end 
        break
    end 
end

%read vertices
for ii=1:length(mesh1)
    if strcmp(mesh1(ii),'nvertex')
        nvertex = str2double(mesh3(ii));
        vertices = zeros(nvertex,2);
        rowc = ii;
        for vv=1:nvertex
            rowc = rowc + 1;
            vertices(vv,1) = str2double(mesh1(rowc));
            vertices(vv,2) = str2double(mesh2(rowc));
        end
    end
end

%build mesh
mesh.ncell = ncell;
mesh.nvertex = nvertex;
mesh.cells = cells;
mesh.vertices = vertices;


%build full cells
cell_maxnvtx = 0;
for cc=1:mesh.ncell 
    if mesh.cells{cc}.nedge > cell_maxnvtx
        cell_maxnvtx = mesh.cells{cc}.nedge;
    end
end
cells_full = zeros(mesh.ncell,cell_maxnvtx);
for cc=1:mesh.ncell 
    cells_full(cc,1:mesh.cells{cc}.nedge) = mesh.cells{cc}.edges(:,1);
end
cells_full(cells_full == 0) = nan;




% %get cell volumes 
% for cc=1:mesh.ncell 
% 
% 
% end
% 

%initialise
cla reset
hold on


% %plot edge normals 
% efins = 0;
% evins = 0;
% enface = zeros(mesh.ncell,2);
% envtx = zeros(2*mesh.ncell,2);
% for cc=1:mesh.ncell 
%     for ee=1:cells{cc}.nedge
%         efins = efins + 1;
%         evins = evins + 1;
%         enface(efins,1) = evins;
%         v1 = mesh.cells{cc}.edges(ee,1);
%         v2 = mesh.cells{cc}.edges(ee,2);
%         emidv =  0.5*(mesh.vertices(v1,:) + mesh.vertices(v2,:));
%         envtx(evins,:) = emidv;
%         dx = mesh.vertices(v2,1) - mesh.vertices(v1,1);
%         dy = mesh.vertices(v2,2) - mesh.vertices(v1,2);
%         evins = evins + 1;
%         enface(efins,2) = evins;
%         envtx(evins,:) = emidv;
%         envtx(evins,1) = envtx(evins,1) + 0.125*dy;
%         envtx(evins,2) = envtx(evins,2) - 0.125*dx;
%     end
% end 
% enface(enface == 0) = nan;
% patch('faces',enface,'vertices',envtx,'edgecolor','r')




%plot mesh 
patch('Faces',cells_full,'Vertices',mesh.vertices,'edgecolor','k','facecolor','none')

plot(mesh.vertices(:,1),mesh.vertices(:,2),'r.','markersize',10)

% patch('Faces',cells_full(4249,:),'Vertices',mesh.vertices,'edgecolor','k','facecolor','r')

%format
hold off
axis equal
axis tight
box on 

% axis([0.8125    0.8127    0.0599    0.0601])