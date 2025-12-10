%Import hex cell mesh function 

%Version = 0.1
%Updated = 10-12-25

%Max Wood 2025
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [mesh,cells_full] = import_hex_cell_mesh(filename)

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
end 