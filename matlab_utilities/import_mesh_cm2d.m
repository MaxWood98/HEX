%Function to import cell_mesh2d volume mesh
%Max Wood - mw16116@bristol.ac.uk
%Univeristy of Bristol - Department of Aerospace Engineering

%Version 1.2
%Updated 22-03-2024

%Function -----------------------------------------------------------------
function [Ncell,Nedge,Nvtx,edge,vtx,cell_lr] = import_mesh_cm2d(filename)

    %Load mesh
    fid = fopen(filename);
    mesh = textscan(fid,'%f %f %f %f');
    fclose(fid);
    
    %Extract mesh
    mesh1 = mesh{1};
    mesh2 = mesh{2};
    mesh3 = mesh{3};
    mesh4 = mesh{4};
    
    %Extract quantities
    Ncell = mesh1(1);
    Nedge = mesh2(1);
    Nvtx = mesh3(1);
    
    %Extract mesh edges
    edge = zeros(Nedge,2);
    cell_lr = zeros(Nedge,2);
    edge(:,1) = mesh1(2:Nedge+1);
    edge(:,2) = mesh2(2:Nedge+1);
    cell_lr(:,1) = mesh3(2:Nedge+1); %left
    cell_lr(:,2) = mesh4(2:Nedge+1); %right
    vtx = zeros(Nvtx,2);
    vtx(:,1) = mesh2(Nedge+2:Nedge+Nvtx+1);
    vtx(:,2) = mesh3(Nedge+2:Nedge+Nvtx+1);
end