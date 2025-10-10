%Build Full Faces Function 

%Version = 0.2
%Updated = 01-05-24

%Max Wood 2024
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [mesh] = build_full_faces(mesh)
    if mesh.ndim == 2
        mesh.faces_full = zeros(mesh.nedge,2);
        mesh.faces_full(:,:) = mesh.edges(:,:);
    elseif mesh.ndim == 3
        mesh.faces_full = zeros(mesh.nface,mesh.max_vtx_inface);
        for ii=1:mesh.nface
            mesh.faces_full(ii,1:mesh.faces{ii}.nvertex) = mesh.faces{ii}.vertices(:);
        end
        mesh.faces_full(mesh.faces_full == 0) = nan;
    end
end