%Import Face-Vertex (.fv) mesh Function 

%Version = 0.1
%Updated = 18-12-24

%Max Wood 2024
%University of Bristol
%Department of Aerospace Engineering

%Function 
function [mesh] = import_fv_mesh(filename)

    %Open
    fid = fopen(filename);
    
    %Read contents as text
    filedatain = textscan(fid,'%s','Delimiter','','MultipleDelimsAsOne',1); 
    filedata = filedatain{1,1};
    
    %Close
    fclose(fid);
    
    %Set length 
    filelen = size(filedata,1);
    
    %Initialise structure
    mesh.ndim = 0;
    mesh.nvertex = 0;
    mesh.nedge = 0;
    mesh.nface = 0;
    mesh.max_vtx_inface = 0;
    mesh.vertices = [];
    mesh.edges = [];
    mesh.faces = [];
    
    %Find and read number of dimensions
    for ii=1:filelen
        if length(filedata{ii,1}) >= 4
            if strcmp(filedata{ii,1}(1:4),'ndim')
                line = filedata{ii,1};
                lines = split(line,'=');
                mesh.ndim = str2double(lines{2});
                break
            end
        end
    end
    
    %Find and read vertices
    for ii=1:filelen
        if length(filedata{ii,1}) >= 7
            if strcmp(filedata{ii,1}(1:7),'nvertex')
        
                %Read number of vertices
                line = filedata{ii,1};
                lines = split(line,'=');
                mesh.nvertex = str2double(lines{2});
                
                %Read vertices
                mesh.vertices = zeros(mesh.nvertex,mesh.ndim);
                vertexdat = filedata(ii+1:ii+mesh.nvertex,1);
                for jj=1:mesh.nvertex
                    line = vertexdat{jj,1};
                    lines = strsplit(line,' ');
                    mesh.vertices(jj,:) = str2double(lines(1:mesh.ndim));
                end
                break
            end
        end
    end
    
    %Find and read edges
    for ii=1:filelen
        if length(filedata{ii,1}) >= 5
            if strcmp(filedata{ii,1}(1:5),'nedge')
            
                %Read number of edges
                line = filedata{ii,1};
                lines = split(line,'=');
                mesh.nedge = str2double(lines{2});
                
                %Read edges
                mesh.edges = zeros(mesh.nedge,2);
                edgedat = filedata(ii+1:ii+mesh.nedge,1);
                for jj=1:mesh.nedge
                    line = edgedat{jj,1};
                    lines = strsplit(line,' ');
                    mesh.edges(jj,:) = str2double(lines(1:2));
                end
                break
            end
        end
    end
    
    
    %Find and read faces
    if mesh.ndim == 3 

        %Load faces
        for ii=1:filelen
            if length(filedata{ii,1}) >= 5
                if strcmp(filedata{ii,1}(1:5),'nface')
                
                    %Read number of faces
                    line = filedata{ii,1};
                    lines = split(line,'=');
                    mesh.nface = str2double(lines{2});
                    
                    %Read faces
                    mesh.faces = cell(mesh.nface,1);
            
                    facedat = filedata(ii+1:ii+mesh.nface,1);
                    for jj=1:mesh.nface
                        line = facedat{jj,1};
                        lines = strsplit(line,' ');
                        mesh.faces{jj}.nvertex = str2double(lines(1));
                        mesh.faces{jj}.vertices = zeros(mesh.faces{jj}.nvertex,1);
                        mesh.faces{jj}.vertices(:) = str2double(lines(2:mesh.faces{jj}.nvertex+1));
                    end
                    break
                end
            end
        end

        %Find maximum number of vertices in a face
        mesh.max_vtx_inface = 0;
        for ii=1:mesh.nface
            mesh.max_vtx_inface = max([mesh.max_vtx_inface,mesh.faces{ii}.nvertex]);
        end
    end 
end