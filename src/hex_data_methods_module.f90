!hex data and methods module 
!max wood
!version : 0.0.7
!updated : 25-03-25

!module 
module hex_data_methods

!integer data types 
use ISO_FORTRAN_ENV, only: in16=>int16
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in64=>int64

!real data types 
use ISO_FORTRAN_ENV, only: sp=>real32 
use ISO_FORTRAN_ENV, only: dp=>real64

!options type
type hex_options
    logical :: cdisplay,allow_postprocess,geom_from_commandline,options_from_commandline,simplify_surfaces
    character(len=:), allocatable :: geomname,geompath,optionsname,optionspath,gradientpath,gradientname
    character(len=:), allocatable :: meshpath,meshname,gradient_method
    character(len=:), allocatable :: mesh_in_out,mode,mesh_treetype,tree_mesh_relation,geom_cell_link_method
    integer(in64) :: nrefine,nflood_coarse,nflood_mid,nflood_fine
    integer(in64) :: ncell_max,nsmooth_interlayer,nsmooth_farfield,nrefzone,nclipplane,nbczone,nbcremzone
    integer(in64), dimension(:), allocatable :: boundarycondition_remove_zones
    real(dp) :: farfield_r,edgelength_min,cellvol_min,vtx_proximity_tolerance
    real(dp) :: tree_offset(3)
    type(hex_refzone), dimension(:), allocatable :: refinement_zones
    type(hex_clipplane), dimension(:), allocatable :: clip_planes
    type(hex_bczone), dimension(:), allocatable :: boundarycondition_zones
end type hex_options 

!refinement zone type
type hex_refzone
    character(len=:), allocatable :: type 
    integer(in64) :: rlevel 
    real(dp) :: xmid,ymid,zmid,xmin,ymin,zmin,xmax,ymax,zmax,radius
    real(dp) :: v1(3),v2(3)
end type hex_refzone

!clipplane plane type
type hex_clipplane
    real(dp) :: v1(3),v2(3)
end type hex_clipplane

!boundary condition zone type
type hex_bczone
    integer(in64) :: bctag
    real(dp) :: xmin,ymin,zmin,xmax,ymax,zmax
end type hex_bczone

!mesh vertex type 
type hex_vertex
    logical :: external,flag
    integer(in64) :: index,tag,idata 
    integer(in64), dimension(:), allocatable :: ivdata
    real(dp) :: rdata
    real(dp) :: gradient(3)
    real(dp) :: coordinate(3)
    real(dp), dimension(:,:), allocatable :: mesh_gradient
    contains
        procedure :: is_on_surface => is_on_surface_vertex
        procedure :: get_surface_edges => get_surface_edges_vertex
end type hex_vertex

!mesh edge type 
type hex_edge
    logical :: external,flag
    integer(in64) :: index,v1,v2,cell1,cell2,tag,idata 
    integer(in64), dimension(:), allocatable :: ivdata
    real(dp) :: rdata
    real(dp) :: direction(3)
    type(hex_vertex), pointer :: vertex1,vertex2
    contains
        procedure :: split => split_hex_edge
        procedure :: area => hex_edge_area
        procedure :: is_boundary => is_boundary_edge
end type hex_edge

!mesh cell type 
type hex_cell
    logical :: flag
    integer(in64) :: index,nedge,tag
    integer(in64), dimension(:), allocatable :: edges 
    real(dp) :: volume 
    real(dp) :: midpoint(3)
    contains 
        procedure :: get_midpoint => get_midpoint_cell
        procedure :: get_edge_loop => get_edge_loop_cell
end type hex_cell

!mesh type
type hex_mesh
    logical :: isvalid
    integer(in64) :: nedge,nvertex,ncell,max_valence,max_cvalence,nvertex_surfint
    integer(in64), dimension(:), allocatable :: valence
    integer(in64), dimension(:,:), allocatable :: v2e,v2c
    type(hex_cell), dimension(:), allocatable :: cell 
    type(hex_edge), dimension(:), allocatable :: edge
    type(hex_vertex), dimension(:), allocatable :: vertex
    contains 
        procedure :: extend => extend_hex_mesh
        procedure :: get_v2e
        procedure :: get_v2c
        procedure :: get_valence
        procedure :: get_cell_areas
        procedure :: index_vertices
        procedure :: index_edges
        procedure :: index_cells
        procedure :: get_cell_edges
        procedure :: get_cells
        procedure :: get_edges_shared_vertex
end type hex_mesh



!methods ==================================================
contains 

!get vertex shared between edges =========================
function get_edges_shared_vertex(self,edge1,edge2) result(vertex)
implicit none 

!mesh class
class(hex_mesh), target :: self

!variables - inout
integer(in64) :: edge1,edge2,vertex

!select vertex
if ((self%edge(edge1)%vertex1%index == self%edge(edge2)%vertex1%index) .OR. (self%edge(edge1)%vertex1%index == self%edge(edge2)%vertex2%index)) then 
    vertex = self%edge(edge1)%vertex1%index 
else
    vertex = self%edge(edge1)%vertex2%index 
end if 
return 
end function get_edges_shared_vertex


!edge is on boundary condition =========================
function is_boundary_edge(self,boundary_tag) result(is_on_btag)
implicit none 

!edge class
class(hex_edge), target :: self

!variables - inout
logical :: is_on_btag
integer(in64) :: boundary_tag

!check
if ((self%cell1 == boundary_tag) .OR. (self%cell2 == boundary_tag)) then 
    is_on_btag = .true.
else
    is_on_btag = .false.
end if 
return 
end function is_boundary_edge


!get the edges on this vertex that are surface edges =========================
function get_surface_edges_vertex(self,mesh) result(surface_edges)
implicit none 

!vertex class
class(hex_vertex), target :: self

!variables - inout
integer(in64), dimension(:), allocatable :: surface_edges
class(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii
integer(in64) :: nsurface

!search edges 
nsurface = 0
do ii=1,mesh%max_valence
    if (mesh%v2e(self%index,ii) .GT. 0) then 
        if ((mesh%edge(mesh%v2e(self%index,ii))%cell1 == -1) .OR. (mesh%edge(mesh%v2e(self%index,ii))%cell2 == -1)) then 
            nsurface = nsurface + 1
        end if 
    end if 
end do 

!accumulate edges 
allocate(surface_edges(nsurface))
surface_edges(:) = 0 
nsurface = 0
do ii=1,mesh%max_valence
    if (mesh%v2e(self%index,ii) .GT. 0) then 
        if ((mesh%edge(mesh%v2e(self%index,ii))%cell1 == -1) .OR. (mesh%edge(mesh%v2e(self%index,ii))%cell2 == -1)) then 
            nsurface = nsurface + 1
            surface_edges(nsurface) = mesh%v2e(self%index,ii)
        end if 
    end if 
end do 
return 
end function get_surface_edges_vertex 


!check if a vertex is on the geometry surface =========================
function is_on_surface_vertex(self,mesh) result(on_surface)
implicit none 

!vertex class
class(hex_vertex), target :: self

!variables - inout
logical :: on_surface
class(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii

!search edges 
on_surface = .false.
do ii=1,mesh%max_valence
    if (mesh%v2e(self%index,ii) .GT. 0) then 
        if ((mesh%edge(mesh%v2e(self%index,ii))%cell1 == -1) .OR. (mesh%edge(mesh%v2e(self%index,ii))%cell2 == -1)) then 
            on_surface = .true.
            exit
        end if 
    end if 
end do 
return 
end function is_on_surface_vertex 


!get cell edge loop =========================
function get_edge_loop_cell(self,mesh,start_with_non_bc_edge) result(loop)
implicit none 

!cell class
class(hex_cell), target :: self

!variables - inout
logical :: start_with_non_bc_edge
integer(in64) :: loop(self%nedge)
class(hex_mesh), target :: mesh

!variables - local 
logical :: nbcefound
integer(in64) :: ii,jj
integer(in64) :: lins,edgec,edgen,eadj
type(hex_vertex), pointer :: vertexc

!build loop 
lins = 2
loop(:) = 0 
if (start_with_non_bc_edge) then !find non boundary condition edge to start from 
    nbcefound = .false.
    do ii=1,self%nedge
        if ((mesh%edge(self%edges(ii))%cell1 .GT. 0) .AND. (mesh%edge(self%edges(ii))%cell2 .GT. 0)) then 
            loop(1) = self%edges(ii)
            nbcefound = .true.
            exit
        end if 
    end do 
    if (.NOT.nbcefound) then 
        loop(1) = self%edges(1)
    end if 
else !start from first edge
    loop(1) = self%edges(1)
end if 
edgec = loop(1) !self%edges(1)
do ii=1,self%nedge

    !find next edge on this cell 
    edgen = -1 
    if (mesh%edge(edgec)%cell1 == self%index) then 
        vertexc => mesh%edge(edgec)%vertex2
    else
        vertexc => mesh%edge(edgec)%vertex1
    end if 
    do jj=1,mesh%valence(vertexc%index)
        eadj = mesh%v2e(vertexc%index,jj)
        if (eadj .NE. edgec) then
            if ((mesh%edge(eadj)%cell1 == self%index) .OR. (mesh%edge(eadj)%cell2 == self%index)) then
                edgen = eadj
                exit
            end if 
        end if 
    end do 

    !exit if loop closed 
    if (edgen == loop(1)) then 
        exit 
    end if 
    edgec = edgen

    !add to loop 
    loop(lins) = edgec
    lins = lins + 1 
end do 
return 
end function get_edge_loop_cell


!get cell midpoint =========================
function get_midpoint_cell(self,mesh) result(midpoint)
implicit none 

!variables - inout
real(dp) :: midpoint(3)
class(hex_cell), target :: self
class(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii
real(dp) :: eweight,weight_total
real(dp) :: coordinate(3)

!evaluate midpoint
midpoint(:) = 0.0d0
weight_total = 0.0d0  
do ii=1,self%nedge
    coordinate = 0.5d0*(mesh%edge(self%edges(ii))%vertex1%coordinate + mesh%edge(self%edges(ii))%vertex2%coordinate)
    eweight = norm2(mesh%edge(self%edges(ii))%vertex1%coordinate - mesh%edge(self%edges(ii))%vertex2%coordinate)
    midpoint = midpoint + coordinate*eweight
    weight_total = weight_total + eweight
end do 
midpoint = midpoint/weight_total
return 
end function get_midpoint_cell


!get cell edges =========================
subroutine get_cell_edges(self)
implicit none 

!variables - import
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii

!assign edges to current cells
do ii=1,self%ncell
    self%cell(ii)%nedge = 0 
end do 
do ii=1,self%nedge
    if (self%edge(ii)%cell1 .GT. 0) then 
        self%cell(self%edge(ii)%cell1)%nedge = self%cell(self%edge(ii)%cell1)%nedge + 1
    end if 
    if (self%edge(ii)%cell2 .GT. 0) then 
        self%cell(self%edge(ii)%cell2)%nedge = self%cell(self%edge(ii)%cell2)%nedge + 1
    end if   
end do 
do ii=1,self%ncell
    if (allocated(self%cell(ii)%edges)) then 
        deallocate(self%cell(ii)%edges)
    end if 
    allocate(self%cell(ii)%edges(self%cell(ii)%nedge))
    self%cell(ii)%nedge = 0 
end do 
do ii=1,self%nedge
    if (self%edge(ii)%cell1 .GT. 0) then 
        self%cell(self%edge(ii)%cell1)%nedge = self%cell(self%edge(ii)%cell1)%nedge + 1
        self%cell(self%edge(ii)%cell1)%edges(self%cell(self%edge(ii)%cell1)%nedge) = self%edge(ii)%index
    end if 
    if (self%edge(ii)%cell2 .GT. 0) then 
        self%cell(self%edge(ii)%cell2)%nedge = self%cell(self%edge(ii)%cell2)%nedge + 1
        self%cell(self%edge(ii)%cell2)%edges(self%cell(self%edge(ii)%cell2)%nedge) = self%edge(ii)%index
    end if   
end do 
return
end subroutine get_cell_edges


!evaluate cell areas =========================
subroutine get_cell_areas(self)
implicit none

!variables - import
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii
real(dp) :: earea

!allocate cells 
if (.NOT. allocated(self%cell)) then 
    allocate(self%cell(self%ncell))
    do ii=1,self%ncell
        self%cell(ii)%flag = .false.
    end do 
elseif (size(self%cell,dim=1) .NE. self%ncell) then 
    deallocate(self%cell)
    allocate(self%cell(self%ncell))
    do ii=1,self%ncell
        self%cell(ii)%flag = .false.
    end do 
end if 

!evaluate volumes
do ii=1,self%ncell
    self%cell(ii)%volume = 0.0d0 
end do 
do ii=1,self%nedge
    earea = self%edge(ii)%area()
    if (self%edge(ii)%cell1 .GT. 0) then 
        self%cell(self%edge(ii)%cell1)%volume = self%cell(self%edge(ii)%cell1)%volume + earea
    end if
    if (self%edge(ii)%cell2 .GT. 0) then 
        self%cell(self%edge(ii)%cell2)%volume = self%cell(self%edge(ii)%cell2)%volume - earea
    end if  
end do 
return 
end subroutine get_cell_areas


!set vertex indecies ========================= 
subroutine index_vertices(self)
implicit none 

!variables - import
class(hex_mesh), target :: self
    
!variables - local 
integer(in64) :: ii

!index vertices
do ii=1,self%nvertex
    self%vertex(ii)%index = ii 
end do 
return 
end subroutine index_vertices


!set edge indecies =========================
subroutine index_edges(self)
implicit none 

!variables - import
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii

!index edges
do ii=1,self%nedge
    self%edge(ii)%index = ii 
end do 
return 
end subroutine index_edges


!set cell indecies =========================
subroutine index_cells(self)
implicit none 

!variables - inout 
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii  
integer(in64) :: ncellc,cins 
integer(in64), dimension(:), allocatable :: cell_map

!find current maximum cell count
ncellc = 0 
do ii=1,self%nedge
    if (self%edge(ii)%cell1 .GT. ncellc) then 
        ncellc = self%edge(ii)%cell1
    end if 
    if (self%edge(ii)%cell2 .GT. ncellc) then 
        ncellc = self%edge(ii)%cell2
    end if 
end do 

!map cells
cins = 0 
allocate(cell_map(ncellc))
cell_map(:) = 0 
do ii=1,self%nedge
    if (self%edge(ii)%cell1 .GT. 0) then 
        if (cell_map(self%edge(ii)%cell1) == 0) then 
            cins = cins + 1
            cell_map(self%edge(ii)%cell1) = cins
        end if 
    end if 
    if (self%edge(ii)%cell2 .GT. 0) then 
        if (cell_map(self%edge(ii)%cell2) == 0) then 
            cins = cins + 1
            cell_map(self%edge(ii)%cell2) = cins
        end if 
    end if 
end do 
do ii=1,self%nedge
    if (self%edge(ii)%cell1 .GT. 0) then 
        self%edge(ii)%cell1 = cell_map(self%edge(ii)%cell1)
    end if 
    if (self%edge(ii)%cell2 .GT. 0) then 
        self%edge(ii)%cell2 = cell_map(self%edge(ii)%cell2)
    end if 
end do 

!set the number of mesh cells 
self%ncell = cins
return 
end subroutine index_cells


!get cells =========================
subroutine get_cells(self)
implicit none 

!variables - import
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii  

!index cells 
call self%index_cells()

!allocate cells 
if (.NOT. allocated(self%cell)) then 
    allocate(self%cell(self%ncell))
    do ii=1,self%ncell
        self%cell(ii)%flag = .false.
    end do 
elseif (size(self%cell,dim=1) .NE. self%ncell) then 
    deallocate(self%cell)
    allocate(self%cell(self%ncell))
    do ii=1,self%ncell
        self%cell(ii)%flag = .false.
    end do 
end if 

!set cell indecies
do ii=1,self%nedge
    if (self%edge(ii)%cell1 .GT. 0) then 
        self%cell(self%edge(ii)%cell1)%index = self%edge(ii)%cell1
    end if 
    if (self%edge(ii)%cell2 .GT. 0) then 
        self%cell(self%edge(ii)%cell2)%index = self%edge(ii)%cell2
    end if 
end do 

!get cell edges 
call self%get_cell_edges()
return 
end subroutine get_cells


!get v2c =========================
subroutine get_v2c(self)
implicit none 

!variables - import
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii,jj,ee
integer(in64) :: vtgt
integer(in64) :: valence(self%nvertex)

!index cell 
call self%index_cells()

!assign edges to current cells 
call self%get_cell_edges()

!set vertex flags
do ii=1,self%nvertex
    self%vertex(ii)%flag = .false.
end do 

!evaluate maximum valence
self%max_cvalence = 0 
valence(:) = 0 
do ii=1,self%ncell
    do ee=1,self%cell(ii)%nedge
        if (.NOT. self%edge(self%cell(ii)%edges(ee))%vertex1%flag) then 
            vtgt = self%edge(self%cell(ii)%edges(ee))%vertex1%index
            valence(vtgt) = valence(vtgt) + 1
            self%edge(self%cell(ii)%edges(ee))%vertex1%flag = .true.
        end if 
        if (.NOT. self%edge(self%cell(ii)%edges(ee))%vertex2%flag) then 
            vtgt = self%edge(self%cell(ii)%edges(ee))%vertex2%index
            valence(vtgt) = valence(vtgt) + 1
            self%edge(self%cell(ii)%edges(ee))%vertex2%flag = .true.
        end if 
    end do 
    do ee=1,self%cell(ii)%nedge
        self%edge(self%cell(ii)%edges(ee))%vertex1%flag = .false.
        self%edge(self%cell(ii)%edges(ee))%vertex2%flag = .false.
    end do 
end do 
self%max_cvalence = maxval(valence)

!build v2c
if (allocated(self%v2c)) then 
    deallocate(self%v2c)
end if 
allocate(self%v2c(self%nvertex,self%max_cvalence))
self%v2c(:,:) = 0 
do ii=1,self%ncell
    do ee=1,self%cell(ii)%nedge
        if (.NOT. self%edge(self%cell(ii)%edges(ee))%vertex1%flag) then 
            vtgt = self%edge(self%cell(ii)%edges(ee))%vertex1%index
            do jj=1,self%max_cvalence
                if (self%v2c(vtgt,jj) == 0) then 
                    self%v2c(vtgt,jj) = ii
                    exit 
                end if 
            end do 
            self%edge(self%cell(ii)%edges(ee))%vertex1%flag = .true.
        end if 
        if (.NOT. self%edge(self%cell(ii)%edges(ee))%vertex2%flag) then 
            vtgt = self%edge(self%cell(ii)%edges(ee))%vertex2%index
            do jj=1,self%max_cvalence
                if (self%v2c(vtgt,jj) == 0) then 
                    self%v2c(vtgt,jj) = ii
                    exit 
                end if 
            end do 
            self%edge(self%cell(ii)%edges(ee))%vertex2%flag = .true.
        end if 
    end do 
    do ee=1,self%cell(ii)%nedge
        self%edge(self%cell(ii)%edges(ee))%vertex1%flag = .false.
        self%edge(self%cell(ii)%edges(ee))%vertex2%flag = .false.
    end do 
end do 
return 
end subroutine get_v2c


!get v2e =========================
subroutine get_v2e(self)
implicit none 

!variables - import
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii,ee

!evaluate maximum valence
call self%get_valence()

!build v2e
if (allocated(self%v2e)) then 
    deallocate(self%v2e)
end if 
allocate(self%v2e(self%nvertex,self%max_valence))
self%v2e(:,:) = 0 
do ii=1,self%nedge
    do ee=1,self%max_valence
        if (self%v2e(self%edge(ii)%vertex1%index,ee) == 0) then 
            self%v2e(self%edge(ii)%vertex1%index,ee) = self%edge(ii)%index
            exit 
        end if 
    end do 
    do ee=1,self%max_valence
        if (self%v2e(self%edge(ii)%vertex2%index,ee) == 0) then 
            self%v2e(self%edge(ii)%vertex2%index,ee) = self%edge(ii)%index
            exit
        end if 
    end do 
end do 
return 
end subroutine get_v2e


!get valence =========================
subroutine get_valence(self)
implicit none 

!variables - import
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii

!evaluate valence
self%max_valence = 0 
if (allocated(self%valence)) then 
    deallocate(self%valence)
end if 
allocate(self%valence(self%nvertex))
self%valence(:) = 0
do ii=1,self%nedge
    self%valence(self%edge(ii)%vertex1%index) = self%valence(self%edge(ii)%vertex1%index) + 1
    self%valence(self%edge(ii)%vertex2%index) = self%valence(self%edge(ii)%vertex2%index) + 1
end do 
self%max_valence = maxval(self%valence)
return 
end subroutine get_valence


!edge area contribution =========================
function hex_edge_area(self) result(area)
implicit none 

!variables - inout
real(dp) :: area
class(hex_edge), target :: self

!variables - local
real(dp) :: v1(2),v2(2)

!evaluate
v1 = self%vertex1%coordinate(1:2)
v2 = self%vertex2%coordinate(1:2)
area = 0.5d0*(v1(1)*v2(2) - v2(1)*v1(2))
return 
end function hex_edge_area


!split edge =========================
subroutine split_hex_edge(self,mesh,vn,enew,vnew)
implicit none 

!variables - import
integer(in64) :: enew,vnew
class(hex_edge), target :: self 
class(hex_mesh), target :: mesh
real(dp), dimension(:) :: vn 

!add new vertex
mesh%nvertex = mesh%nvertex + 1
mesh%vertex(mesh%nvertex)%tag = 0 
mesh%vertex(mesh%nvertex)%coordinate = vn
mesh%vertex(mesh%nvertex)%index = mesh%nvertex
vnew = mesh%nvertex

!add new edge 
mesh%nedge = mesh%nedge + 1
mesh%edge(mesh%nedge)%tag = 0 
mesh%edge(mesh%nedge)%cell1 = self%cell1
mesh%edge(mesh%nedge)%cell2 = self%cell2
mesh%edge(mesh%nedge)%vertex1 => mesh%vertex(mesh%nvertex)
mesh%edge(mesh%nedge)%vertex2 => self%vertex2 
mesh%edge(mesh%nedge)%index = mesh%nedge
mesh%edge(mesh%nedge)%direction = self%direction
enew = mesh%nedge

!update current edge endpoint 
self%vertex2 => mesh%vertex(mesh%nvertex)
return 
end subroutine split_hex_edge


!extend mesh =========================
subroutine extend_hex_mesh(self,size_set)  
implicit none 

!variables - import
integer(in64) :: size_set
class(hex_mesh), target :: self

!variables - local 
integer(in64) :: ii,len_vertex,len_edge
integer(in64) :: edge_vindex(self%nedge,2)
type(hex_vertex), dimension(:), allocatable :: vertex_temp 
type(hex_edge), dimension(:), allocatable :: edge_temp 

!index edge-vertex links
do ii=1,self%nedge
    edge_vindex(ii,1) = self%edge(ii)%vertex1%index
    edge_vindex(ii,2) = self%edge(ii)%vertex2%index
end do 

!extend self
if (size_set == 0) then 
    len_vertex = 2*size(self%vertex,dim=1)
    len_edge = 2*size(self%edge,dim=1)
else
    len_vertex = size_set
    len_edge = size_set
end if 
allocate(vertex_temp(self%nvertex))
allocate(edge_temp(self%nedge))
vertex_temp = self%vertex
edge_temp = self%edge
deallocate(self%vertex)
allocate(self%vertex(len_vertex))
self%vertex(1:self%nvertex) = vertex_temp
deallocate(self%edge)
allocate(self%edge(len_edge))
self%edge(1:self%nedge) = edge_temp

!reassign edge-vertex pointers
do ii=1,len_edge
    self%edge(ii)%vertex1 => null()
    self%edge(ii)%vertex2 => null()
end do 
do ii=1,self%nedge
    self%edge(ii)%vertex1 => self%vertex(edge_vindex(ii,1))
    self%edge(ii)%vertex2 => self%vertex(edge_vindex(ii,2))
end do 
return 
end subroutine extend_hex_mesh


end module hex_data_methods