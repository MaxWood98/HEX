!halfedge mesh structure main class
!max wood
!version : 0.0.9
!updated : 17-02-25

!module
module halfedge_mesh

!dependancies 
use facevertex_mesh

!Integer data types 
use ISO_FORTRAN_ENV, only: in16=>int16
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in64=>int64

!Real data types 
use ISO_FORTRAN_ENV, only: sp=>real32 
use ISO_FORTRAN_ENV, only: dp=>real64

!items 
implicit none 

!constants 
integer(in64), parameter :: face_max_edges_he = 1000000

!vertex type 
type vertex
    integer(in64) :: index,tag
    real(dp) :: rdata
    real(dp), dimension(:), allocatable :: coordinate  
    type(edge), pointer :: edge 
    contains 
        procedure :: get_valence => get_valence_v
        procedure :: get_connected_edges => get_connected_edges_v
end type vertex

!edge type 
type edge
    logical :: flag 
    integer(in64) :: index,tag
    real(dp) :: rdata
    real(dp), dimension(:), allocatable :: direction,normal
    type(vertex), pointer :: origin,origin0
    type(edge), pointer :: opposite,previous,next
    type(face), pointer :: adjacent 
    contains 
        procedure :: split => split_e
end type edge

!face type 
type face
    integer(in64) :: index,tag
    real(dp) :: rdata
    type(edge), pointer :: edge   
    contains 
        procedure :: get_normal => get_normal_f
        procedure :: get_area => get_area_f
        procedure :: get_edges => get_edges_f
        procedure :: get_number_of_edges => get_number_of_edges_f
end type face

!halfedge type 
type halfedge 
    integer(in64) :: nvertex,nedge,nface,ndim
    type(vertex), dimension(:), allocatable :: vertex
    type(edge), dimension(:), allocatable :: edge 
    type(face), dimension(:), allocatable :: face 
    contains 
        procedure :: build_from_fv
        procedure :: set_indecies
        procedure :: extend => extend_halfedge
end type halfedge 


!methods ==================================================
contains 


!general methods ==============


!build from face vertex .fv file ========================= ***** check coordinate dimensions 2D/3D??????????????
subroutine build_from_fv(self,filename) 
implicit none

!variables - import
class(halfedge), target :: self
character(*), intent(in) :: filename

!variables - local
integer(in64) :: ii,ee,nhalfedge
integer(in64) :: ev1,ev2,eidx,efidx,eidx_p,he1,he2
integer(in64), dimension(:,:), allocatable :: edge_idx,edge_idx_f
type(facevertex) :: mesh

!read .fv mesh 
call mesh%read_fv_file(filename)
if (mesh%nedge == 0) then !build edges if not supplied
    call mesh%build_edges_from_faces()
end if 
if (mesh%ndim == 2) then !initialise pseudo-faces if mesh is 2d
    call mesh%build_edge_loop_faces()
end if 

!allocate halfedge structure 
allocate(self%vertex(mesh%nvertex))
self%nvertex = mesh%nvertex
allocate(self%edge(2*mesh%nedge))
self%nedge = 2*mesh%nedge
allocate(self%face(mesh%nface))
self%nface = mesh%nface
self%ndim = mesh%ndim

!initialise vertices 
do ii=1,self%nvertex
    allocate(self%vertex(ii)%coordinate(3)) !force to 3D
    self%vertex(ii)%coordinate = mesh%vertices(ii,:)
end do 

!initialise edges 
do ii=1,self%nedge
    allocate(self%edge(ii)%direction(3)) !force to 3D
    self%edge(ii)%direction = 0.0d0 
    allocate(self%edge(ii)%normal(3)) !force to 3D
    self%edge(ii)%normal = 0.0d0
end do 

!set indecies 
call self%set_indecies()

!build connectivity
call mesh%get_connectivity()

!construct half edges 
nhalfedge = 0 
allocate(edge_idx(mesh%nedge,2))
allocate(edge_idx_f(mesh%nedge,2))
edge_idx(:,:) = 0
edge_idx_f(:,:) = 0 
do ii=1,mesh%nface !build edges

    !add half edges if required in this face
    do ee=1,mesh%faces(ii)%nvertex
        ev1 = mesh%faces(ii)%vertices(ee)
        ev2 = mesh%faces(ii)%vertices(mod(ee,mesh%faces(ii)%nvertex) + 1)
        eidx = mesh%faces(ii)%edges(ee)
        nhalfedge = nhalfedge + 1
        self%edge(nhalfedge)%origin => self%vertex(ev1)
        self%edge(nhalfedge)%origin0 => self%vertex(ev1)
        self%edge(nhalfedge)%adjacent => self%face(ii)
        if (edge_idx(eidx,1) == 0) then 
            edge_idx(eidx,1) = nhalfedge
            edge_idx_f(eidx,1) = ii
        else
            edge_idx(eidx,2) = nhalfedge
            edge_idx_f(eidx,2) = ii
        end if 
        if (ee == 1) then !store halfedge for this face
            efidx = nhalfedge
        end if 
    end do 

    !assign edge to this face
    self%face(ii)%edge => self%edge(efidx)

    !link next and previous in this face
    do ee=1,mesh%faces(ii)%nvertex
        ev1 = ee
        ev2 = mod(ee,mesh%faces(ii)%nvertex) + 1
        if (edge_idx_f(mesh%faces(ii)%edges(ev1),1) == ii) then 
            he1 = edge_idx(mesh%faces(ii)%edges(ev1),1)
        else
            he1 = edge_idx(mesh%faces(ii)%edges(ev1),2)
        end if
        if (edge_idx_f(mesh%faces(ii)%edges(ev2),1) == ii) then 
            he2 = edge_idx(mesh%faces(ii)%edges(ev2),1)
        else
            he2 = edge_idx(mesh%faces(ii)%edges(ev2),2)
        end if
        self%edge(he1)%next => self%edge(he2)
        self%edge(he2)%previous => self%edge(he1)
    end do 
end do 
do ii=1,mesh%nedge !build halfedges on any shell edges 
    if (edge_idx(ii,2) == 0) then 

        !current half edge
        he1 = edge_idx(ii,1)

        !select the correct origin for the new halfedge
        if (self%edge(he1)%origin%index == mesh%edges(ii,1)) then 
            ev1 = mesh%edges(ii,2)
        else
            ev1 = mesh%edges(ii,1)
        end if 

        !build the new halfedge 
        nhalfedge = nhalfedge + 1
        self%edge(nhalfedge)%origin => self%vertex(ev1)
        self%edge(nhalfedge)%origin0 => self%vertex(ev1)
        self%edge(nhalfedge)%adjacent => null()
        
        !store in index array 
        edge_idx(ii,2) = -nhalfedge
    end if 
end do 
do ii=1,mesh%nedge !link previous-next for any halfedges on shell edges 
    if (edge_idx(ii,2) < 0) then 

        !origin index
        ev1 = self%edge(abs(edge_idx(ii,2)))%origin%index

        !find the other connected shell edge
        eidx = 0 
        do ee=1,mesh%valence(ev1)
            if (mesh%v2e(ev1,ee) .NE. ii) then
                if ((mesh%e2f(mesh%v2e(ev1,ee),1) == 0) .OR. (mesh%e2f(mesh%v2e(ev1,ee),2) == 0)) then
                    eidx = mesh%v2e(ev1,ee)
                    exit
                end if 
            end if 
        end do 
        if (eidx == 0) then 
            print *, '** unlinked shell edge : ', abs(edge_idx(ii,2))
        end if 

        !find the shell edge half edge on eidx
        if (edge_idx(eidx,1) < 0) then 
            eidx_p = abs(edge_idx(eidx,1))
        else
            eidx_p = abs(edge_idx(eidx,2))
        end if 

        !link the edges
        self%edge(abs(edge_idx(ii,2)))%previous => self%edge(eidx_p)
        self%edge(eidx_p)%next => self%edge(abs(edge_idx(ii,2)))
    end if 
end do 
do ii=1,mesh%nedge !link opposites
    he1 = abs(edge_idx(ii,1))
    he2 = abs(edge_idx(ii,2))
    self%edge(he1)%opposite => self%edge(he2)
    self%edge(he2)%opposite => self%edge(he1)
end do 

!assign edge to each vertex
do ii=1,self%nedge
    ev1 = self%edge(ii)%origin%index 
    self%vertex(ev1)%edge => self%edge(ii)
end do 

!set the direction of each edge 
do ii=1,self%nedge
    self%edge(ii)%direction = self%edge(ii)%opposite%origin%coordinate - self%edge(ii)%origin%coordinate
    self%edge(ii)%direction = self%edge(ii)%direction/norm2(self%edge(ii)%direction)
end do 
return 
end subroutine build_from_fv


!set indecies =========================
subroutine set_indecies(self) 
implicit none 

!variables - import
class(halfedge) :: self 

!variables - local
integer(in64) :: ii 

!set vertices
do ii=1,self%nvertex
    self%vertex(ii)%index = ii 
end do 

!set edges
do ii=1,self%nedge
    self%edge(ii)%index = ii 
end do 

!set faces
do ii=1,self%nface
    self%face(ii)%index = ii 
end do 
return 
end subroutine set_indecies


!extend (double the length of the vertex - edge - face arrays) =========================
subroutine extend_halfedge(self,size_set)
implicit none 

!variables - import
integer(in64) :: size_set
class(halfedge), target :: self

!variables - local 
integer(in64) :: ii,lencrd,lenitem
integer(in64), dimension(:), allocatable :: vedge,fedge
integer(in64), dimension(:,:), allocatable :: edge_data
real(dp), dimension(:,:), allocatable :: vcoordinates,edirs,enrms

!set indexing  
call self%set_indecies()

!get coordinate dimension 
lencrd = size(self%vertex(1)%coordinate,dim=1)

!store edge data
allocate(edge_data(self%nedge,6)) !origin | origin0 | opposite | previous | next | adjacent 
allocate(edirs(self%nedge,lencrd))
allocate(enrms(self%nedge,lencrd))
do ii=1,self%nedge
    edge_data(ii,1) = self%edge(ii)%origin%index
    edge_data(ii,2) = self%edge(ii)%origin0%index
    edge_data(ii,3) = self%edge(ii)%opposite%index
    edge_data(ii,4) = self%edge(ii)%previous%index
    edge_data(ii,5) = self%edge(ii)%next%index
    if (associated(self%edge(ii)%adjacent)) then 
        edge_data(ii,6) = self%edge(ii)%adjacent%index
    else
        edge_data(ii,6) = -1 !tag as unassoiated
    end if 
    edirs(ii,:) = self%edge(ii)%direction
    enrms(ii,:) = self%edge(ii)%normal
end do 

!store vertex data
allocate(vedge(self%nvertex))
allocate(vcoordinates(self%nvertex,lencrd))
do ii=1,self%nvertex
    vedge(ii) = self%vertex(ii)%edge%index
    vcoordinates(ii,:) = self%vertex(ii)%coordinate
end do 

!store face data
allocate(fedge(self%nface))
do ii=1,self%nface
    fedge(ii) = self%face(ii)%edge%index
end do 

!extend edge
if (size_set == 0) then  
    lenitem = 2*size(self%edge,dim=1)
else
    lenitem = size_set
end if 
deallocate(self%edge)
allocate(self%edge(lenitem))
do ii=1,lenitem
    self%edge(ii)%origin => null()
    self%edge(ii)%origin0 => null()
    self%edge(ii)%opposite => null()
    self%edge(ii)%previous => null()
    self%edge(ii)%next => null()
    self%edge(ii)%adjacent => null()
    allocate(self%edge(ii)%direction(lencrd))
    self%edge(ii)%direction = 0.0d0 
    allocate(self%edge(ii)%normal(lencrd))
    self%edge(ii)%normal = 0.0d0 
end do 

!extend vertex
if (size_set == 0) then  
    lenitem = 2*size(self%vertex,dim=1)
else
    lenitem = size_set
end if 
deallocate(self%vertex)
allocate(self%vertex(lenitem))
do ii=1,lenitem
    self%vertex(ii)%edge => null()
    allocate(self%vertex(ii)%coordinate(lencrd))
    self%vertex(ii)%coordinate(:) = 0.0d0 
end do 

!extend face 
if (size_set == 0) then  
    lenitem = 2*size(self%face,dim=1)
else
    lenitem = size_set
end if 
deallocate(self%face)
allocate(self%face(lenitem))
do ii=1,lenitem
    self%face(ii)%edge => null()
end do 

!reassign edge data 
do ii=1,self%nedge
    self%edge(ii)%origin => self%vertex(edge_data(ii,1))
    self%edge(ii)%origin0 => self%vertex(edge_data(ii,2))
    self%edge(ii)%opposite => self%edge(edge_data(ii,3))
    self%edge(ii)%previous => self%edge(edge_data(ii,4))
    self%edge(ii)%next => self%edge(edge_data(ii,5))
    if (edge_data(ii,6) == -1) then 
        self%edge(ii)%adjacent => null()
    else
        self%edge(ii)%adjacent => self%face(edge_data(ii,6))
    end if 
    self%edge(ii)%direction = edirs(ii,:)
    self%edge(ii)%normal = enrms(ii,:)
end do 

!reassign vertex data 
do ii=1,self%nvertex
    self%vertex(ii)%edge => self%edge(vedge(ii))
    self%vertex(ii)%coordinate = vcoordinates(ii,:)
end do 

!reassign face data 
do ii=1,self%nface 
    self%face(ii)%edge => self%edge(fedge(ii))
end do 

!reconstruct indexing 
call self%set_indecies()
return 
end subroutine extend_halfedge


!general methods ==============
!=======================================================
!vertex methods ===============


!get connected edges =========================
subroutine get_connected_edges_v(self,edges)
implicit none 

!variables - import
class(vertex), intent(in) :: self
type(edge), dimension(:), allocatable :: edges

!variables - local
integer(in64) :: ii 
integer(in64) :: valence 
type(edge), pointer :: edgec

!count edges 
valence = self%get_valence()

!accumulate edges
if (allocated(edges)) then 
    deallocate(edges)
end if 
allocate(edges(valence))
edgec => self%edge
do ii=1,valence
    edges(ii) = edgec
    edgec => edgec%opposite%next
end do 
return 
end subroutine get_connected_edges_v


!get valence =========================
function get_valence_v(self) result(valence)
implicit none 

!variables - import
class(vertex), intent(in) :: self
integer(in64) :: valence

!variables - local
integer(in64) :: ii 
type(edge), pointer :: edgec

!count edges 
ii = 0 
valence = 0
edgec => self%edge
do while (ii == 0)
    valence = valence + 1 
    edgec => edgec%opposite%next
    if (associated(edgec,self%edge)) then
        ii = 1
        exit 
    end if 
end do 
return 
end function get_valence_v


!vertex methods ===============
!=======================================================
!edge methods =================


!split edge =========================
subroutine split_e(self,mesh,vn,nvidx,neidx) 
implicit none  

!variables - import
class(edge), target :: self 
class(halfedge), target :: mesh
integer(in64) :: nvidx
integer(in64) :: neidx(2)
real(dp), dimension(:) :: vn

!index of the new edges 
neidx(1) = mesh%nedge + 1 !main 
neidx(2) = neidx(1) + 1 !opposite

!index of the new vertex
nvidx = mesh%nvertex + 1

!create the new vertex
mesh%vertex(nvidx)%tag = 0 
mesh%vertex(nvidx)%index = nvidx
mesh%vertex(nvidx)%coordinate = vn

!create the new edge
mesh%edge(neidx(1))%index = neidx(1)
mesh%edge(neidx(1))%previous => self
mesh%edge(neidx(1))%next => self%next 
mesh%edge(neidx(1))%adjacent => self%adjacent 
mesh%edge(neidx(1))%opposite => mesh%edge(neidx(2))
mesh%edge(neidx(1))%origin => mesh%vertex(nvidx)
mesh%edge(neidx(1))%origin0 => self%origin0
mesh%edge(neidx(1))%direction = self%direction
mesh%edge(neidx(1))%normal = self%normal
mesh%edge(neidx(1))%tag = self%tag

!create the new opposite edge
mesh%edge(neidx(2))%index = neidx(2)
mesh%edge(neidx(2))%previous => self%opposite%previous 
mesh%edge(neidx(2))%next => self%opposite
mesh%edge(neidx(2))%adjacent => self%opposite%adjacent 
mesh%edge(neidx(2))%opposite => mesh%edge(neidx(1))
mesh%edge(neidx(2))%origin => self%opposite%origin 
mesh%edge(neidx(2))%origin0 => self%opposite%origin0 
mesh%edge(neidx(2))%direction = self%opposite%direction
mesh%edge(neidx(2))%normal = self%opposite%normal
mesh%edge(neidx(2))%tag = self%opposite%tag

!update the surrounding connectivity
self%next%previous => mesh%edge(neidx(1))
self%opposite%previous%next => mesh%edge(neidx(2))

!update the connecetivity of the original edge and its opposite
self%next => mesh%edge(neidx(1))
self%opposite%previous => mesh%edge(neidx(2))

!update the origin of the original opposite edge 
self%opposite%origin => mesh%vertex(nvidx)

!add edge to the new vertex
mesh%vertex(nvidx)%edge => mesh%edge(neidx(1))

!increment the edge and vertex counts 
mesh%nvertex = mesh%nvertex + 1
mesh%nedge = mesh%nedge + 2
return 
end subroutine split_e


!edge methods =================
!=======================================================
!face methods =================


!get area =========================
function get_area_f(self) result(area)
implicit none 

!variables - import
class(face), intent(in) :: self
real(dp) :: area

!evaluate area
area = norm2(self%get_normal(unit_length = .False.))
return 
end function get_area_f


!get normal =========================
function get_normal_f(self,unit_length) result(normal)
implicit none 

!variables - import
class(face), intent(in) :: self
real(dp) :: normal(size(self%edge%origin%coordinate,1))
logical, optional :: unit_length

!variables - local 
integer(in64) :: ii
real(dp) :: vtxC(size(self%edge%origin%coordinate,1)),vtxN(size(self%edge%origin%coordinate,1))
type(edge), allocatable, dimension(:) :: edges

!set default unit_length
if (.NOT. present(unit_length)) then
    unit_length = .FALSE.
end if 

!get edges of this face
call self%get_edges(edges)

!accumulate normal
normal(:) = 0.0d0 
do ii=1,size(edges,1)
    vtxC(:) = edges(ii)%origin%coordinate
    vtxN(:) = edges(ii)%opposite%origin%coordinate
    normal(1) = normal(1) - 0.5d0*(vtxN(3) + vtxC(3))*(vtxN(2) - vtxC(2))
    normal(2) = normal(2) - 0.5d0*(vtxN(1) + vtxC(1))*(vtxN(3) - vtxC(3))
    normal(3) = normal(3) - 0.5d0*(vtxN(2) + vtxC(2))*(vtxN(1) - vtxC(1))
end do 
if (unit_length) then 
    normal(:) = normal(:)/norm2(normal)
end if 
return 
end function get_normal_f


!get edges =========================
subroutine get_edges_f(self,edges) 
implicit none 

!variables - import
class(face), intent(in) :: self
type(edge), allocatable, dimension(:) :: edges

!variables - local
integer(in64) :: ii 
integer(in64) :: nedge
type(edge), pointer :: edgec 

!count and accumulate edges 
nedge = self%get_number_of_edges() 
if (allocated(edges)) then
    deallocate(edges)
end if 
allocate(edges(nedge))
edgec => self%edge
do ii=1,face_max_edges_he
    edges(ii) = edgec
    edgec => edgec%next
    if (associated(edgec,self%edge)) then
        exit  
    end if 
end do 
return 
end subroutine get_edges_f


!get number of edges =========================
function get_number_of_edges_f(self) result(nedge)
implicit none 

!variables - import 
integer(in64) :: nedge
class(face), intent(in) :: self

!variables - local
integer(in64) :: ii 
type(edge), pointer :: edgec 

!count edges
nedge = 0 
edgec => self%edge
do ii=1,face_max_edges_he
    nedge = nedge + 1
    edgec => edgec%next
    if (associated(edgec,self%edge)) then
        exit  
    end if 
end do 
return 
end function get_number_of_edges_f


!face methods =================


end module halfedge_mesh 
