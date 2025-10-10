!hex data module 
!max wood
!version : 0.0.2
!updated : 19-01-25

!module 
module hex_data

!integer data types 
use ISO_FORTRAN_ENV, only: in16=>int16
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in64=>int64

!real data types 
use ISO_FORTRAN_ENV, only: sp=>real32 
use ISO_FORTRAN_ENV, only: dp=>real64

!options type
type hex_options
    logical :: cdisplay
    character(len=:), allocatable :: geomname 
    integer(in64) :: nrefine,nflood
    integer(in64) :: ncell_max
    real(dp) :: farfield_r
end type hex_options 

!mesh vertex type 
type hex_vertex
    logical :: external
    integer(in64) :: index,tag 
    real(dp) :: coordinate(3)
end type hex_vertex

!mesh edge type 
type hex_edge
    logical :: external
    integer(in64) :: index,cell1,cell2,tag 
    type(hex_vertex), pointer :: vertex1,vertex2
    contains
        procedure :: split => split_hex_edge
end type hex_edge

!mesh type
type hex_mesh
    logical :: isvalid
    integer(in64) :: nedge,nvertex
    type(hex_edge), dimension(:), allocatable :: edge
    type(hex_vertex), dimension(:), allocatable :: vertex
    contains 
        procedure :: extend => extend_hex_mesh
end type hex_mesh

!tritree type
type tritree
    integer(in64) :: nnode,nvertex
    integer(in64), dimension(:), allocatable :: vertex_tag
    real(dp), dimension(:,:), allocatable :: vertices 
    type(trinode), dimension(:), allocatable :: node
    contains 
        procedure :: refine_node
end type tritree

!tri node type
type trinode
    logical :: flag = .false.,visited = .false.
    integer(in64) :: level = 1 
    integer(in64) :: index = 0
    integer(in64) :: vertex1 = 0,vertex2 = 0,vertex3 = 0
    integer(in64) :: evertex1 = 0,evertex2 = 0,evertex3 = 0
    integer(in64), pointer :: aevertex1 => null(),aevertex2 => null(),aevertex3 => null()
    type(trinode), pointer :: parent => null()
    type(trinode), pointer :: child1 => null()
    type(trinode), pointer :: child2 => null()
    type(trinode), pointer :: child3 => null()
    type(trinode), pointer :: child4 => null()
    type(trinode), pointer :: adjacent1 => null()
    type(trinode), pointer :: adjacent2 => null()
    type(trinode), pointer :: adjacent3 => null()
    contains 
        procedure :: map_adjacent_node
        procedure :: cascade_adjacency
end type trinode


!methods ==================================================
contains 


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
enew = mesh%nedge

!update current edge endpoint 
self%vertex2 => mesh%vertex(mesh%nvertex)
return 
end subroutine split_hex_edge


!extend mesh =========================
subroutine extend_hex_mesh(self)  
implicit none 

!variables - import
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
len_vertex = size(self%vertex,dim=1)
len_edge = size(self%edge,dim=1)
allocate(vertex_temp(self%nvertex))
allocate(edge_temp(self%nedge))
vertex_temp = self%vertex
edge_temp = self%edge
deallocate(self%vertex)
allocate(self%vertex(2*len_vertex))
self%vertex(1:self%nvertex) = vertex_temp
deallocate(self%edge)
allocate(self%edge(2*len_edge))
self%edge(1:self%nedge) = edge_temp

!reassign edge-vertex pointers
do ii=1,2*len_edge
    self%edge(ii)%vertex1 => null()
    self%edge(ii)%vertex2 => null()
end do 
do ii=1,self%nedge
    self%edge(ii)%vertex1 => self%vertex(edge_vindex(ii,1))
    self%edge(ii)%vertex2 => self%vertex(edge_vindex(ii,2))
end do 
return 
end subroutine extend_hex_mesh


!refine trinode =========================
subroutine refine_node(tri_tree,node)
implicit none 

!variables - import 
class(tritree), target :: tri_tree
class(trinode), target :: node

!variables - local 
type(trinode), pointer :: child1,child2,child3,child4 
type(trinode) :: node_map1,node_map2,node_map3

!build adjacent node maps
if (associated(node%adjacent1)) then 
    call node%map_adjacent_node(node%adjacent1,node_map1)
end if 
if (associated(node%adjacent2)) then 
    call node%map_adjacent_node(node%adjacent2,node_map2)
end if 
if (associated(node%adjacent3)) then 
    call node%map_adjacent_node(node%adjacent3,node_map3)
end if 

!build new edge vertices 
if (node%evertex1 == 0) then 
    tri_tree%nvertex = tri_tree%nvertex + 1
    tri_tree%vertices(tri_tree%nvertex,:) = 0.5d0*(tri_tree%vertices(node%vertex1,:) + tri_tree%vertices(node%vertex2,:))
    node%evertex1 = tri_tree%nvertex
    if (associated(node_map1%aevertex3)) then 
        node_map1%aevertex3 = tri_tree%nvertex
    end if 
end if 
if (node%evertex2 == 0) then 
    tri_tree%nvertex = tri_tree%nvertex + 1
    tri_tree%vertices(tri_tree%nvertex,:) = 0.5d0*(tri_tree%vertices(node%vertex2,:) + tri_tree%vertices(node%vertex3,:))
    node%evertex2 = tri_tree%nvertex
    if (associated(node_map2%aevertex3)) then 
        node_map2%aevertex3 = tri_tree%nvertex
    end if 
end if 
if (node%evertex3 == 0) then 
    tri_tree%nvertex = tri_tree%nvertex + 1
    tri_tree%vertices(tri_tree%nvertex,:) = 0.5d0*(tri_tree%vertices(node%vertex3,:) + tri_tree%vertices(node%vertex1,:))
    node%evertex3 = tri_tree%nvertex
    if (associated(node_map3%aevertex3)) then 
        node_map3%aevertex3 = tri_tree%nvertex
    end if 
end if 

!build child nodes
tri_tree%nnode = tri_tree%nnode + 1
child1 => tri_tree%node(tri_tree%nnode)
tri_tree%nnode = tri_tree%nnode + 1
child2 => tri_tree%node(tri_tree%nnode)
tri_tree%nnode = tri_tree%nnode + 1
child3 => tri_tree%node(tri_tree%nnode)
tri_tree%nnode = tri_tree%nnode + 1
child4 => tri_tree%node(tri_tree%nnode)

!set child 1
node%child1 => child1
child1%level = node%level + 1
child1%parent => node  
child1%vertex1 = node%vertex1
child1%vertex2 = node%evertex1
child1%vertex3 = node%evertex3
child1%adjacent2 => child4 
if (associated(node_map1%child1)) then 
    child1%adjacent1 => node_map1%child1
    call node_map1%child1%cascade_adjacency(child1%parent,child1)
else
    child1%adjacent1 => node%adjacent1
end if 
if (associated(node_map3%child3)) then 
    child1%adjacent3 => node_map3%child3
    call node_map3%child3%cascade_adjacency(child1%parent,child1)
else
    child1%adjacent3 => node%adjacent3
end if 

!set child 2
node%child2 => child2
child2%level = node%level + 1
child2%parent => node  
child2%vertex1 = node%evertex1
child2%vertex2 = node%vertex2
child2%vertex3 = node%evertex2
child2%adjacent3 => child4
if (associated(node_map1%child3)) then 
    child2%adjacent1 => node_map1%child3
    call node_map1%child3%cascade_adjacency(child2%parent,child2)
else
    child2%adjacent1 => node%adjacent1
end if 
if (associated(node_map2%child1)) then 
    child2%adjacent2 => node_map2%child1
    call node_map2%child1%cascade_adjacency(child2%parent,child2)
else
    child2%adjacent2 => node%adjacent2
end if 

!set child 3
node%child3 => child3
child3%level = node%level + 1
child3%parent => node  
child3%vertex1 = node%evertex3
child3%vertex2 = node%evertex2
child3%vertex3 = node%vertex3
child3%adjacent1 => child4
if (associated(node_map2%child3)) then 
    child3%adjacent2 => node_map2%child3
    call node_map2%child3%cascade_adjacency(child3%parent,child3)
else
    child3%adjacent2 => node%adjacent2
end if 
if (associated(node_map3%child1)) then 
    child3%adjacent3 => node_map3%child1
    call node_map3%child1%cascade_adjacency(child3%parent,child3)
else
    child3%adjacent3 => node%adjacent3
end if 

!set child 4
node%child4 => child4
child4%level = node%level + 1
child4%parent => node  
child4%vertex1 = node%evertex2
child4%vertex2 = node%evertex3
child4%vertex3 = node%evertex1
child4%adjacent1 => child3
child4%adjacent2 => child1
child4%adjacent3 => child2
return 
end subroutine refine_node


!build oriented adjacent node map =========================
subroutine map_adjacent_node(self,adjnode,node_map)
implicit none 

!variables - import 
class(trinode), target :: self
type(trinode), pointer :: adjnode
type(trinode) :: node_map

!find relative orientation to build map
if (associated(adjnode)) then 
    if (associated(adjnode%adjacent1,self)) then 
        node_map%child1 => adjnode%child2
        node_map%child2 => adjnode%child3
        node_map%child3 => adjnode%child1
        node_map%aevertex1 => adjnode%evertex2
        node_map%aevertex2 => adjnode%evertex3
        node_map%aevertex3 => adjnode%evertex1
    elseif (associated(adjnode%adjacent2,self)) then 
        node_map%child1 => adjnode%child3
        node_map%child2 => adjnode%child1
        node_map%child3 => adjnode%child2
        node_map%aevertex1 => adjnode%evertex3
        node_map%aevertex2 => adjnode%evertex1
        node_map%aevertex3 => adjnode%evertex2
    elseif (associated(adjnode%adjacent3,self)) then 
        node_map%child1 => adjnode%child1
        node_map%child2 => adjnode%child2
        node_map%child3 => adjnode%child3
        node_map%aevertex1 => adjnode%evertex1
        node_map%aevertex2 => adjnode%evertex2
        node_map%aevertex3 => adjnode%evertex3
    end if 
else
    node_map%child1 => null()
    node_map%child2 => null()
    node_map%child3 => null()
    node_map%aevertex1 => null()
    node_map%aevertex2 => null()
    node_map%aevertex3 => null()
end if 
return 
end subroutine map_adjacent_node


!cascade adjacency update =========================
recursive subroutine cascade_adjacency(self,adj_tgt,adj_new)
implicit none

!variables - import 
class(trinode) :: self
type(trinode), pointer :: adj_tgt,adj_new

!variables - local 
logical :: updated

!update if present 
updated = .false.
if (associated(self%adjacent1,adj_tgt)) then 
    self%adjacent1 => adj_new
end if 
if (associated(self%adjacent2,adj_tgt)) then 
    self%adjacent2 => adj_new
end if 
if (associated(self%adjacent3,adj_tgt)) then 
    self%adjacent3 => adj_new
end if 

!call on children if an update is made
if (updated) then 
    if (associated(self%child1)) then 
        call self%child1%cascade_adjacency(adj_tgt,adj_new)
    end if 
    if (associated(self%child2)) then 
        call self%child2%cascade_adjacency(adj_tgt,adj_new)
    end if 
    if (associated(self%child3)) then 
        call self%child3%cascade_adjacency(adj_tgt,adj_new)
    end if 
end if 
return 
end subroutine cascade_adjacency   


end module hex_data