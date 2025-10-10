!geometry kdtree module
!max wood
!version : 0.0.1
!updated : 17-01-25

!module
module geometry_kdtree
use halfedge_mesh 

!Integer data types 
use ISO_FORTRAN_ENV, only: in16=>int16
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in64=>int64

!Real data types 
use ISO_FORTRAN_ENV, only: sp=>real32 
use ISO_FORTRAN_ENV, only: dp=>real64
    
!items ==============
implicit none 

!edge pointer array 
type edgepointer
    type(edge), pointer :: edge
end type edgepointer

!face pointer array 
type facepointer
    type(face), pointer :: face
end type facepointer

!node type
type node
    integer(in64) :: nitem,depth,dim_div,index
    real(dp) :: median 
    real(dp), dimension(:), allocatable :: mini,maxi,box_min,box_max
    type(edgepointer), dimension(:), allocatable :: edge
    type(facepointer), dimension(:), allocatable :: face
    type(node), pointer :: parent,child1,child2
    contains 
        procedure :: add_item_2d
end type node 

!tree type 
type gkdtree
    logical :: display
    integer(in64) :: nnode 
    integer(in64) :: ndim_mesh,ndim_tree,depth,min_node_divsize
    real(dp), dimension(:,:), allocatable :: coordinates 
    type(node), dimension(:), allocatable :: node
    contains 
        procedure :: build_tree
        procedure :: nodes_overlapping_region_2d
end type gkdtree

!methods ==================================================
contains 


!build tree subroutine ==============
subroutine build_tree(self,mesh,maxdepth,min_node_divsize,display)
implicit none 

!variables - import 
logical :: display
integer(in64) :: maxdepth,min_node_divsize
class(gkdtree), target :: self
type(halfedge), target :: mesh

!variables - local 
integer(in64) :: ii,jj,rr,nn
integer(in64) :: depth,dimidx,nchild1,nchild2,nchildn,nnode0,nincomplete,nodecount
integer(in64), dimension(:), allocatable :: child_sort 
real(dp) :: avnodesize

!set dimension 
self%ndim_mesh = mesh%ndim
self%ndim_tree = 2*mesh%ndim

!build 2*mesh%ndim tree bounding box coordinates for items in the mesh 
if (self%ndim_mesh == 2) then 
    allocate(self%coordinates(mesh%nedge,4))
    do ii=1,mesh%nedge
        self%coordinates(ii,1) = min(mesh%edge(ii)%origin%coordinate(1),mesh%edge(ii)%opposite%origin%coordinate(1)) !xmin
        self%coordinates(ii,2) = min(mesh%edge(ii)%origin%coordinate(2),mesh%edge(ii)%opposite%origin%coordinate(2)) !ymin
        self%coordinates(ii,3) = max(mesh%edge(ii)%origin%coordinate(1),mesh%edge(ii)%opposite%origin%coordinate(1)) !xmin
        self%coordinates(ii,4) = max(mesh%edge(ii)%origin%coordinate(2),mesh%edge(ii)%opposite%origin%coordinate(2)) !ymin
    end do 
elseif (self%ndim_mesh == 3) then 

end if 

!initialise the tree
self%display = display
self%nnode = 1
self%min_node_divsize = min_node_divsize
allocate(self%node(mesh%nvertex*maxdepth))
do ii=1,mesh%nvertex*maxdepth
    self%node(ii)%nitem = 0 
    self%node(ii)%depth = -1 
    self%node(ii)%parent => null() 
    self%node(ii)%child1 => null() 
    self%node(ii)%child2 => null() 
    self%node(ii)%index = -1
end do 

!initialise the root node 
allocate(self%node(1)%edge(mesh%nedge))
do ii=1,mesh%nedge
    self%node(1)%edge(ii)%edge => mesh%edge(ii)
end do 
self%node(1)%nitem = mesh%nedge
allocate(self%node(1)%mini(self%ndim_tree))
allocate(self%node(1)%maxi(self%ndim_tree))
do ii=1,self%ndim_tree
    self%node(1)%mini(ii) = minval(self%coordinates(:,ii))
    self%node(1)%maxi(ii) = maxval(self%coordinates(:,ii))
end do 
allocate(self%node(1)%box_min(self%ndim_mesh))
allocate(self%node(1)%box_max(self%ndim_mesh))
do ii=1,self%ndim_mesh
    self%node(1)%box_min(ii) = self%node(1)%mini(ii)
    self%node(1)%box_max(ii) = self%node(1)%maxi(ii+self%ndim_mesh)
end do 
self%node(1)%depth = 0 
self%node(1)%index = 1

!construct the tree
depth = 1
if (self%ndim_mesh == 2) then 
    allocate(child_sort(mesh%nedge))
elseif (self%ndim_mesh == 3) then 
    
end if 
child_sort(:) = 0 
do rr=1,maxdepth*self%ndim_tree

    !target dimension index
    dimidx = 1 + mod(depth-1,self%ndim_tree)

    !divide each node 
    nchildn = 0 
    nincomplete = 0 
    nnode0 = self%nnode
    do nn=1,nnode0
        ! if ((self%node(nn)%depth == depth - 1) .AND. (self%node(nn)%nitem > min_node_divsize)) then 
        if ((.NOT. associated(self%node(nn)%child1)) .AND. (self%node(nn)%nitem > min_node_divsize)) then 

            !median of the is node in the target coordinate direction
            self%node(nn)%median = 0.0 
            do ii=1,self%node(nn)%nitem
                self%node(nn)%median = self%node(nn)%median + self%coordinates(self%node(nn)%edge(ii)%edge%index,dimidx)
            end do 
            self%node(nn)%median = self%node(nn)%median/real(self%node(nn)%nitem,dp)
            self%node(nn)%dim_div = dimidx

            !sort items into child nodes 
            nchild1 = 0
            nchild2 = 0 
            child_sort(1:self%node(nn)%nitem) = 0 
            do ii=1,self%node(nn)%nitem
                if (self%coordinates(self%node(nn)%edge(ii)%edge%index,dimidx) .GE. self%node(nn)%median) then 
                    nchild1 = nchild1 + 1
                    child_sort(ii) = 1
                else
                    nchild2 = nchild2 + 1
                    child_sort(ii) = 2
                end if 
            end do 

            !build chidren if valid divide
            if ((nchild1 > 0) .AND. (nchild2 > 0)) then 
                if (self%ndim_mesh == 2) then 

                    !build child node 1
                    self%nnode = self%nnode + 1
                    self%node(self%nnode)%parent => self%node(nn)
                    self%node(nn)%child1 => self%node(self%nnode)
                    self%node(nn)%child1%nitem = nchild1
                    self%node(nn)%child1%depth = depth
                    allocate(self%node(nn)%child1%edge(nchild1))
                    allocate(self%node(nn)%child1%mini(self%ndim_tree))
                    allocate(self%node(nn)%child1%maxi(self%ndim_tree))
                    nchild1 = 0 
                    do ii=1,self%node(nn)%nitem
                        if (child_sort(ii) == 1) then 
                            nchild1 = nchild1 + 1
                            self%node(nn)%child1%edge(nchild1)%edge => self%node(nn)%edge(ii)%edge
                        end if 
                    end do 
                    self%node(nn)%child1%mini(:) = huge(0.0d0)
                    self%node(nn)%child1%maxi(:) = -huge(0.0d0)
                    do ii=1,self%ndim_tree
                        do jj=1,self%node(nn)%child1%nitem
                            if (self%coordinates(self%node(nn)%child1%edge(jj)%edge%index,ii) .LT. &
                                self%node(nn)%child1%mini(ii)) then
                                self%node(nn)%child1%mini(ii) = self%coordinates(self%node(nn)%child1%edge(jj)%edge%index,ii)
                            end if 
                            if (self%coordinates(self%node(nn)%child1%edge(jj)%edge%index,ii) .GT. &
                                self%node(nn)%child1%maxi(ii)) then
                                self%node(nn)%child1%maxi(ii) = self%coordinates(self%node(nn)%child1%edge(jj)%edge%index,ii)
                            end if 
                        end do 
                    end do 
                    allocate(self%node(nn)%child1%box_min(self%ndim_mesh))
                    allocate(self%node(nn)%child1%box_max(self%ndim_mesh))
                    do ii=1,self%ndim_mesh
                        self%node(nn)%child1%box_min(ii) = self%node(nn)%child1%mini(ii)
                        self%node(nn)%child1%box_max(ii) = self%node(nn)%child1%maxi(ii+self%ndim_mesh)
                    end do 
                    self%node(nn)%child1%index = self%nnode
                    
                    !build child node 2
                    self%nnode = self%nnode + 1
                    self%node(self%nnode)%parent => self%node(nn)
                    self%node(nn)%child2 => self%node(self%nnode)
                    self%node(nn)%child2%nitem = nchild2
                    self%node(nn)%child2%depth = depth
                    allocate(self%node(nn)%child2%edge(nchild2))
                    allocate(self%node(nn)%child2%mini(self%ndim_tree))
                    allocate(self%node(nn)%child2%maxi(self%ndim_tree))
                    nchild2 = 0 
                    do ii=1,self%node(nn)%nitem
                        if (child_sort(ii) == 2) then 
                            nchild2 = nchild2 + 1
                            self%node(nn)%child2%edge(nchild2)%edge => self%node(nn)%edge(ii)%edge
                        end if 
                    end do 
                    self%node(nn)%child2%mini(:) = huge(0.0d0)
                    self%node(nn)%child2%maxi(:) = -huge(0.0d0)
                    do ii=1,self%ndim_tree
                        do jj=1,self%node(nn)%child2%nitem
                            if (self%coordinates(self%node(nn)%child2%edge(jj)%edge%index,ii) .LT. &
                                self%node(nn)%child2%mini(ii)) then
                                self%node(nn)%child2%mini(ii) = self%coordinates(self%node(nn)%child2%edge(jj)%edge%index,ii)
                            end if 
                            if (self%coordinates(self%node(nn)%child2%edge(jj)%edge%index,ii) .GT. &
                                self%node(nn)%child2%maxi(ii)) then
                                self%node(nn)%child2%maxi(ii) = self%coordinates(self%node(nn)%child2%edge(jj)%edge%index,ii)
                            end if 
                        end do 
                    end do 
                    allocate(self%node(nn)%child2%box_min(self%ndim_mesh))
                    allocate(self%node(nn)%child2%box_max(self%ndim_mesh))
                    do ii=1,self%ndim_mesh
                        self%node(nn)%child2%box_min(ii) = self%node(nn)%child2%mini(ii)
                        self%node(nn)%child2%box_max(ii) = self%node(nn)%child2%maxi(ii+self%ndim_mesh)
                    end do 
                    self%node(nn)%child2%index = self%nnode
                elseif (self%ndim_mesh == 3) then 


                end if 

                !count new children 
                nchildn = nchildn + 2
            end if 

            !count incomplete divisions **************************************
            nincomplete = nincomplete + 1
        end if 
    end do 

    !exit if no new nodes are constructed
    if (nchildn == 0) then 
        exit
    end if 

    !find average node size
    nodecount = 0 
    avnodesize = 0.0d0 
    do nn=1,self%nnode
        ! if (self%node(nn)%depth == depth) then 
        if (.NOT. associated(self%node(nn)%child1)) then 
            avnodesize = avnodesize + real(self%node(nn)%nitem,dp)
            nodecount = nodecount + 1
        end if 
    end do 
    avnodesize = avnodesize/real(nodecount,dp)

    !display
    if (self%display) then 
        write(*,'(A,I0,A,I0,A,I0,A,F0.4,A)') '    level -> ', rr ,' {dimension ',dimidx,'} nnode -> ',&
        self%nnode,' {av size ',avnodesize,'}'
    end if 
    
    !increment the tree depth 
    depth = depth + 1
end do 

!tag each edge with its leaf node
if (self%ndim_mesh == 2) then 
    do nn=1,self%nnode
        if (.NOT. associated(self%node(nn)%child1)) then 
            do ii=1,self%node(nn)%nitem
                mesh%edge(self%node(nn)%edge(ii)%edge%index)%tag = nn 
            end do 
        end if 
    end do 
elseif (self%ndim_mesh == 3) then 


end if 
return 
end subroutine build_tree


!add items to node ==============
subroutine add_item_2d(self,nedge)
implicit none 

!variables - import 
class(node) :: self
type(edge), target :: nedge

!variables - local 
integer(in64) :: ii 
type(edgepointer), dimension(:), allocatable :: etemp

!add item 
if (self%nitem + 1 .GT. size(self%edge,dim = 1)) then 
    allocate(etemp(self%nitem))
    do ii=1,self%nitem
        etemp(ii)%edge => self%edge(ii)%edge 
    end do 
    deallocate(self%edge)
    allocate(self%edge(2*self%nitem))
    do ii=1,self%nitem
        self%edge(ii)%edge => etemp(ii)%edge
    end do 
end if 
self%nitem = self%nitem + 1
self%edge(self%nitem)%edge => nedge
return 
end subroutine add_item_2d


!get nodes overlapping region 2d subroutine ==============
subroutine nodes_overlapping_region_2d(self,xmin,xmax,ymin,ymax,ntgtnode,tgt_nodes)
implicit none 

!variables - import 
integer(in64) :: ntgtnode
real(dp) :: xmin,xmax,ymin,ymax
class(gkdtree), target :: self
type(node), dimension(:) :: tgt_nodes

!variables - local 
logical :: overlap 
integer(in64) :: nn,ii 
integer(in64) :: nsearch,nsearchN,ntgt
integer(in64) :: node_search(self%nnode),node_searchN(self%nnode)

!search for nodes overlapping the region xmin,xmax,ymin,ymax
ntgtnode = 0
nsearch = 1
node_search(1) = 1 
do nn=1,self%nnode 

    !search current nodes 
    nsearchN = 0 
    do ii=1,nsearch

        !target node
        ntgt = node_search(ii)

        !check if overlap 
        overlap = .false.
        if ((self%node(ntgt)%box_max(1) .GE. xmin) .AND. (self%node(ntgt)%box_min(1) .LE. xmax)) then 
            if ((self%node(ntgt)%box_max(2) .GE. ymin) .AND. (self%node(ntgt)%box_min(2) .LE. ymax)) then 
                overlap = .true.
            end if 
        end if 

        !process if overlap 
        if (overlap) then 
            if (.NOT.associated(self%node(ntgt)%child1)) then !has no children so return this node
                ntgtnode = ntgtnode + 1
                tgt_nodes(ntgtnode) = self%node(ntgt)
            else !add children to new search list
                nsearchN = nsearchN + 1
                node_searchN(nsearchN) = self%node(ntgt)%child1%index   
                nsearchN = nsearchN + 1
                node_searchN(nsearchN) = self%node(ntgt)%child2%index  
            end if
        end if 
    end do

    !exit if no new nodes found 
    if (nsearchN == 0) then
        exit 
    end if 

    !update search list
    nsearch = nsearchN
    node_search(1:nsearchN) = node_searchN(1:nsearchN)
end do 
return 
end subroutine nodes_overlapping_region_2d


end module geometry_kdtree