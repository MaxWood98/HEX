!face-vertex mesh structure main class
!max wood
!version : 0.0.3
!updated : 09-02-25

!module
module facevertex_mesh
use io_utilities

!integer data types 
use ISO_FORTRAN_ENV, only: in16=>int16
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in64=>int64

!real data types 
use ISO_FORTRAN_ENV, only: sp=>real32 
use ISO_FORTRAN_ENV, only: dp=>real64

!items 
implicit none 

!constants 
integer(in64), parameter :: face_max_edges_fv = 10000

!facevertex face type
type face_fv
    integer(in64) :: nvertex 
    integer(in64), dimension(:), allocatable :: vertices,edges
end type face_fv

!facevertex type 
type facevertex 
    integer(in64) :: nvertex,nedge,nface,ndim,maxvalence
    integer(in64), dimension(:), allocatable :: valence
    integer(in64), dimension(:,:), allocatable :: edges
    integer(in64), dimension(:,:), allocatable :: v2v,v2e,v2f,e2f
    real(dp), dimension(:,:), allocatable :: vertices
    type(face_fv), dimension(:), allocatable :: faces
    contains 
        procedure :: read_fv_file
        procedure :: build_edges_from_faces 
        procedure :: get_valence_2d
        procedure :: get_valence_3d
        procedure :: get_valence
        procedure :: build_v2v_v2e
        procedure :: build_v2f
        procedure :: build_f2e
        procedure :: build_e2f
        procedure :: get_connectivity
        procedure :: build_edge_loop_faces
end type facevertex 


!methods ==================================================
contains 


!read .fv mesh file =========================
subroutine read_fv_file(self,filename)
implicit none 

!variables - import
class(facevertex) :: self
character(*), intent(in) :: filename

!variables - local
integer(in32) :: fh 
integer(in64) :: ii,jj
integer(in64) :: iostatus,idxsp1
real(dp) :: edgeline(2)
character(len=face_max_edges_fv) :: rtemp 

!check if file exists 
if (.NOT. file_exists(filename)) then 
    write(*,'(A)') '** cannot locate geometry file: '//trim(filename)
    stop 
end if 

!open mesh file 
fh = 11
open(fh,file=filename)

!initialise maxvalence 
self%maxvalence = 0 

!read number of dimensions
self%ndim = 0 
iostatus = 0 
do while (iostatus == 0)
    read(fh,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:4) == 'ndim') then 
        read(rtemp(6:len_trim(rtemp)),*) self%ndim
        exit 
    end if 
end do
rewind(fh)

!read vertices 
self%nvertex = 0 
iostatus = 0 
do while (iostatus == 0)
    read(fh,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:7) == 'nvertex') then 
        read(rtemp(9:len_trim(rtemp)),*) self%nvertex
        allocate(self%vertices(self%nvertex,3)) !force always to 3d (in 2d the z value is set to zero)
        self%vertices(:,:) = 0.0d0 
        do ii=1,self%nvertex
            read(fh,*) self%vertices(ii,1:self%ndim) 
        end do 
        exit 
    end if 
end do
rewind(fh)

!read edges 
self%nedge = 0 
iostatus = 0 
do while (iostatus == 0)
    read(fh,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:5) == 'nedge') then 
        read(rtemp(7:len_trim(rtemp)),*) self%nedge
        if (self%nedge .GT. 0) then 
            allocate(self%edges(self%nedge,2))
            do ii=1,self%nedge
                read(fh,*) edgeline 
                self%edges(ii,:) = int(edgeline,in64)
            end do 
        end if 
        exit 
    end if 
end do
rewind(fh)

!read faces 
self%nface = 0 
iostatus = 0 
do while (iostatus == 0)
    read(fh,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:5) == 'nface') then 
        read(rtemp(7:len_trim(rtemp)),*) self%nface
        if (self%nface .GT. 0) then 
            allocate(self%faces(self%nface))
            do ii=1,self%nface
                read(fh,'(A)') rtemp 
                idxsp1 = 0 
                do jj=1,face_max_edges_fv
                    if (rtemp(jj:jj) == ' ') then 
                        read(rtemp(1:jj-1),*) self%faces(ii)%nvertex 
                        idxsp1 = jj
                        exit 
                    end if 
                end do 
                allocate(self%faces(ii)%vertices(self%faces(ii)%nvertex))
                read(rtemp(idxsp1+1:len_trim(rtemp)),*) self%faces(ii)%vertices(:)
            end do 
        end if 
        exit 
    end if 
end do

!close file 
close(fh)
return 
end subroutine read_fv_file


!build edges from faces subroutine =========================
subroutine build_edges_from_faces(self)
implicit none 

!variables - Import
class(facevertex) :: self

!variables - Local 
integer(in64) :: ff,ee,vv 
integer(in64) :: ev1,ev2,evalid,edge_idx
integer(in64), dimension(:,:), allocatable :: vconnect,edgeidx 

!evaluate valence
call self%get_valence()

!initialise
allocate(vconnect(self%nvertex,self%maxvalence))
allocate(edgeidx(self%nvertex,self%maxvalence))
vconnect(:,:) = 0
edgeidx(:,:) = 0

!index edges 
self%nedge = 0 
do ff=1,self%nface
    do ee=1,self%faces(ff)%nvertex

        !edge end vertices
        ev1 = ee
        ev2 = mod(ee,self%faces(ff)%nvertex) + 1
        ev1 = self%faces(ff)%vertices(ev1)
        ev2 = self%faces(ff)%vertices(ev2)

        !check against vconnect 
        evalid = 1 
        do vv=1,self%maxvalence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit 
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !add if valid 
        if (evalid == 1) then 

            !increment edge count 
            self%nedge = self%nedge + 1

            !add edge to connection structure     
            do vv=1,self%maxvalence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    edgeidx(ev1,vv) = self%nedge
                    exit 
                end if 
            end do 
            do vv=1,self%maxvalence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    edgeidx(ev2,vv) = self%nedge
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!build edges 
if (allocated(self%edges)) then 
    deallocate(self%edges)
end if 
allocate(self%edges(self%nedge,2))
self%edges(:,:) = 0 
do ff=1,self%nface
    do ee=1,self%faces(ff)%nvertex

        !edge end vertices
        ev1 = ee
        ev2 = mod(ee,self%faces(ff)%nvertex) + 1
        ev1 = self%faces(ff)%vertices(ev1)
        ev2 = self%faces(ff)%vertices(ev2)

        !check index of this edge 
        edge_idx = 0 
        do vv=1,self%maxvalence
            if (vconnect(ev1,vv) == ev2) then 
                edge_idx = edgeidx(ev1,vv)
                exit 
            end if 
        end do 

        !construct
        if (edge_idx .GT. 0) then !if valid index then remove from edgeidx and construct 

            !remove edge index 
            do vv=1,self%maxvalence
                if (edgeidx(ev1,vv) == edge_idx) then 
                    edgeidx(ev1,vv) = -1*edgeidx(ev1,vv)
                    exit 
                end if 
            end do 
            do vv=1,self%maxvalence
                if (edgeidx(ev2,vv) == edge_idx) then 
                    edgeidx(ev2,vv) = -1*edgeidx(ev2,vv)
                    exit 
                end if 
            end do 

            !build edge 
            self%edges(edge_idx,1) = ev1
            self%edges(edge_idx,2) = ev2
        end if 
    end do 
end do 
return 
end subroutine build_edges_from_faces


!get valence subroutine =========================
subroutine get_valence(self)
implicit none 

!variables - import
class(facevertex) :: self

!switch on ndim 
if (self%ndim == 2) then 
    call self%get_valence_2d()
elseif (self%ndim == 3) then 
    call self%get_valence_3d()
end if 
return 
end subroutine get_valence


!get valence 2d subroutine =========================
subroutine get_valence_2d(self)
implicit none 

!variables - import
class(facevertex) :: self

!variables - local
integer(in64) :: ee

!initialse valence array
if (allocated(self%valence)) then 
    deallocate(self%valence)
end if 
allocate(self%valence(self%nvertex))
self%valence(:) = 0 

!evaluate valence 
do ee=1,self%nedge
    self%valence(self%edges(ee,1)) = self%valence(self%edges(ee,1)) + 1
    self%valence(self%edges(ee,2)) = self%valence(self%edges(ee,2)) + 1
end do 

!set maximum valence
self%maxvalence = maxval(self%valence(:))
return 
end subroutine get_valence_2d


!get valence 3d subroutine =========================
subroutine get_valence_3d(self)
implicit none 

!variables - import
class(facevertex) :: self

!variables - local
integer(in64) :: ff,ee,vv
integer(in64) :: ev1,ev2,ubValence,evalid
integer(in64), dimension(:,:), allocatable :: vconnect

!initialse valence array
if (allocated(self%valence)) then 
    deallocate(self%valence)
end if 
allocate(self%valence(self%nvertex)) 

!upper bound of maximum valence 
self%valence(:) = 0 
do ff=1,self%nface
    do ee=1,self%faces(ff)%nvertex

        !edge end vertices
        ev1 = ee
        ev2 = mod(ee,self%faces(ff)%nvertex) + 1
        ev1 = self%faces(ff)%vertices(ev1) 
        ev2 = self%faces(ff)%vertices(ev2) 

        !accumulate valence 
        self%valence(ev1) = self%valence(ev1) + 1
        self%valence(ev2) = self%valence(ev2) + 1
    end do 
end do 
ubValence = 4*maxval(self%valence)

!construct actual valence of each vertex
allocate(vconnect(self%nvertex,ubValence))
vconnect(:,:) = 0 
self%valence(:) = 0
do ff=1,self%nface
    do ee=1,self%faces(ff)%nvertex

        !edge end vertices
        ev1 = ee
        ev2 = mod(ee,self%faces(ff)%nvertex) + 1
        ev1 = self%faces(ff)%vertices(ev1) 
        ev2 = self%faces(ff)%vertices(ev2)  

        !check against vconnect 
        evalid = 1
        do vv=1,ubValence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !add valence if new edge
        if (evalid == 1) then 
            
            !increment valence on each vertex
            self%valence(ev1) = self%valence(ev1) + 1
            self%valence(ev2) = self%valence(ev2) + 1

            !update vconnect
            do vv=1,ubValence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    exit 
                end if 
            end do 
            do vv=1,ubValence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!set maximum valence
self%maxvalence = maxval(self%valence(:))
return 
end subroutine get_valence_3d


!get connectivity subroutine =========================
subroutine get_connectivity(self)
implicit none 

!variables - import
class(facevertex) :: self

!get valence
call self%get_valence()

!build v2v and v2e
call self%build_v2v_v2e()

!build v2f
call self%build_v2f()

!build f2e
call self%build_f2e()

!build e2f 
call self%build_e2f()
return 
end subroutine get_connectivity


!build v2v and v2e =========================
subroutine build_v2v_v2e(self)
implicit none 

!variables - import
class(facevertex) :: self

!variables - local 
integer(in64) :: ee
integer(in64) :: v1,v2
integer(in64) :: vins(self%nvertex)

!build v2v and v2e
vins(:) = 0 
allocate(self%v2v(self%nvertex,self%maxvalence))
allocate(self%v2e(self%nvertex,self%maxvalence))
self%v2v(:,:) = 0 
self%v2e(:,:) = 0 
do ee=1,self%nedge

    !edge vertices 
    v1 = self%edges(ee,1)
    v2 = self%edges(ee,2)

    !increment counts
    vins(v1) = vins(v1) + 1
    vins(v2) = vins(v2) + 1

    !add to v2v
    self%v2v(v1,vins(v1)) = v2
    self%v2v(v2,vins(v2)) = v1

    !add to v2e
    self%v2e(v1,vins(v1)) = ee
    self%v2e(v2,vins(v2)) = ee
end do
return 
end subroutine build_v2v_v2e


!Build v2f subroutine =========================
subroutine build_v2f(self)
implicit none 

!variables - import
class(facevertex) :: self

!variables - local 
integer(in64) :: ff,vv
integer(in64) :: v1
integer(in64) :: vins(self%nvertex)

!build v2f
vins(:) = 0 
allocate(self%v2f(self%nvertex,self%maxvalence))
self%v2f(:,:) = 0 
do ff=1,self%nface
    do vv=1,self%faces(ff)%nvertex
        v1 = self%faces(ff)%vertices(vv)
        vins(v1) = vins(v1) + 1
        self%v2f(v1,vins(v1)) = ff
    end do 
end do 
return
end subroutine build_v2f


!Build f2e subroutine =========================
subroutine build_f2e(self)
implicit none 

!variables - import
class(facevertex) :: self

!variables - local 
integer(in64) :: ff,vv,aa
integer(in64) :: v1,v2,etgt

!build f2e
do ff=1,self%nface
    allocate(self%faces(ff)%edges(self%faces(ff)%nvertex))
    self%faces(ff)%edges(:) = 0 
    do vv=1,self%faces(ff)%nvertex

        !end vertices of this edge 
        v1 = vv 
        v2 = mod(vv,self%faces(ff)%nvertex) + 1
        v1 = self%faces(ff)%vertices(v1)
        v2 = self%faces(ff)%vertices(v2)

        !search v2e for connection 
        do aa=1,self%valence(v1)
            etgt = self%v2e(v1,aa)
            if ((self%edges(etgt,1) == v2) .OR. (self%edges(etgt,2) == v2)) then
                self%faces(ff)%edges(vv) = etgt 
                exit 
            end if 
        end do 
    end do 
end do 
return 
end subroutine build_f2e


!build e2f subroutine =========================
subroutine build_e2f(self)
implicit none 

!variables - import
class(facevertex) :: self

!variables - Local 
integer(in64) :: ff,vv
integer(in64) :: etgt
integer(in64) :: eins(self%nedge)

!build e2f 
eins(:) = 0 
allocate(self%e2f(self%nedge,2))
self%e2f(:,:) = 0 
do ff=1,self%nface
    do vv=1,self%faces(ff)%nvertex
        etgt = self%faces(ff)%edges(vv)
        eins(etgt) = eins(etgt) + 1
        self%e2f(etgt,eins(etgt)) = ff
    end do 
end do 
return 
end subroutine build_e2f


!build edge loop faces =========================
subroutine build_edge_loop_faces(self)
implicit none 

!variables - import
class(facevertex) :: self

!variables - Local 
integer(in64) :: ii,ff
integer(in64) :: nloop,edge0,edgeC
integer(in64) :: v2e(self%nvertex,2),edge_tag(self%nedge),nedge_loop(self%nedge)

!build edge-vertex connectivity 
v2e(:,:) = 0 
do ii=1,self%nedge
    v2e(self%edges(ii,1),2) = ii 
    v2e(self%edges(ii,2),1) = ii 
end do 

!flood edge loops 
nloop = 0
edge_tag(:) = 0 
nedge_loop(:) = 0 
do ff=1,self%nedge

    !find new edge
    edge0 = 0
    do ii=1,self%nedge
        if (edge_tag(ii) == 0) then 
            edge0 = ii
            exit 
        end if 
    end do 

    !exit if no new edges 
    if (edge0 == 0) then 
        exit 
    end if 

    !increment loop 
    nloop = nloop + 1

    !flood loop to build face
    edgeC = edge0
    do ii=1,self%nedge
        edge_tag(edgeC) = nloop
        nedge_loop(nloop) = nedge_loop(nloop) + 1
        edgeC = v2e(self%edges(edgeC,2),2)
        if (edgeC == edge0) then 
            exit 
        end if 
    end do 
end do 

!build faces
self%nface = nloop
allocate(self%faces(nloop))
do ff=1,nloop
    self%faces(ff)%nvertex = nedge_loop(ff)
    allocate(self%faces(ff)%vertices(self%faces(ff)%nvertex))
    edge0 = 0
    do ii=1,self%nedge
        if (edge_tag(ii) == ff) then 
            edge0 = ii
            exit 
        end if 
    end do 
    edgeC = edge0
    do ii=1,self%nedge
        self%faces(ff)%vertices(ii) = self%edges(edgeC,1)
        edgeC = v2e(self%edges(edgeC,2),2)
        if (edgeC == edge0) then 
            exit 
        end if 
    end do 
end do 
return 
end subroutine build_edge_loop_faces


end module facevertex_mesh