!hex mesh generator utilities module 
!max wood
!version : 0.0.1
!updated : 10-10-25

!module 
module hex_utilities
use halfedge_mesh
use facevertex_mesh
use hex_data_methods
contains 

!set 2d harfedge edge normals from the corresponding facevertex =========================
subroutine set_2d_halfedge_normals_from_facevertex(geometry_he,geometry_fv)
implicit none

!variables - import 
type(halfedge) :: geometry_he
type(facevertex) :: geometry_fv

!variables local 
integer(in64) :: ii 
integer(in64) :: v1,v2
real(dp) :: fv_normals(geometry_fv%nedge,2) 

!get the facevertex edge normals 
do ii=1,geometry_fv%nedge
    v1 = geometry_fv%edges(ii,1)
    v2 = geometry_fv%edges(ii,2)
    fv_normals(ii,1) = geometry_fv%vertices(v1,2) - geometry_fv%vertices(v2,2)
    fv_normals(ii,2) = geometry_fv%vertices(v2,1) - geometry_fv%vertices(v1,1)
    fv_normals(ii,:) = fv_normals(ii,:)/norm2(fv_normals(ii,:))
end do 

!set these on the halfedge mesh
do ii=1,geometry_he%nedge
    geometry_he%edge(ii)%normal(1) = fv_normals(geometry_he%edge(ii)%tag,1)
    geometry_he%edge(ii)%normal(2) = fv_normals(geometry_he%edge(ii)%tag,2)
    geometry_he%edge(ii)%normal(3) = 0.0d0 
end do 
return 
end subroutine set_2d_halfedge_normals_from_facevertex

end module hex_utilities
