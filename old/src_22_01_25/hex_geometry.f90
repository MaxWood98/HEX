!hex geometry module 
!max wood
!version : 0.0.1
!updated : 19-01-25

!module 
module hex_geometry
use hex_data
use, intrinsic :: ieee_arithmetic

!routines 
contains 


!is point within a line segment ==============
function is_in_line_segment(v1,v2,vp) result(in_line)
implicit none 

!variables - inout
logical :: in_line
real(dp) :: v1(2),v2(2),vp(2)

!variables - local
real(dp) :: f

!check fraction
f = dot_product((v2 - v1),(vp - v1))/(norm2((v2 - v1))**2)
if ((f .GE. 0.0d0) .AND. (f .LE. 1.0d0)) then 
    in_line = .true.
else
    in_line = .false.
end if 
return 
end function is_in_line_segment


!line-line intersection -> L1 = v1 v2 | L2 = v3 v4 ==============
function line_line_intersect(v1,v2,v3,v4) result(vi)
implicit none 

!variables - inout
real(dp) :: vi(2)
real(dp) :: v1(2),v2(2),v3(2),v4(2)

!variables - local
real(dp) :: dval,t

!Intersection denominator 
dval = (v1(1) - v2(1))*(v3(2) - v4(2)) - (v1(2) - v2(2))*(v3(1) - v4(1))

!Test paralell
if (dval .NE. 0.0d0) then !If not parallel -> find intersection location
    t = ((v1(1) - v3(1))*(v3(2) - v4(2)) - (v1(2) - v3(2))*(v3(1) - v4(1)))/dval !L1 parameter
    vi(1) = (v1(1) + t*(v2(1) - v1(1)))
    vi(2) = (v1(2) + t*(v2(2) - v1(2)))
else !Take as not intersecting as parallel  
    vi(:) = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if
return 
end function line_line_intersect


end module hex_geometry