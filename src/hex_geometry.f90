!hex geometry module 
!max wood
!version : 0.0.3
!updated : 22-06-25

!module 
module hex_geometry
use hex_data_methods
use, intrinsic :: ieee_arithmetic

!routines 
contains 

!build point normal plane from two points (2d) ==============
subroutine get_point_normal_plane_from_two_points(vp1,vp2,pref,pnorm)
implicit none 

!variables - inout
real(dp) :: vp1(2),vp2(2),pref(2),pnorm(2)

!build the plane 
pref = 0.5d0*(vp1 + vp2)
pnorm(1) = vp2(2) - vp1(2)
pnorm(2) = vp1(1) - vp2(1)
pnorm = pnorm/(norm2(pnorm) + 1e-12_dp)
return
end subroutine get_point_normal_plane_from_two_points


!angle between two vectors ==============
function vec2vec_angle(vec1,vec2) result(angle)
implicit none 

!variables - import
real(dp) :: vec1(3),vec2(3)

!result
real(dp) :: angle

!evaluate (should be unit vectors ideally)
angle = atan2(norm2(cross_product(vec1,vec2)),dot_product(vec1,vec2))
return 
end function vec2vec_angle


!cross product ===========================
function cross_product(a,b) result(n)
implicit none

!variables - import
real(dp) :: n(3),a(3),b(3)

!cross product
n(1) = a(2)*b(3) - a(3)*b(2)
n(2) = -(a(1)*b(3) - a(3)*b(1))
n(3) = a(1)*b(2) - a(2)*b(1)
return
end function cross_product


!line segment 2d area contribution ==============
function asegment(v1,v2) result(Aseg) !+ve area for CCW oriented shapes
implicit none 

!variables - import
real(dp) :: v1(2),v2(2)

!result
real(dp) :: aseg

!segment area
aseg = 0.5d0*(v1(1)*v2(2) - v2(1)*v1(2))
return 
end function asegment


!fraction along line segment ==============
function fraction_along_line(v1,v2,vp) result(f)
implicit none 

!variables - inout
real(dp) :: f
real(dp) :: v1(2),v2(2),vp(2)

!evaluate fraction
f = dot_product((v2 - v1),(vp - v1))/(norm2((v2 - v1))**2)
return 
end function fraction_along_line


!is point within a line ==============
function is_in_line_segment(v1,v2,vp) result(in_line)
implicit none 

!variables - inout
logical :: in_line
real(dp) :: v1(2),v2(2),vp(2)

!variables - local
real(dp) :: f

!check fraction
f = fraction_along_line(v1,v2,vp)
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

!intersection denominator 
dval = (v1(1) - v2(1))*(v3(2) - v4(2)) - (v1(2) - v2(2))*(v3(1) - v4(1))

!test paralell
if (abs(dval) .NE. 0.0d0) then !if not parallel -> find intersection location
    t = ((v1(1) - v3(1))*(v3(2) - v4(2)) - (v1(2) - v3(2))*(v3(1) - v4(1)))/dval !L1 parameter
    vi(1) = (v1(1) + t*(v2(1) - v1(1)))
    vi(2) = (v1(2) + t*(v2(2) - v1(2)))
else !take as not intersecting as parallel  
    vi(:) = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if
return 
end function line_line_intersect


!line-line intersection (complex) -> L1 = v1 v2 | L2 = v3 v4 ==============
function line_line_intersect_cpx(v1,v2,v3,v4) result(vi)
implicit none 

!variables - inout
complex(dp) :: vi(2)
complex(dp) :: v1(2),v2(2),v3(2),v4(2)

!variables - local
complex(dp) :: dval,t

!intersection denominator 
dval = (v1(1) - v2(1))*(v3(2) - v4(2)) - (v1(2) - v2(2))*(v3(1) - v4(1))

!test paralell
if (abs(dval) .NE. 0.0d0) then !if not parallel -> find intersection location
    t = ((v1(1) - v3(1))*(v3(2) - v4(2)) - (v1(2) - v3(2))*(v3(1) - v4(1)))/dval !L1 parameter
    vi(1) = (v1(1) + t*(v2(1) - v1(1)))
    vi(2) = (v1(2) + t*(v2(2) - v1(2)))
else !take as not intersecting as parallel  
    vi(:) = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if
return 
end function line_line_intersect_cpx


!line-line intersection in lengths boolean ==============
function line_line_intersect_internal_bool(v1,v2,v3,v4) result(intersect)
implicit none 

!variables - inout
logical :: intersect 
real(dp) :: v1(2),v2(2),v3(2),v4(2)

!Variables - Local 
real(dp) :: vi(2)

!intersect lines 
vi = line_line_intersect(v1,v2,v3,v4) 
if (ieee_is_nan(vi(1))) then 
    intersect = .false.
else
    intersect = .false.
    if (is_in_line_segment(v1,v2,vi)) then 
        if (is_in_line_segment(v3,v4,vi)) then 
            intersect = .true.
        end if
    end if 
end if 
return 
end function line_line_intersect_internal_bool


!minimum distance point to line ==============
function min_dist_point_to_line(vl1,vl2,vp) result(vid)
implicit none 

!Variables - Import
real(dp) :: vp(2),vl1(2),vl2(2),vid(3)

!Variables - Local 
real(dp) :: t,f
real(dp) :: edir(2)

!Evaluate closest point 
edir = vl2 - vl1
edir = edir/norm2(edir)
t = dot_product(vp - vl1,edir)
vid(1:2) = vl1 + t*edir
f = dot_product(vid(1:2) - vl1,edir)/norm2(vl2 - vl1)
if (f .GT. 1.0d0) then 
    vid(1:2) = vl2
elseif (f .LT. 0.0d0) then 
    vid(1:2) = vl1
end if 
vid(3) = norm2(vid(1:2) - vp)
return 
end function min_dist_point_to_line
    

!closest point on a line to a point ==============
function closest_point_on_line_to_point(vl1,vl2,vp) result(vid)
implicit none 

!Variables - Import
real(dp) :: vp(2),vl1(2),vl2(2),vid(2)

!Variables - Local 
real(dp) :: t
real(dp) :: edir(2)

!evaluate
edir = vl2 - vl1
edir = edir/norm2(edir)
t = dot_product(vp - vl1,edir)
vid = vl1 + t*edir
return 
end function closest_point_on_line_to_point


!line-line intersection position gradient L1 = v1 v2 | L2 = v3 v4 (gradient wrt l1 vertices) ==============
! returns (x,y) gradient of the intersection location wrt the location of v1 and v2
subroutine line_line_intersect_gradient(v1,v2,v3,v4,grad_v1,grad_v2)
implicit none 

!variables - inout
real(dp) :: grad_v1(4),grad_v2(4) !dx x |dx y | dy x | dy y
real(dp) :: v1(2),v2(2),v3(2),v4(2)

!variables - local 
integer(in64) :: ii 
real(dp) :: h
real(dp) :: vgrad(2)
complex(dp) :: cpx_step
complex(dp) :: v1c(2),v2c(2),v3c(2),v4c(2),vp(2),vic(2)

!set the complex step stepsize
h = 1.0d0/10.0d0**40
cpx_step = cmplx(0.0d0,y=h,kind=dp)

!cast v1 v2 v3 v4 to complex
v1c(1) = cmplx(v1(1),y=0.0d0,kind=dp)
v1c(2) = cmplx(v1(2),y=0.0d0,kind=dp)
v2c(1) = cmplx(v2(1),y=0.0d0,kind=dp)
v2c(2) = cmplx(v2(2),y=0.0d0,kind=dp)
v3c(1) = cmplx(v3(1),y=0.0d0,kind=dp)
v3c(2) = cmplx(v3(2),y=0.0d0,kind=dp)
v4c(1) = cmplx(v4(1),y=0.0d0,kind=dp)
v4c(2) = cmplx(v4(2),y=0.0d0,kind=dp)

!grad_v1 x
vp = v1c
vp(1) = vp(1) + cpx_step
vic = line_line_intersect_cpx(vp,v2c,v3c,v4c)
vgrad(1) = aimag(vic(1))/h
vgrad(2) = aimag(vic(2))/h
grad_v1(1:2) = vgrad !dx | dy of intersection wrt moving v1 in x

!grad_v1 y
vp = v1c
vp(2) = vp(2) + cpx_step
vic = line_line_intersect_cpx(vp,v2c,v3c,v4c)
vgrad(1) = aimag(vic(1))/h
vgrad(2) = aimag(vic(2))/h
grad_v1(3:4) = vgrad !dx | dy of intersection wrt moving v1 in y

!grad_v2 x
vp = v2c
vp(1) = vp(1) + cpx_step
vic = line_line_intersect_cpx(v1c,vp,v3c,v4c)
vgrad(1) = aimag(vic(1))/h
vgrad(2) = aimag(vic(2))/h
grad_v2(1:2) = vgrad !dx | dy of intersection wrt moving v2 in x

!grad_v2 y
vp = v2c
vp(2) = vp(2) + cpx_step
vic = line_line_intersect_cpx(v1c,vp,v3c,v4c)
vgrad(1) = aimag(vic(1))/h
vgrad(2) = aimag(vic(2))/h
grad_v2(3:4) = vgrad !dx | dy of intersection wrt moving v2 in y

!set any nan values to zero 
do ii=1,4
    if (isnan(grad_v1(ii))) then 
        grad_v1(ii) = 0.0d0 
    end if 
    if (isnan(grad_v2(ii))) then 
        grad_v2(ii) = 0.0d0 
    end if 
end do 
return 
end subroutine line_line_intersect_gradient


!line triangle intersect boolean function ==============
function line_triangle_intersect_bool(vl1,vl2,vt1,vt2,vt3) result(intersect)
implicit none 

!variables - inout
logical :: intersect
real(dp) :: vl1(2),vl2(2),vt1(2),vt2(2),vt3(2)

!variables - local 
real(dp) :: at

!initialise 
intersect = .false.

!check intersects 
if (line_line_intersect_internal_bool(vl1,vl2,vt1,vt2)) then 
    intersect = .true.
    return 
elseif (line_line_intersect_internal_bool(vl1,vl2,vt2,vt3)) then 
    intersect = .true.
    return 
elseif (line_line_intersect_internal_bool(vl1,vl2,vt3,vt1)) then 
    intersect = .true.
    return 
end if 

!check containments 
at = triangle_area_2d(vt1,vt2,vt3)
if (is_in_triangle(vl1,vt1,vt2,vt3,at,-1.0d0)) then 
    intersect = .true.
    return 
elseif (is_in_triangle(vl2,vt1,vt2,vt3,at,-1.0d0)) then 
    intersect = .true.
    return 
end if 
return 
end function line_triangle_intersect_bool


!is in triangle bool ==============
function is_in_triangle(vp,vt1,vt2,vt3,at,sign_offset) result(inside)
implicit none 

!variables - inout
logical :: inside
real(dp) :: at,sign_offset
real(dp) :: vp(2),vt1(2),vt2(2),vt3(2)

!variables - local 
real(dp) :: vb(3)

!get barycentric coordinates
vb = cartesian2barycentric(vp,vt1,vt2,vt3,at)
vb = vb*sign_offset

!set containment 
inside = .false.
if ((vb(1) .GE. 0.0d0) .AND. (vb(1) .LE. 1.0d0)) then 
    if ((vb(2) .GE. 0.0d0) .AND. (vb(2) .LE. 1.0d0)) then 
        if ((vb(3) .GE. 0.0d0) .AND. (vb(3) .LE. 1.0d0)) then 
            inside = .true.
        end if 
    end if 
end if 
return 
end function is_in_triangle


!cartesian to barycentric ==============
function cartesian2barycentric(vp,vt1,vt2,vt3,at) result(vb)
implicit none 

!variables - inout
real(dp) :: at 
real(dp) :: vb(3)
real(dp) :: vp(2),vt1(2),vt2(2),vt3(2)

!evaluate
vb(1) = vt2(1)*vt3(2) - vt3(1)*vt2(1) + (vt2(2) - vt3(2))*vp(1) + (vt3(1) - vt2(1))*vp(2)
vb(2) = vt3(1)*vt1(2) - vt1(1)*vt3(1) + (vt3(2) - vt1(2))*vp(1) + (vt1(1) - vt3(1))*vp(2)
vb(3) = vt1(1)*vt2(2) - vt2(1)*vt1(1) + (vt1(2) - vt2(2))*vp(1) + (vt2(1) - vt1(1))*vp(2)
vb = vb/at
return 
end function cartesian2barycentric


!triangle area 2d ==============
function triangle_area_2d(vt1,vt2,vt3) result(at)
implicit none 

!variables - inout
real(dp) :: at 
real(dp) :: vt1(2),vt2(2),vt3(2)

!variables - local 
real(dp) :: a,b,c 

!side lengths
a = norm2(vt2 - vt1)
b = norm2(vt3 - vt2)
c = norm2(vt1 - vt3)

!area
at = 0.25d0*sqrt(4.0d0*a*a*b*b - (a*a + b*b - c*c)**2)
return 
end function triangle_area_2d


end module hex_geometry