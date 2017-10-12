! Mar 2009 D.Pietrobon

!!$ subroutine find_index(n, x, xm, ixm)
!!$
!!$   implicit none
!!$
!!$   integer, intent(in) :: n
!!$   integer, intent(out) :: ixm
!!$   double precision, dimension(n), intent(in) :: x
!!$   double precision, intent(in) :: xm
!!$   integer i
!!$   double precision, dimension(n) :: ii
!!$
!!$!real :: t1, t2
!!$
!!$   integer imask(n)
!!$!call cpu_time(t1)
!!$
!!$   imask = 0.
!!$
!!$   ii = abs(x - xm)
!!$   where (ii .EQ. minval(ii) ) imask = 1.
!!$   do i = 1,n 
!!$      ii(i) = i
!!$   enddo
!!$   ixm = dot_product(imask, ii)
!!$
!!$!call cpu_time(t2)
!!$!print*, 'find_index', t2-t1
!!$
!!$   return
!!$
!!$ end subroutine find_index
 
! ---------------------------
!!$ subroutine grhoDE_splint(xa, ya, y2a, n, x, y)
 subroutine grhoDE_splint(x, y)

! Original NR
   use ModelParams, only: ailn, afln, nstep_D, grhoDE_lna, lna_i, lna_f, dlna
   use Precision

   implicit none

!   integer n
!   double precision x, y, xa(n), ya(n), y2a(n),
   real(dl), INTENT(in) :: x
   real(dl), INTENT(out) :: y
   real(dl) ::  xx!,lna_i, lna_f, dlna
   integer k, khi, klo
   real(dl) a, b, h

   khi = nstep_D
   klo = 1 

!   lna_i = ailn * LOG(10.d0)
!   lna_f = afln * LOG(10.d0)
!   dlna = (lna_f - lna_i) / (nstep_D-1)

   xx = abs(x - lna_i) / dlna

   k = INT(xx)

!if (x == 0.) print*,'k =',k

   if (grhoDE_lna(k,1) .gt. x) then
      khi = k
      klo = k-1
   else
      klo = k
      khi = k+1
   endif

   h = grhoDE_lna(khi,1)-grhoDE_lna(klo,1)

   if (h .eq.0) then 
      print*, 'step h vanishing in mysplint', k
      pause
   endif

   a = (grhoDE_lna(khi,1)-x)/h
   b = (x-grhoDE_lna(klo,1))/h

   y = a*grhoDE_lna(klo,2) + b*grhoDE_lna(khi,2) + ((a**3-a)*grhoDE_lna(klo,3) + (b**3-b)*grhoDE_lna(khi,3)) * (h**2)/6.

   return

 end subroutine grhoDE_splint


! ---------------------------
