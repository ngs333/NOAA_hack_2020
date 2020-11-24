module fyppm_mod

 implicit none
 private
 public fyppm

 integer :: grid_type = 2

 real, parameter:: r3 = 1./3.
 real, parameter:: near_zero = 1.E-25

#ifdef WAVE_FORM
! Suresh & Huynh scheme 2.2 (purtabation form)
! The wave-form is more diffusive than scheme 2.1
 real, parameter:: b1 =   0.0375
 real, parameter:: b2 =  -7./30.
 real, parameter:: b3 =  -23./120.
 real, parameter:: b4 =  13./30.
 real, parameter:: b5 = -11./240.
#else
! scheme 2.1: perturbation form
 real, parameter:: b1 =   1./30.
 real, parameter:: b2 = -13./60.
 real, parameter:: b3 = -13./60.
 real, parameter:: b4 =  0.45
 real, parameter:: b5 = -0.05
#endif

 real, parameter:: t11 = 27./28., t12 = -13./28., t13=3./7.
 real, parameter:: s11 = 11./14., s14 = 4./7.,    s15=3./14.
!----------------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
!----------------------------------------------------
! Non-monotonic
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.

 contains

 subroutine fyppm(c,  q,  flux, jord, ifirst, ilast, jfirst, jlast, npy, dm, ppm_limiter, &
                  dya, isd, ied, jsd, jed, js, je, ng, km, al, bl, br, dq )
 integer, INTENT(IN) :: jfirst, jlast, ifirst, ilast
 integer, INTENT(IN) :: isd, ied
 integer, INTENT(IN) :: jsd, jed
 integer, INTENT(IN) :: js, je
 integer, INTENT(IN) :: ng
 integer, INTENT(IN) :: jord
 integer, INTENT(IN) :: npy, km
 real   , intent(IN) :: ppm_limiter
 real   , INTENT(IN) :: q(ifirst:ilast,jfirst-ng:jlast+ng,km)
 real   , intent(in) :: c(isd:ied,js:je+1,km )                 ! Courant number
 real   , INTENT(OUT):: flux(ifirst:ilast,jfirst:jlast+1,km)   ! Flux
 real   , INTENT(OUT):: dm(ifirst:ilast,jfirst-2:jlast+2,km)
 real   , intent(in) :: dya(isd:ied,js:je)
 real   , intent(out):: al(ifirst:ilast,jfirst-1:jlast+2,km)
 real   , intent(out):: bl(ifirst:ilast,jfirst-1:jlast+1,km)
 real   , intent(out):: br(ifirst:ilast,jfirst-1:jlast+1,km)
 real   , intent(out):: dq(ifirst:ilast,jfirst-3:jlast+2,km)


! Local:
 real dl, dr, pmp, lac, ct, qe
 real pmp_1, lac_1, pmp_2, lac_2
 real xt, x0, x1
 integer i, j, js3, je3, jt, k

 if( jord<=10 ) then    ! jord=8, 9, 10
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
    do k = 1, km
       do j=js-2,je+2
          do i=ifirst,ilast
             xt = 0.25*(q(i,j+1,k) - q(i,j-1,k))
             dm(i,j,k) = sign(min(abs(xt), max(q(i,j-1,k), q(i,j,k), q(i,j+1,k)) - q(i,j,k),   &
                  q(i,j,k) - min(q(i,j-1,k), q(i,j,k), q(i,j+1,k))), xt)
          enddo
       enddo
    enddo
!$ACC end kernels
    if (grid_type < 3) then
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
       do k = 1, km
          do j=max(3,js-1),min(npy-2,je+2)
             do i=ifirst,ilast
                al(i,j,k) = 0.5*(q(i,j-1,k)+q(i,j,k)) + r3*(dm(i,j-1,k) - dm(i,j,k))
             enddo
          enddo
       enddo
!$ACC end kernels
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
       do k = 1, km
          do j=js-3,je+2
             do i=ifirst,ilast
                dq(i,j,k) = q(i,j+1,k) - q(i,j,k)
             enddo
          enddo
       enddo
!$ACC end kernels
       if ( jord==8 ) then
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=max(3,js-1),min(npy-3,je+1)
                do i=ifirst,ilast
                   xt = 2.*dm(i,j,k)
                   bl(i,j,k) = -sign(min(abs(xt), abs(al(i,j,k)-q(i,j,k))),   xt)
                   br(i,j,k) =  sign(min(abs(xt), abs(al(i,j+1,k)-q(i,j,k))), xt)
                enddo
             enddo
          enddo
!$ACC end kernels
       elseif( jord==9 ) then
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=max(3,js-1),min(npy-3,je+1)
                do i=ifirst,ilast
                   pmp_1 = -2.*dq(i,j,k)
                   lac_1 = pmp_1 + 1.5*dq(i,j+1,k)
                   bl(i,j,k) = min(max(0., pmp_1, lac_1), max(al(i,j,k)-q(i,j,k), min(0., pmp_1, lac_1)))
                   pmp_2 = 2.*dq(i,j-1,k)
                   lac_2 = pmp_2 - 1.5*dq(i,j-2,k)
                   br(i,j,k) = min(max(0., pmp_2, lac_2), max(al(i,j+1,k)-q(i,j,k), min(0., pmp_2, lac_2)))
                enddo
             enddo
          enddo
!$ACC end kernels
       else    ! jord=10
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=max(3,js-1),min(npy-3,je+1)
                do i=ifirst,ilast
                   bl(i,j,k) = al(i,j,k) - q(i,j,k)
                   br(i,j,k) = al(i,j+1,k) - q(i,j,k)
                   if ( abs(dm(i,j-1,k))+abs(dm(i,j,k))+abs(dm(i,j+1,k)) < near_zero ) then
                      bl(i,j,k) = 0.
                      br(i,j,k) = 0.
                   elseif( abs(3.*(bl(i,j,k)+br(i,j,k))) > abs(bl(i,j,k)-br(i,j,k)) ) then
                      pmp_1 = -2.*dq(i,j,k)
                      lac_1 = pmp_1 + 1.5*dq(i,j+1,k)
                      bl(i,j,k) = min(max(0.,pmp_1,lac_1), max(bl(i,j,k), min(0.,pmp_1,lac_1)))
                      pmp_2 = 2.*dq(i,j-1,k)
                      lac_2 = pmp_2 - 1.5*dq(i,j-2,k)
                      br(i,j,k) = min(max(0.,pmp_2,lac_2), max(br(i,j,k), min(0.,pmp_2,lac_2)))
                   endif
                enddo
             enddo
          enddo
!$ACC end kernels
       endif

          !--------------
          ! Fix the edges:
          !--------------
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
       if( js==1 ) then
           do k = 1, km
             do i=ifirst,ilast
                br(i,2,k) = al(i,3,k) - q(i,2,k)
                !           xt = t11*(q(i,0)+q(i,1)) + t12*(q(i,-1)+q(i,2))   &
                !                                  + t13*(dm(i,2)-dm(i,-1))
                xt = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0,k)+q(i,1,k))   &
                     -dya(i,1)*(q(i,-1,k)+q(i,2,k))) / ( dya(i,1)+dya(i,2) )
                bl(i,1,k) = xt - q(i,1,k)
                br(i,0,k) = xt - q(i,0,k)
                xt = s14*dm(i,-1,k) - s11*dq(i,-1,k) + q(i,0,k)

                !           xt = min( xt, max(q(i,-1), q(i,0)) )
                !           xt = max( xt, min(q(i,-1), q(i,0)) )

                bl(i,0,k) = xt - q(i,0,k)
                xt = s15*q(i,1,k) + s11*q(i,2,k) - s14*dm(i,2,k)

                !           xt = min( xt, max(q(i,1), q(i,2)) )
                !           xt = max( xt, min(q(i,1), q(i,2)) )

                br(i,1,k) = xt - q(i,1,k)
                bl(i,2,k) = xt - q(i,2,k)
             enddo
             !         if ( jord<=9 ) then
             do j=0,2
                call pert_ppm(ilast-ifirst+1, q(ifirst,j,k), bl(ifirst,j,k), br(ifirst,j,k), 1)
             enddo
             !         endif
          enddo
       endif

       if( (je+1)==npy ) then
          do k = 1, km
             do i=ifirst,ilast
                bl(i,npy-2,k) = al(i,npy-2,k) - q(i,npy-2,k)
                !           xt = t11*( q(i,npy-1)+q(i,npy)) + t12*(q(i,npy-2)+q(i,npy+1))   &
                !                                         + t13*(dm(i,npy+1)-dm(i,npy-2))
                xt = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1,k)+q(i,npy,k))  &
                     -dya(i,npy-1)*(q(i,npy-2,k)+q(i,npy+1,k)))/(dya(i,npy-1)+dya(i,npy-2))
                br(i,npy-1,k) = xt - q(i,npy-1,k)
                bl(i,npy,k  ) = xt - q(i,npy,k)
                xt = s11*dq(i,npy,k) - s14*dm(i,npy+1,k) + q(i,npy,k)

                !           xt = min( xt, max( q(i,npy), q(i,npy+1)) )
                !           xt = max( xt, min( q(i,npy), q(i,npy+1)) )

                br(i,npy,k) = xt - q(i,npy,k)
                xt = s15*q(i,npy-1,k) + s11*q(i,npy-2,k) + s14*dm(i,npy-2,k)

                !           xt = min( xt, max( q(i,npy-2), q(i,npy-1)) )
                !           xt = max( xt, min( q(i,npy-2), q(i,npy-1)) )

                br(i,npy-2,k) = xt - q(i,npy-2,k)
                bl(i,npy-1,k) = xt - q(i,npy-1,k)
             enddo
             !         if ( jord<=9 ) then
             do j=npy-2,npy
                call pert_ppm(ilast-ifirst+1, q(ifirst,j,k), bl(ifirst,j,k), br(ifirst,j,k), 1)
             enddo
             !         endif
          enddo
       endif
!$ACC end kernels

    else
!          !---------------
!          ! grid_type == 4
!          !---------------
!       do k = 1, km
!          do j=jfirst-1,jlast+2
!             do i=ifirst,ilast
!                al(i,j,k) = 0.5*(q(i,j-1,k)+q(i,j,k)) + r3*(dm(i,j-1,k) - dm(i,j,k))
!             enddo
!          enddo
!       enddo
!       do k = 1, km
!          do j=jfirst-3,jlast+2
!             do i=ifirst,ilast
!                dq(i,j,k) = q(i,j+1,k) - q(i,j,k)
!             enddo
!          enddo
!       enddo
!       do k = 1, km
!          do j=jfirst-1,jlast+1
!             do i=ifirst,ilast
!                pmp = -2.*dq(i,j,k)
!                lac = pmp + 1.5*dq(i,j+1,k)
!                bl(i,j,k) = min(max(0.,pmp,lac), max(al(i,j,k)-q(i,j,k), min(0.,pmp,lac)))
!                pmp = 2.*dq(i,j-1,k)
!                lac = pmp - 1.5*dq(i,j-2,k)
!                br(i,j,k) = min(max(0.,pmp,lac), max(al(i,j+1,k)-q(i,j,k), min(0.,pmp,lac)))
!             enddo
!          enddo
!       enddo
    endif
!$ACC kernels loop collapse(3) async
    do k = 1, km
       do j=jfirst,jlast+1
          do i=ifirst,ilast
#ifdef INTEL_OPT
!             ct = c(i,j,k)
!             if( ct>0. ) then
!                jt = j-1
!                qe = br(i,j-1,k)
!             else
!                jt = j
!                qe = bl(i,j,k)
!             endif
!             ct = -abs(ct)
!             flux(i,j,k) = q(i,jt,k) + (1.+ct)*( qe + ct*(bl(i,jt)+br(i,jt)) )
#else
             if( c(i,j,k)>0. ) then
                flux(i,j,k) = q(i,j-1,k) + (1.-c(i,j,k))*(br(i,j-1,k)-c(i,j,k)*(bl(i,j-1,k)+br(i,j-1,k)))
             else
                flux(i,j,k) = q(i,j  ,k) + (1.+c(i,j,k))*(bl(i,j,k)+c(i,j,k)*(bl(i,j,k)+br(i,j,k)))
             endif
#endif
          enddo
       enddo
    enddo
!$ACC end kernels
!!$ACC wait
 else
       !-------------------------------
       ! For positive definite tracers:
       !-------------------------------
       ! jord=11: PPM mono constraint (Lin 2004)
       ! jord=12: Huynh 2nd constraint (Lin 2004) + positive definite (Lin & Rood 1996)
       ! jord>12: positive definite only (Lin & Rood 1996)
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
    do k = 1, km
       do j=js-2,je+2
          do i=ifirst,ilast
             xt = 0.25*(q(i,j+1,k) - q(i,j-1,k))
             dm(i,j,k) = sign(min(abs(xt), max(q(i,j-1,k), q(i,j,k), q(i,j+1,k)) - q(i,j,k),   &
                  q(i,j,k) - min(q(i,j-1,k), q(i,j,k), q(i,j+1,k))), xt)
          enddo
       enddo
    enddo
!$ACC end kernels
    if (grid_type < 3) then

       js3 = max(3,js-1); je3 = min(npy-3,je+1)
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
       do k = 1, km
          do j=js3,min(npy-2,je+2)
             do i=ifirst,ilast
                al(i,j,k) = 0.5*(q(i,j-1,k)+q(i,j,k)) + r3*(dm(i,j-1,k) - dm(i,j,k))
             enddo
          enddo
       enddo
!$ACC end kernels
       if ( jord==11 ) then
!          do k = 1, km
!             do j=js3,je3
!                do i=ifirst,ilast
!                   xt = 2.*dm(i,j,k)
!                   bl(i,j,k) = -sign(min(abs(xt), abs(al(i,j,k)-q(i,j,k))), xt)
!                   br(i,j,k) =  sign(min(abs(xt), abs(al(i,j+1,k)-q(i,j,k))), xt)
!                enddo
!             enddo
!          enddo
       elseif( jord==12 ) then
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=js-3,je+2
                do i=ifirst,ilast
                   dq(i,j,k) = q(i,j+1,k) - q(i,j,k)
                enddo
             enddo
          enddo
!$ACC end kernels
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=js3,je3
                do i=ifirst,ilast
                   pmp = -2.*dq(i,j,k)
                   lac = pmp + 1.5*dq(i,j+1,k)
                   bl(i,j,k) = min(max(0.,pmp,lac), max(al(i,j,k)-q(i,j,k), min(0.,pmp,lac)))
                   pmp = 2.*dq(i,j-1,k)
                   lac = pmp - 1.5*dq(i,j-2,k)
                   br(i,j,k) = min(max(0.,pmp,lac), max(al(i,j+1,k)-q(i,j,k), min(0.,pmp,lac)))
                enddo
             enddo
          enddo
!$ACC end kernels
       else  ! jord=13
!          do k = 1, km
!             do j=js-3,je+2
!                do i=ifirst,ilast
!                   dq(i,j,k) = q(i,j+1,k) - q(i,j,k)
!                enddo
!             enddo
!          enddo
!          do k = 1, km
!             do j=js3,je3
!                do i=ifirst,ilast
!                   if ( abs(dm(i,j-1,k))+abs(dm(i,j,k))+abs(dm(i,j+1,k)) < near_zero ) then
!                      bl(i,j,k) = 0.
!                      br(i,j,k) = 0.
!                   else
!                      bl(i,j,k) = al(i,j,k) - q(i,j,k)
!                      br(i,j,k) = al(i,j+1,k) - q(i,j,k)
!                      if( abs(3.*(bl(i,j,k)+br(i,j,k))) > abs(bl(i,j,k)-br(i,j,k)) ) then
!                         pmp_1 = -2.*dq(i,j,k)
!                         lac_1 = pmp_1 + 1.5*dq(i,j+1,k)
!                         bl(i,j,k) = min(max(0., pmp_1, lac_1), max(bl(i,j,k), min(0., pmp_1, lac_1)))
!                         pmp_2 = 2.*dq(i,j-1,k)
!                         lac_2 = pmp_2 - 1.5*dq(i,j-2,k)
!                         br(i,j,k) = min(max(0., pmp_2, lac_2), max(br(i,j,k), min(0., pmp_2, lac_2)))
!                      endif
!                   endif
!                enddo
!             enddo
!          enddo
       endif

       if ( jord/=11 ) then
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
             ! Positive definite constraint:
          do k = 1, km
             do j=js3,je3
                call pert_ppm(ilast-ifirst+1, q(ifirst,j,k), bl(ifirst,j,k), br(ifirst,j,k), 0)
             enddo
          enddo
!$ACC end kernels
       endif

!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          !--------------
          ! Fix the edges:
          !--------------
       if( js==1 ) then
           do k = 1, km
             do i=ifirst,ilast
                br(i,2,k) = al(i,3,k) - q(i,2,k)
                !           xt = t11*(q(i,0)+q(i,1)) + t12*(q(i,-1)+q(i,2))   &
                !              + t13*(dm(i,2)-dm(i,-1))
!!!         xt = 0.75*(q(i,0)+q(i,1)) - 0.25*(q(i,-1)+q(i,2))
                xt = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0,k)+q(i,1,k))  &
                     -dya(i,1)*(q(i,-1,k)+q(i,2,k))) / (dya(i,1)+dya(i,2))
                xt = max(0., xt)
                bl(i,1,k) = xt - q(i,1,k)
                br(i,0,k) = xt - q(i,0,k)
                xt = 4./7.*dm(i,-1,k) + 11./14.*q(i,-1,k) + 3./14.*q(i,0,k)
                xt = max(0., xt)
                bl(i,0,k) = xt - q(i,0,k)

                xt = 3./14.*q(i,1,k) + 11./14.*q(i,2,k) - 4./7.*dm(i,2,k)
                xt = max(0., xt)
                br(i,1,k) = xt - q(i,1,k)
                bl(i,2,k) = xt - q(i,2,k)
             enddo
             do j=0,2
                call pert_ppm(ilast-ifirst+1, q(ifirst,j,k), bl(ifirst,j,k), br(ifirst,j,k), 1)
             enddo
          enddo
       endif

       if( (je+1)==npy ) then
          do k = 1, km
             do i=ifirst,ilast
                bl(i,npy-2,k) = al(i,npy-2,k) - q(i,npy-2,k)
                !           xt = t11*(q(i,npy-1)+q(i,npy)) + t12*(q(i,npy-2)+q(i,npy+1))   &
                !               + t13*(dm(i,npy+1)-dm(i,npy-2))
!!!         xt = 0.75*(q(i,npy-1)+q(i,npy)) - 0.25*(q(i,npy-2)+q(i,npy+1))
                xt = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1,k)+q(i,npy,k)) &
                     - dya(i,npy-1)*(q(i,npy-2,k)+q(i,npy+1,k)))  &
                     / ( dya(i,npy-1)+dya(i,npy-2) )
                xt = max(0., xt)
                br(i,npy-1,k) = xt - q(i,npy-1,k)
                bl(i,npy,k  ) = xt - q(i,npy,k)
                xt = 3./14.*q(i,npy,k) + 11./14.*q(i,npy+1,k) - 4./7.*dm(i,npy+1,k)
                xt = max(0., xt)
                br(i,npy,k) = xt - q(i,npy,k)
                xt = 3./14.*q(i,npy-1,k) + 11./14.*q(i,npy-2,k) + 4./7.*dm(i,npy-2,k)
                xt = max(0., xt)
                br(i,npy-2,k) = xt - q(i,npy-2,k)
                bl(i,npy-1,k) = xt - q(i,npy-1,k)
             enddo
             do j=npy-2,npy
                call pert_ppm(ilast-ifirst+1, q(ifirst,j,k), bl(ifirst,j,k), br(ifirst,j,k), 1)
             enddo
          enddo
       endif
!$ACC end kernels

    else
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
       do k = 1, km
          do j=js-1,je+2
             do i=ifirst,ilast
                al(i,j,k) = 0.5*(q(i,j-1,k)+q(i,j,k)) + r3*(dm(i,j-1,k) - dm(i,j,k))
             enddo
          enddo
       enddo
!$ACC end kernels

       if ( jord==11 ) then
!          do k = 1, km
!             do j=js-1,je+1
!                do i=ifirst,ilast
!                   xt = 2.*dm(i,j,k)
!                   bl(i,j,k) = -sign(min(abs(xt), abs(al(i,j,k)-q(i,j,k))), xt)
!                   br(i,j,k) =  sign(min(abs(xt), abs(al(i,j+1,k)-q(i,j,k))), xt)
!                enddo
!             enddo
!          enddo
       elseif( jord==12 ) then
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=js-3,je+2
                do i=ifirst,ilast
                   dq(i,j,k) = q(i,j+1,k) - q(i,j,k)
                enddo
             enddo
          enddo
!$ACC end kernels
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=js-1,je+1
                do i=ifirst,ilast
                   pmp = -2.*dq(i,j,k)
                   lac = pmp + 1.5*dq(i,j+1,k)
                   bl(i,j,k) = min(max(0.,pmp,lac), max(al(i,j,k)-q(i,j,k), min(0.,pmp,lac)))
                   pmp = 2.*dq(i,j-1,k)
                   lac = pmp - 1.5*dq(i,j-2,k)
                   br(i,j,k) = min(max(0.,pmp,lac), max(al(i,j+1,k)-q(i,j,k), min(0.,pmp,lac)))
                enddo
             enddo
          enddo
!$ACC end kernels
       else
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=js-1,je+1
                do i=ifirst,ilast
                   bl(i,j,k) = al(i,j,k) - q(i,j,k)
                   br(i,j,k) = al(i,j+1,k) - q(i,j,k)
                enddo
             enddo
          enddo
!$ACC end kernels
       endif

       if ( jord/=11 ) then
             ! Positive definite constraint:
!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
          do k = 1, km
             do j=js-1,je+1
                call pert_ppm(ilast-ifirst+1, q(ifirst,j,k), bl(ifirst,j,k), br(ifirst,j,k), 0)
             enddo
          enddo
!$ACC end kernels
       endif

    endif

!$ACC kernels present(flux,dm,q,c,dya,al,bl,br,dq) async
    do k = 1, km
       do j=js,je+1
          do i=ifirst,ilast
             if( c(i,j,k)>0. ) then
                flux(i,j,k) = q(i,j-1,k) + (1.-c(i,j,k))*(br(i,j-1,k)-c(i,j,k)*(bl(i,j-1,k)+br(i,j-1,k)))
             else
                flux(i,j,k) = q(i,j,k  ) + (1.+c(i,j,k))*(bl(i,j,k)+c(i,j,k)*(bl(i,j,k)+br(i,j,k)))
             endif
          enddo
       enddo
    enddo
!$ACC end kernels
!!$ACC wait
 end if

 end subroutine fyppm


 subroutine pert_ppm(im2, a0, al, ar, iv)
!$acc routine vector
 integer, intent(in):: im2
 integer, intent(in):: iv
 real, intent(in)   :: a0(im2)
 real, intent(inout):: al(im2), ar(im2)
! Local:
 real a4, da1, da2, a6da, fmin
 integer i
 real, parameter:: r12 = 1./12.
 integer, parameter :: im = 56

!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
        a4 = -3.*(ar(i) + al(i))
       da1 =      ar(i) - al(i)
      if( abs(da1) < -a4 ) then
         fmin = a0(i) + 0.25/a4*da1**2 + a4*r12
         if( fmin < 0. ) then
             if( ar(i)>0. .and. al(i)>0. ) then
                 ar(i) = 0.
                 al(i) = 0.
             elseif( da1 > 0. ) then
                 ar(i) = -2.*al(i)
             else
                 al(i) = -2.*ar(i)
             endif
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
       if ( al(i)*ar(i) < 0. ) then
            da1 = al(i) - ar(i)
            da2 = da1**2
            a6da = 3.*(al(i)+ar(i))*da1
            if( a6da < -da2 ) then
                ar(i) = -2.*al(i)
            elseif( a6da > da2 ) then
                al(i) = -2.*ar(i)
            endif
       else
! effect of dm=0 included here
            al(i) = 0.
            ar(i) = 0.
       endif
  enddo
 endif

 end subroutine pert_ppm

end module fyppm_mod

