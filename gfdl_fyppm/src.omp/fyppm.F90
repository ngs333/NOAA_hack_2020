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

 integer :: ifirst = 1, ilast = 96
!$acc declare copyin(ifirst,ilast,grid_type)

 contains

 subroutine fyppm(c,  q,  flux, jord, jfirst, jlast, npy, dm, ppm_limiter, &
                  dya, isd, ied, jsd, jed, js, je, ng )
!$acc routine seq
 integer, INTENT(IN) :: jfirst, jlast               !  Y-Dir strip
 integer, INTENT(IN) :: isd, ied
 integer, INTENT(IN) :: jsd, jed
 integer, INTENT(IN) :: js, je
 integer, INTENT(IN) :: ng
 integer, INTENT(IN) :: jord
 integer, INTENT(IN) :: npy
 real   , intent(IN) :: ppm_limiter
 real   , INTENT(IN) :: q(ifirst:ilast,jfirst-ng:jlast+ng)
 real   , intent(in) :: c(isd:ied,js:je+1 )                 ! Courant number
 real   , INTENT(OUT):: flux(ifirst:ilast,jfirst:jlast+1)   ! Flux
 real   , INTENT(OUT):: dm(ifirst:ilast,jfirst-2:jlast+2)
 real   , intent(in) :: dya(isd:ied,js:je ) 

! Local:
 logical extm(ifirst:ilast,jfirst-2:jlast+2)
 real al(ifirst:ilast,jfirst-1:jlast+2)
 real bl(ifirst:ilast,jfirst-1:jlast+1)
 real br(ifirst:ilast,jfirst-1:jlast+1)
 real dq(ifirst:ilast,jfirst-3:jlast+2)
 real dl, dr, pmp, lac, ct, qe
 real pmp_1, lac_1, pmp_2, lac_2
 real xt, x0, x1
 integer i, j, js3, je3, jt

 if (jord<=4) then
   do j=js-2,je+2
      do i=ifirst,ilast
         xt = 0.25*(q(i,j+1) - q(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                            q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
      enddo
   enddo

  if (grid_type < 3) then
   do j=max(3,js-1),min(npy-2,je+2)
      do i=ifirst,ilast
         al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo
!--------------
! Fix the edges:
!--------------
      if( js==1 ) then
         do i=ifirst,ilast
            x0 = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0)+q(i,1))   &
               -dya(i,1)*(q(i,-1)+q(i,2))) / ( dya(i,1)+dya(i,2) )
            al(i,1) = x0
            x1 = s15*q(i,0) + s11*q(i,-1) + s14*dm(i,-1)
            dm(i,0) = 0.5*(x0 - x1)
            dm(i,0) = sign(min(abs(dm(i,0)), max(q(i,0), x0, x1) - q(i,0),   &
                          q(i,0) - min(q(i,0), x0, x1)), dm(i,0))
            al(i,0) = 0.5*(q(i,-1)+q(i,0)) + r3*(dm(i,-1) - dm(i,0))
!
                 x1 = s15*q(i,1) + s11*q(i,2) - s14*dm(i,2)
            dm(i,1) = 0.5*(x1 - x0)
            dm(i,1) = sign(min(abs(dm(i,1)), max(q(i,1), x0, x1) - q(i,1),    &
                                    q(i,1) - min(q(i,1), x0, x1)), dm(i,1))
            al(i,2) = 0.5*(q(i,1)+q(i,2)) + r3*(dm(i,1) - dm(i,2))
         enddo
      endif

      if( (je+1)==npy ) then
         do i=ifirst,ilast
            x0 = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1)+q(i,npy))  &
               -dya(i,npy-1)*(q(i,npy-2)+q(i,npy+1)))/(dya(i,npy-1)+dya(i,npy-2))
            al(i,npy) = x0
            x1 = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*dm(i,npy-2)
            dm(i,npy-1) = 0.5*(x0 - x1)
            dm(i,npy-1) = sign(min(abs(dm(i,npy-1)), max(q(i,npy-1), x0, x1) - q(i,npy-1),  &
                                        q(i,npy-1) - min(q(i,npy-1), x0, x1)), dm(i,npy-1))
            al(i,npy-1) = 0.5*(q(i,npy-2)+q(i,npy-1)) + r3*(dm(i,npy-2) - dm(i,npy-1))
!
            x1 = s15*q(i,npy) + s11*q(i,npy+1) - s14*dm(i,npy+1)
            dm(i,npy) = 0.5*(x1 - x0)
            dm(i,npy) = sign(min(abs(dm(i,npy)), max(q(i,npy), x0, x1) - q(i,npy),   &
                                      q(i,npy) - min(q(i,npy), x0, x1)), dm(i,npy))
            al(i,npy+1) = 0.5*(q(i,npy)+q(i,npy+1)) + r3*(dm(i,npy) - dm(i,npy+1))
         enddo
      endif
  else
! Doubly periodic BC:
      do j=js-1,je+2
         do i=ifirst,ilast
            al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
         enddo
      enddo
  endif

  if ( jord==3 ) then
      do j=js-1,je+1
         do i=ifirst,ilast
            bl(i,j) = al(i,j  ) - q(i,j)
            br(i,j) = al(i,j+1) - q(i,j)
         enddo
         call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
      enddo
      do j=js,je+1
         do i=ifirst,ilast
         if( c(i,j)>0. ) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
         enddo
      enddo
  else
! Inlined limiter
!!!!$acc loop
   do j=js,je+1
      do i=ifirst,ilast
         if( c(i,j)>0. ) then
             xt = ppm_limiter*dm(i,j-1)
             dl = sign(min(abs(xt), abs(al(i,j-1)-q(i,j-1))), xt)
             dr = sign(min(abs(xt), abs(al(i,j)-q(i,j-1))),   xt)
             flux(i,j) = q(i,j-1) + (1.-c(i,j))*(c(i,j)*(dl-dr)+dr)
         else
             xt = ppm_limiter*dm(i,j)
             dl = sign(min(abs(xt), abs(al(i,j)-q(i,j))),   xt)
             dr = sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
             flux(i,j) = q(i,j) - (1.+c(i,j))*(c(i,j)*(dl-dr)+dl)
         endif
      enddo
   enddo
  endif

 elseif (jord==5) then
! PPM with Hunyh's 2nd constraint

   do j=jfirst-3, jlast+2
      do i=ifirst,ilast
         dq(i,j) = q(i,j+1) - q(i,j)
      enddo
   enddo

   do j=jfirst-2,jlast+2
      do i=ifirst,ilast
         xt = 0.25*(q(i,j+1) - q(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                            q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
      enddo
   enddo

   do j=jfirst-1,jlast+2
      do i=ifirst,ilast
         al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      enddo
   enddo

   do j=jfirst-1,jlast+1
      do i=ifirst,ilast
            pmp = -2.*dq(i,j) 
            lac = pmp + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0.,pmp,lac), max(al(i,j)-q(i,j), min(0.,pmp,lac)))
            pmp = 2.*dq(i,j-1)
            lac = pmp - 1.5*dq(i,j-2)
            br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
      enddo
   enddo

   do j=jfirst,jlast+1
      do i=ifirst,ilast
         if(c(i,j)>0.) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo

 elseif( jord==6 .or. jord==7 ) then

   do j=js-3,je+2
      do i=ifirst,ilast
         dq(i,j) = q(i,j+1) - q(i,j)
      enddo
   enddo

   do j=js-2,je+2
      do i=ifirst,ilast
         if ( dq(i,j-1)*dq(i,j) > 0. ) then
              extm(i,j) = .false.
         else
              extm(i,j) = .true.
         endif
      enddo
   enddo

   do j=max(3,js-1),min(npy-3,je+1)
      do i=ifirst,ilast
         if ( extm(i,j) .and. (extm(i,j-1) .or. extm(i,j+1)) ) then
! 2-delta-wave filter
              bl(i,j) = 0.
              br(i,j) = 0.
         else
              bl(i,j) = b5*q(i,j-2) + b4*q(i,j-1) + b3*q(i,j) + b2*q(i,j+1) + b1*q(i,j+2)
              br(i,j) = b1*q(i,j-2) + b2*q(i,j-1) + b3*q(i,j) + b4*q(i,j+1) + b5*q(i,j+2)
         endif
      enddo
   enddo

   if( js==1 ) then
         do i=ifirst,ilast
!           br(i,2) = al(i,3) - q(i,2)
            br(i,2) = p1*(q(i,2)+q(i,3)) + p2*(q(i,1)+q(i,4)) - q(i,2)
            xt = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0)+q(i,1))   &
               -dya(i,1)*(q(i,-1)+q(i,2))) / ( dya(i,1)+dya(i,2) )
            bl(i,1) = xt - q(i,1)
            br(i,0) = xt - q(i,0)

!           xt = s14*0.25*(q(i,0)-q(i,-2)) - s11*(q(i,0)-q(i,-1)) + q(i,0)
            xt = c1*q(i,-2) + c2*q(i,-1) + c3*q(i,0)
            xt = min( xt, max(q(i,-1), q(i,0)) )
            xt = max( xt, min(q(i,-1), q(i,0)) )
            bl(i,0) = xt - q(i,0)

!           xt = s15*q(i,1) + s11*q(i,2) - s14*0.25*(q(i,3)-q(i,1))
            xt = c3*q(i,1) + c2*q(i,2) + c1*q(i,3)
            xt = min( xt, max(q(i,1), q(i,2)) )
            xt = max( xt, min(q(i,1), q(i,2)) )
            br(i,1) = xt - q(i,1)
            bl(i,2) = xt - q(i,2)
         enddo
   endif

   if( (je+1)==npy ) then
         do i=ifirst,ilast
!           bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
            bl(i,npy-2) = p1*(q(i,npy-3)+q(i,npy-2)) + p2*(q(i,npy-4)+q(i,npy-1)) - q(i,npy-2)
            xt = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1)+q(i,npy))  &
               -dya(i,npy-1)*(q(i,npy-2)+q(i,npy+1)))/(dya(i,npy-1)+dya(i,npy-2))
            br(i,npy-1) = xt - q(i,npy-1)
            bl(i,npy  ) = xt - q(i,npy)

!           xt = s11*(q(i,npy+1)-q(i,npy)) - s14*0.25*(q(i,npy+2)-q(i,npy)) + q(i,npy)
            xt = c3*q(i,npy) + c2*q(i,npy+1) + c1*q(i,npy+2)
            xt = min( xt, max(q(i,npy), q(i,npy+1)) )
            xt = max( xt, min(q(i,npy), q(i,npy+1)) )
            br(i,npy) = xt - q(i,npy)

!           xt = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*0.25*(q(i,npy-1)-q(i,npy-3))
            xt = c1*q(i,npy-3) + c2*q(i,npy-2) + c3*q(i,npy-1)
            xt = min( xt, max(q(i,npy-2), q(i,npy-1)) )
            xt = max( xt, min(q(i,npy-2), q(i,npy-1)) )
            br(i,npy-2) = xt - q(i,npy-2)
            bl(i,npy-1) = xt - q(i,npy-1)
         enddo
   endif

! Positive definite constraint:
   if ( jord==7 ) then
        do j=jfirst-1,jlast+1
           call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 0)
        enddo
   endif

   do j=jfirst,jlast+1
      do i=ifirst,ilast
         if(c(i,j)>0.) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo

 elseif( jord<=10 ) then    ! jord=8, 9, 10

   do j=js-2,je+2
      do i=ifirst,ilast
              xt = 0.25*(q(i,j+1) - q(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                            q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
      enddo
   enddo

   if (grid_type < 3) then

       do j=max(3,js-1),min(npy-2,je+2)
          do i=ifirst,ilast
             al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
          enddo
       enddo

       do j=js-3,je+2
          do i=ifirst,ilast
             dq(i,j) = q(i,j+1) - q(i,j)
          enddo
       enddo
      
       if ( jord==8 ) then
         do j=max(3,js-1),min(npy-3,je+1)
         do i=ifirst,ilast
            xt = 2.*dm(i,j)
            bl(i,j) = -sign(min(abs(xt), abs(al(i,j)-q(i,j))),   xt)
            br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
         enddo
         enddo
       elseif( jord==9 ) then
         do j=max(3,js-1),min(npy-3,je+1)
         do i=ifirst,ilast
              pmp_1 = -2.*dq(i,j) 
              lac_1 = pmp_1 + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0., pmp_1, lac_1), max(al(i,j  )-q(i,j), min(0., pmp_1, lac_1)))
              pmp_2 = 2.*dq(i,j-1)
              lac_2 = pmp_2 - 1.5*dq(i,j-2)
            br(i,j) = min(max(0., pmp_2, lac_2), max(al(i,j+1)-q(i,j), min(0., pmp_2, lac_2)))
         enddo
         enddo
       else    ! jord=10
         do j=max(3,js-1),min(npy-3,je+1)
            do i=ifirst,ilast
               bl(i,j) = al(i,j  ) - q(i,j)
               br(i,j) = al(i,j+1) - q(i,j)
              if ( abs(dm(i,j-1))+abs(dm(i,j))+abs(dm(i,j+1)) < near_zero ) then
                   bl(i,j) = 0.
                   br(i,j) = 0.
              elseif( abs(3.*(bl(i,j)+br(i,j))) > abs(bl(i,j)-br(i,j)) ) then
                     pmp_1 = -2.*dq(i,j) 
                     lac_1 = pmp_1 + 1.5*dq(i,j+1)
                   bl(i,j) = min(max(0.,pmp_1,lac_1), max(bl(i,j), min(0.,pmp_1,lac_1)))
                     pmp_2 = 2.*dq(i,j-1)
                     lac_2 = pmp_2 - 1.5*dq(i,j-2)
                   br(i,j) = min(max(0.,pmp_2,lac_2), max(br(i,j), min(0.,pmp_2,lac_2)))
            endif
         enddo
         enddo
       endif

!--------------
! Fix the edges:
!--------------
      if( js==1 ) then
         do i=ifirst,ilast
            br(i,2) = al(i,3) - q(i,2)
!           xt = t11*(q(i,0)+q(i,1)) + t12*(q(i,-1)+q(i,2))   &
!                                  + t13*(dm(i,2)-dm(i,-1))
            xt = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0)+q(i,1))   &
               -dya(i,1)*(q(i,-1)+q(i,2))) / ( dya(i,1)+dya(i,2) )
            bl(i,1) = xt - q(i,1)
            br(i,0) = xt - q(i,0)
            xt = s14*dm(i,-1) - s11*dq(i,-1) + q(i,0)

!           xt = min( xt, max(q(i,-1), q(i,0)) )
!           xt = max( xt, min(q(i,-1), q(i,0)) )

            bl(i,0) = xt - q(i,0)
            xt = s15*q(i,1) + s11*q(i,2) - s14*dm(i,2)

!           xt = min( xt, max(q(i,1), q(i,2)) )
!           xt = max( xt, min(q(i,1), q(i,2)) )

            br(i,1) = xt - q(i,1)
            bl(i,2) = xt - q(i,2)
         enddo
!         if ( jord<=9 ) then
            do j=0,2
               call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
            enddo
!         endif
      endif

      if( (je+1)==npy ) then
         do i=ifirst,ilast
            bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
!           xt = t11*( q(i,npy-1)+q(i,npy)) + t12*(q(i,npy-2)+q(i,npy+1))   &
!                                         + t13*(dm(i,npy+1)-dm(i,npy-2))
            xt = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1)+q(i,npy))  &
               -dya(i,npy-1)*(q(i,npy-2)+q(i,npy+1)))/(dya(i,npy-1)+dya(i,npy-2))
            br(i,npy-1) = xt - q(i,npy-1)
            bl(i,npy  ) = xt - q(i,npy)
            xt = s11*dq(i,npy) - s14*dm(i,npy+1) + q(i,npy)

!           xt = min( xt, max( q(i,npy), q(i,npy+1)) )
!           xt = max( xt, min( q(i,npy), q(i,npy+1)) )

            br(i,npy) = xt - q(i,npy)
            xt = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*dm(i,npy-2)

!           xt = min( xt, max( q(i,npy-2), q(i,npy-1)) )
!           xt = max( xt, min( q(i,npy-2), q(i,npy-1)) )

            br(i,npy-2) = xt - q(i,npy-2)
            bl(i,npy-1) = xt - q(i,npy-1)
         enddo
!         if ( jord<=9 ) then
            do j=npy-2,npy
               call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
            enddo
!         endif
      endif

   else
!---------------
! grid_type == 4
!---------------

      do j=jfirst-1,jlast+2
         do i=ifirst,ilast
            al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
         enddo
      enddo

      do j=jfirst-3,jlast+2
         do i=ifirst,ilast
            dq(i,j) = q(i,j+1) - q(i,j)
         enddo
      enddo
      
      do j=jfirst-1,jlast+1
         do i=ifirst,ilast
            pmp = -2.*dq(i,j) 
            lac = pmp + 1.5*dq(i,j+1)
            bl(i,j) = min(max(0.,pmp,lac), max(al(i,j  )-q(i,j), min(0.,pmp,lac)))
            pmp = 2.*dq(i,j-1)
            lac = pmp - 1.5*dq(i,j-2)
            br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
         enddo
      enddo

   endif

   do j=jfirst,jlast+1
      do i=ifirst,ilast
#ifdef INTEL_OPT
         ct = c(i,j)
         if( ct>0. ) then
             jt = j-1
             qe = br(i,j-1) 
         else
             jt = j
             qe = bl(i,j) 
         endif
         ct = -abs(ct)
         flux(i,j) = q(i,jt) + (1.+ct)*( qe + ct*(bl(i,jt)+br(i,jt)) )
#else
         if( c(i,j)>0. ) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
#endif
      enddo
   enddo

 else
!-------------------------------
! For positive definite tracers:
!-------------------------------
! jord=11: PPM mono constraint (Lin 2004)
! jord=12: Huynh 2nd constraint (Lin 2004) + positive definite (Lin & Rood 1996)
! jord>12: positive definite only (Lin & Rood 1996)


   do j=js-2,je+2
      do i=ifirst,ilast
         xt = 0.25*(q(i,j+1) - q(i,j-1))
         dm(i,j) = sign(min(abs(xt), max(q(i,j-1), q(i,j), q(i,j+1)) - q(i,j),   &
                            q(i,j) - min(q(i,j-1), q(i,j), q(i,j+1))), xt)
      enddo
   enddo

   if (grid_type < 3) then

      js3 = max(3,js-1); je3 = min(npy-3,je+1)

      do j=js3,min(npy-2,je+2)
         do i=ifirst,ilast
            al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
         enddo
      enddo

      if ( jord==11 ) then
         do j=js3,je3
            do i=ifirst,ilast
               xt = 2.*dm(i,j)
               bl(i,j) = -sign(min(abs(xt), abs(al(i,j  )-q(i,j))), xt)
               br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
            enddo
         enddo
      elseif( jord==12 ) then
         do j=js-3,je+2
            do i=ifirst,ilast
               dq(i,j) = q(i,j+1) - q(i,j)
            enddo
         enddo
         do j=js3,je3
            do i=ifirst,ilast
               pmp = -2.*dq(i,j) 
               lac = pmp + 1.5*dq(i,j+1)
               bl(i,j) = min(max(0.,pmp,lac), max(al(i,j  )-q(i,j), min(0.,pmp,lac)))
               pmp = 2.*dq(i,j-1)
               lac = pmp - 1.5*dq(i,j-2)
               br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
            enddo
         enddo
      else  ! jord=13

         do j=js-3,je+2
            do i=ifirst,ilast
               dq(i,j) = q(i,j+1) - q(i,j)
            enddo
         enddo
         do j=js3,je3
            do i=ifirst,ilast
               if ( abs(dm(i,j-1))+abs(dm(i,j))+abs(dm(i,j+1)) < near_zero ) then
                    bl(i,j) = 0.
                    br(i,j) = 0.
               else
                 bl(i,j) = al(i,j  ) - q(i,j)
                 br(i,j) = al(i,j+1) - q(i,j)
                 if( abs(3.*(bl(i,j)+br(i,j))) > abs(bl(i,j)-br(i,j)) ) then
                      pmp_1 = -2.*dq(i,j) 
                      lac_1 = pmp_1 + 1.5*dq(i,j+1)
                    bl(i,j) = min(max(0., pmp_1, lac_1), max(bl(i,j), min(0., pmp_1, lac_1)))
                      pmp_2 = 2.*dq(i,j-1)
                      lac_2 = pmp_2 - 1.5*dq(i,j-2)
                    br(i,j) = min(max(0., pmp_2, lac_2), max(br(i,j), min(0., pmp_2, lac_2)))
                 endif
               endif
            enddo
         enddo
      endif
      
      if ( jord/=11 ) then
! Positive definite constraint:
         do j=js3,je3
            call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 0)
         enddo
      endif

!--------------
! Fix the edges:
!--------------
      if( js==1 ) then
         do i=ifirst,ilast
            br(i,2) = al(i,3) - q(i,2)
!           xt = t11*(q(i,0)+q(i,1)) + t12*(q(i,-1)+q(i,2))   &
!              + t13*(dm(i,2)-dm(i,-1))
!!!         xt = 0.75*(q(i,0)+q(i,1)) - 0.25*(q(i,-1)+q(i,2))
            xt = 0.5*((2.*dya(i,1)+dya(i,2))*(q(i,0)+q(i,1))  &
               -dya(i,1)*(q(i,-1)+q(i,2))) / (dya(i,1)+dya(i,2))
            xt = max(0., xt)
            bl(i,1) = xt - q(i,1)
            br(i,0) = xt - q(i,0)
            xt = 4./7.*dm(i,-1) + 11./14.*q(i,-1) + 3./14.*q(i,0)
            xt = max(0., xt)
            bl(i,0) = xt - q(i,0)

            xt = 3./14.*q(i,1) + 11./14.*q(i,2) - 4./7.*dm(i,2)
            xt = max(0., xt)
            br(i,1) = xt - q(i,1)
            bl(i,2) = xt - q(i,2)
         enddo
         do j=0,2
            call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
         enddo
      endif

      if( (je+1)==npy ) then
         do i=ifirst,ilast
            bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
!           xt = t11*(q(i,npy-1)+q(i,npy)) + t12*(q(i,npy-2)+q(i,npy+1))   &
!               + t13*(dm(i,npy+1)-dm(i,npy-2))
!!!         xt = 0.75*(q(i,npy-1)+q(i,npy)) - 0.25*(q(i,npy-2)+q(i,npy+1))
            xt = 0.5*((2.*dya(i,npy-1)+dya(i,npy-2))*(q(i,npy-1)+q(i,npy)) &
               - dya(i,npy-1)*(q(i,npy-2)+q(i,npy+1)))  &
                / ( dya(i,npy-1)+dya(i,npy-2) )
            xt = max(0., xt)
            br(i,npy-1) = xt - q(i,npy-1)
            bl(i,npy  ) = xt - q(i,npy)
            xt = 3./14.*q(i,npy) + 11./14.*q(i,npy+1) - 4./7.*dm(i,npy+1)
            xt = max(0., xt)
            br(i,npy) = xt - q(i,npy)
            xt = 3./14.*q(i,npy-1) + 11./14.*q(i,npy-2) + 4./7.*dm(i,npy-2)
            xt = max(0., xt)
            br(i,npy-2) = xt - q(i,npy-2)
            bl(i,npy-1) = xt - q(i,npy-1)
         enddo
         do j=npy-2,npy
            call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 1)
         enddo
      endif

   else

      do j=js-1,je+2
         do i=ifirst,ilast
            al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
         enddo
      enddo

      if ( jord==11 ) then
         do j=js-1,je+1
            do i=ifirst,ilast
               xt = 2.*dm(i,j)
               bl(i,j) = -sign(min(abs(xt), abs(al(i,j  )-q(i,j))), xt)
               br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j))), xt)
            enddo
         enddo
      elseif( jord==12 ) then
         do j=js-3,je+2
            do i=ifirst,ilast
               dq(i,j) = q(i,j+1) - q(i,j)
            enddo
         enddo
         do j=js-1,je+1
            do i=ifirst,ilast
               pmp = -2.*dq(i,j) 
               lac = pmp + 1.5*dq(i,j+1)
               bl(i,j) = min(max(0.,pmp,lac), max(al(i,j  )-q(i,j), min(0.,pmp,lac)))
               pmp = 2.*dq(i,j-1)
               lac = pmp - 1.5*dq(i,j-2)
               br(i,j) = min(max(0.,pmp,lac), max(al(i,j+1)-q(i,j), min(0.,pmp,lac)))
            enddo
         enddo
      else
         do j=js-1,je+1
            do i=ifirst,ilast
               bl(i,j) = al(i,j  ) - q(i,j)
               br(i,j) = al(i,j+1) - q(i,j)
            enddo
         enddo
      endif

      if ( jord/=11 ) then
! Positive definite constraint:
         do j=js-1,je+1
            call pert_ppm(ilast-ifirst+1, q(ifirst,j), bl(ifirst,j), br(ifirst,j), 0)
         enddo
      endif

   endif

   do j=js,je+1
      do i=ifirst,ilast
         if( c(i,j)>0. ) then
            flux(i,j) = q(i,j-1) + (1.-c(i,j))*(br(i,j-1)-c(i,j)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j) = q(i,j  ) + (1.+c(i,j))*(bl(i,j  )+c(i,j)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo
 endif

 end subroutine fyppm


 subroutine pert_ppm(im2, a0, al, ar, iv)
!$acc routine
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

