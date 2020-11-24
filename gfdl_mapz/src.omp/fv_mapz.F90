!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
! SJL: Apr 12, 2012
! This revision may actually produce rounding level differences due to the elimination of KS to compute
! pressure level for remapping.
module fv_mapz_mod


  implicit none

  real, parameter:: ptop_min=1.d-8
  real, parameter:: RADIUS = 6.371e+6
  real, parameter:: PI = 3.14159265358979323846
  real, parameter:: RVGAS  = 461.50
  real, parameter:: RDGAS  = 287.04
  real, parameter:: GRAV   = 9.8
  real, parameter:: HLV = 2.5e6
  real, parameter:: HLF = 3.34e5
  real, parameter:: KAPPA  = 2.0/7.0
  real, parameter:: CP_AIR = RDGAS/KAPPA
  real, parameter:: CP_VAPOR = 4.0*RVGAS
  real, parameter:: consv_min= 0.001   ! below which no correction applies
  real, parameter:: t_min= 184.   ! below which applies stricter constraint
  real, parameter:: r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real, parameter:: cv_vap = 3.*rvgas  ! 1384.5
  real, parameter:: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68
! real, parameter:: c_ice = 2106.           ! heat capacity of ice at 0.C
  real, parameter:: c_ice = 1972.           ! heat capacity of ice at -15.C
  real, parameter:: c_liq = 4.1855e+3    ! GFS: heat capacity of water at 0C
! real, parameter:: c_liq = 4218.        ! ECMWF-IFS
  real, parameter:: cp_vap = cp_vapor   ! 1846.
  real, parameter:: tice = 273.16

  real(kind=4) :: E_Flux = 0.
  private

  type fv_grid_type
     real, allocatable ::  rsin2(:,:)
     real, allocatable ::  cosa_s(:,:)
     real, allocatable ::  area(:,:)

  end type fv_grid_type

  public Lagrangian_to_Eulerian, fv_grid_type
  public CP_AIR, RVGAS, RDGAS


contains

 subroutine Lagrangian_to_Eulerian(last_step, consv, ps, pe, delp, pkz, pk,   &
                                   pdt, km, is,ie,js,je, isd,ied,jsd,jed,       &
                      nq, nwat, sphum, u, v, w, delz, pt, q, hs, r_vir, cp,  &
                      akap, kord_mt, kord_wz, kord_tr, kord_tm,  peln, te0_2d,        &
                      ng, ua, va, omga, te, ws,   &
                      ptop, ak, bk, gridstruct, &
                      hydrostatic, do_omega)
  logical, intent(in):: last_step
  real,    intent(in):: pdt                   ! phys time step
  integer, intent(in):: km
  integer, intent(in):: nq                    ! number of tracers (including h2o)
  integer, intent(in):: nwat
  integer, intent(in):: sphum                 ! index for water vapor (specific humidity)
  integer, intent(in):: ng
  integer, intent(in):: is,ie,isd,ied         ! starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         ! starting & ending Y-Dir index
  integer, intent(in):: kord_mt               ! Mapping order for the vector winds
  integer, intent(in):: kord_wz               ! Mapping order/option for w
  integer, intent(in):: kord_tr(nq)           ! Mapping order for tracers
  integer, intent(in):: kord_tm               ! Mapping order for thermodynamics

  real, intent(in):: consv                 ! factor for TE conservation
  real, intent(in):: r_vir
  real, intent(in):: cp
  real, intent(in):: akap
  real, intent(in):: hs(isd:ied,jsd:jed)  ! surface geopotential
  real, intent(inout):: te0_2d(is:ie,js:je)
  real, intent(in):: ws(is:ie,js:je)

  logical, intent(in):: do_omega
  real, intent(in) :: ptop
  real, intent(in) :: ak(km+1)
  real, intent(in) :: bk(km+1)
  type(fv_grid_type), intent(IN), target :: gridstruct

! !INPUT/OUTPUT
  real, intent(inout):: pk(is:ie,js:je,km+1) ! pe to the kappa
  real, intent(inout):: q(isd:ied,jsd:jed,km,*)
  real, intent(inout):: delp(isd:ied,jsd:jed,km) ! pressure thickness
  real, intent(inout)::  pe(is-1:ie+1,km+1,js-1:je+1) ! pressure at layer edges
  real, intent(inout):: ps(isd:ied,jsd:jed)      ! surface pressure

! u-wind will be ghosted one latitude to the north upon exit
  real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   ! u-wind (m/s)
  real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   ! v-wind (m/s)
  real, intent(inout)::  w(isd:     ,jsd:     ,1:)   ! vertical velocity (m/s)
  real, intent(inout):: pt(isd:ied  ,jsd:jed  ,km)   ! cp*virtual potential temperature
                                                     ! as input; output: temperature
  real, intent(inout), dimension(isd:,jsd:,1:)::delz
  logical, intent(in):: hydrostatic

  real, intent(inout)::   ua(isd:ied,jsd:jed,km)   ! u-wind (m/s) on physics grid
  real, intent(inout)::   va(isd:ied,jsd:jed,km)   ! v-wind (m/s) on physics grid
  real, intent(inout):: omga(isd:ied,jsd:jed,km)   ! vertical press. velocity (pascal/sec)
  real, intent(inout)::   peln(is:ie,km+1,js:je)     ! log(pe)
  real, intent(out)::    pkz(is:ie,js:je,km)       ! layer-mean pk for converting t to pt
  real, intent(out)::     te(isd:ied,jsd:jed,km)

! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
  real, dimension(is:ie,js:je):: te_2d, zsum0, zsum1, dpln
  real, dimension(is:ie,km)  :: q2, dp2
  real, dimension(is:ie,km+1):: pe1, pe2, pk1, pk2, pn2, phis
  real, dimension(is:ie+1,km+1):: pe0, pe3
  real, dimension(is:ie):: gz, cvm, qv
  real rcp, rg, tmp, tpe, rrg, bkh, dtmp, k1k
  logical:: fast_mp_consv
  integer:: i,j,k
  integer:: nt, iq, n, kmp, kp, k_next

       k1k = rdgas/cv_air   ! akap / (1.-akap) = rg/Cv=0.4
        rg = rdgas
       rcp = 1./ cp
       rrg = -rdgas/grav

!$OMP parallel do default(none) shared(is,ie,js,je,km,pe,ptop,kord_tm,hydrostatic, &
!$OMP                                  pt,pk,rg,peln,q,nwat,    &
!$OMP                                  sphum,r_vir,rcp,k1k,delp, &
!$OMP                                  delz,akap,pkz,te,u,v,ps, gridstruct, last_step, &
!$OMP                                  ak,bk,nq,isd,ied,jsd,jed,kord_tr, &
!$OMP                                  hs,w,ws,kord_wz,do_omega,omga,rrg,kord_mt,ua)    &
!$OMP                          private(qv,gz,cvm,kp,k_next,bkh,dp2,   &
!$OMP                                  pe0,pe1,pe2,pe3,pk1,pk2,pn2,phis,q2)
  do 1000 j=js,je+1

     do k=1,km+1
        do i=is,ie
           pe1(i,k) = pe(i,k,j)
        enddo
     enddo

     do i=is,ie
        pe2(i,   1) = ptop
        pe2(i,km+1) = pe(i,km+1,j)
     enddo

  if ( j /= (je+1) ) then
       if ( kord_tm < 0 ) then
! Note: pt at this stage is Theta_v
             if ( hydrostatic ) then
! Transform virtual pt to virtual Temp
             do k=1,km
                   do i=is,ie
                      pt(i,j,k) = pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
                   enddo
             enddo
             else
! Transform "density pt" to "density temp"
               do k=1,km
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
! Using dry pressure for the definition of the virtual potential temperature
!                    pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*    &
!                                              pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
                  enddo
               enddo
             endif         ! hydro test
       elseif ( hydrostatic ) then
           call pkez(km, is, ie, js, je, j, pe, pk, akap, peln, pkz, ptop)
! Compute cp*T + KE
           do k=1,km
                 do i=is,ie
                    te(i,j,k) = 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                 v(i,j,k)**2+v(i+1,j,k)**2 -  &
                               (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))  &
                              + cp_air*pt(i,j,k)*pkz(i,j,k)
                 enddo
           enddo
       endif

     if ( .not. hydrostatic ) then
           do k=1,km
              do i=is,ie
                 delz(i,j,k) = -delz(i,j,k) / delp(i,j,k) ! ="specific volume"/grav
              enddo
           enddo
      endif

! update ps
      do i=is,ie
         ps(i,j) = pe1(i,km+1)
      enddo
!
! Hybrid sigma-P coordinate:
!
        do k=2,km
           do i=is,ie
              pe2(i,k) = ak(k) + bk(k)*pe(i,km+1,j)
           enddo
        enddo
        do k=1,km
           do i=is,ie
              dp2(i,k) = pe2(i,k+1) - pe2(i,k)
           enddo
        enddo

!------------
! update delp
!------------
      do k=1,km
         do i=is,ie
            delp(i,j,k) = dp2(i,k)
         enddo
      enddo

!------------------
! Compute p**Kappa
!------------------
   do k=1,km+1
      do i=is,ie
         pk1(i,k) = pk(i,j,k)
      enddo
   enddo

   do i=is,ie
      pn2(i,   1) = peln(i,   1,j)
      pn2(i,km+1) = peln(i,km+1,j)
      pk2(i,   1) = pk1(i,   1)
      pk2(i,km+1) = pk1(i,km+1)
   enddo

   do k=2,km
      do i=is,ie
         pn2(i,k) = log(pe2(i,k))
         pk2(i,k) = exp(akap*pn2(i,k))
      enddo
   enddo

   if ( kord_tm<0 ) then
!----------------------------------
! Map t using logp
 !----------------------------------
         call map_scalar(km,  peln(is,1,j),  pt, gz,   &
                         km,  pn2,           pt,              &
                         is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm), t_min)
   else
! Map pt using pe
         call map1_ppm (km,  pe1,  pt,  gz,       &
                        km,  pe2,  pt,                  &
                        is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm))
   endif

!----------------
! Map constituents
!----------------
      if ( nq > 0 ) then
! Remap one tracer at a time
         do iq=1,nq
             call map1_q2(km, pe1, q(isd,jsd,1,iq),     &
                          km, pe2, q2, dp2,             &
                          is, ie, 0, kord_tr(iq), j, isd, ied, jsd, jed, 0.)
            do k=1,km
               do i=is,ie
                  q(i,j,k,iq) = q2(i,k)
               enddo
            enddo
         enddo
      endif

   if ( .not. hydrostatic ) then
! Remap vertical wind:
        call map1_ppm (km,   pe1,  w,  ws(is,j),   &
                       km,   pe2,  w,              &
                       is, ie, j, isd, ied, jsd, jed, -2, kord_wz)
! Remap delz for hybrid sigma-p coordinate
        call map1_ppm (km,   pe1, delz,  gz,   &
                       km,   pe2, delz,              &
                       is, ie, j, isd,  ied,  jsd,  jed,  1, abs(kord_tm))
        do k=1,km
           do i=is,ie
              delz(i,j,k) = -delz(i,j,k)*dp2(i,k)
           enddo
        enddo
   endif

!----------
! Update pk
!----------
   do k=1,km+1
      do i=is,ie
         pk(i,j,k) = pk2(i,k)
      enddo
   enddo

!----------------
   if ( do_omega ) then
! Start do_omega
! Copy omega field to pe3
      do i=is,ie
         pe3(i,1) = 0.
      enddo
      do k=2,km+1
         do i=is,ie
            pe3(i,k) = omga(i,j,k-1)
         enddo
      enddo
   endif

   do k=1,km+1
      do i=is,ie
          pe0(i,k)   = peln(i,k,j)
         peln(i,k,j) =  pn2(i,k)
      enddo
   enddo

!------------
! Compute pkz
!------------
   if ( hydrostatic ) then
      do k=1,km
         do i=is,ie
            pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
         enddo
      enddo
   else
! Note: pt at this stage is T_v or T_m
         do k=1,km
         if ( kord_tm < 0 ) then
           do i=is,ie
              pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
! Using dry pressure for the definition of the virtual potential temperature
!             pkz(i,j,k) = exp(akap*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
           enddo
         else
           do i=is,ie
              pkz(i,j,k) = exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
! Using dry pressure for the definition of the virtual potential temperature
!             pkz(i,j,k) = exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
           enddo
           if ( last_step) then
              do i=is,ie
                 pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)
              enddo
           endif
         endif
         enddo
   endif

! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
   if ( do_omega ) then
   do k=1,km
      do i=is,ie
         dp2(i,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))
      enddo
   enddo
   do i=is,ie
       k_next = 1
       do n=1,km
          kp = k_next
          do k=kp,km
             if( dp2(i,n) <= pe0(i,k+1) .and. dp2(i,n) >= pe0(i,k) ) then
                 omga(i,j,n) = pe3(i,k)  +  (pe3(i,k+1) - pe3(i,k)) *    &
                       (dp2(i,n)-pe0(i,k)) / (pe0(i,k+1)-pe0(i,k) )
                 k_next = k
                 exit
             endif
          enddo
       enddo
   enddo
   endif     ! end do_omega

  endif !(j < je+1)

      do i=is,ie+1
         pe0(i,1) = pe(i,1,j)
      enddo
!------
! map u
!------
      do k=2,km+1
         do i=is,ie
            pe0(i,k) = 0.5*(pe(i,k,j-1)+pe1(i,k))
         enddo
      enddo

      do k=1,km+1
         bkh = 0.5*bk(k)
         do i=is,ie
            pe3(i,k) = ak(k) + bkh*(pe(i,km+1,j-1)+pe1(i,km+1))
         enddo
      enddo

      call map1_ppm( km, pe0(is:ie,:),   u,   gz,   &
                     km, pe3(is:ie,:),   u,               &
                     is, ie, j, isd, ied, jsd, jed+1, -1, kord_mt)

   if (j < je+1) then
!------
! map v
!------
       do i=is,ie+1
          pe3(i,1) = ak(1)
       enddo

       do k=2,km+1
          bkh = 0.5*bk(k)
          do i=is,ie+1
             pe0(i,k) =         0.5*(pe(i-1,k,   j)+pe(i,k,   j))
             pe3(i,k) = ak(k) + bkh*(pe(i-1,km+1,j)+pe(i,km+1,j))
          enddo
       enddo

       call map1_ppm (km, pe0,  v, gz,    &
                      km, pe3,  v, is, ie+1,    &
                      j, isd, ied+1, jsd, jed, -1, kord_mt)
   endif ! (j < je+1)

     do k=1,km
        do i=is,ie
           ua(i,j,k) = pe2(i,k+1)
        enddo
     enddo

1000  continue

!$OMP parallel default(none) shared(is,ie,js,je,km,kmp,ptop,u,v,pe,ua,isd,ied,jsd,jed,kord_mt, &
!$OMP                               te_2d,te,delp,hydrostatic,hs,rg,pt,peln, &
!$OMP                               cp,delz,nwat,       &
!$OMP                               r_vir,sphum,w,pk,pkz,last_step,consv, &
!$OMP                               zsum1,zsum0,te0_2d,               &
!$OMP                               ng,gridstruct,E_Flux,pdt,dtmp,q,      &
!$OMP                               rrg,akap,&
!$OMP                               fast_mp_consv,kord_tm) &
!$OMP                       private(pe0,pe1,pe2,pe3,qv,cvm,gz,phis,tpe,tmp, dpln)

!$OMP do
  do k=2,km
     do j=js,je
        do i=is,ie
           pe(i,k,j) = ua(i,j,k-1)
        enddo
     enddo
  enddo

dtmp = 0.
if( last_step ) then

  if ( consv > consv_min ) then

!$OMP do
    do j=js,je
       if ( hydrostatic ) then
            do i=is,ie
               gz(i) = hs(i,j)
               do k=1,km
                  gz(i) = gz(i) + rg*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
               enddo
            enddo
            do i=is,ie
               te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i)
            enddo

            do k=1,km
            do i=is,ie
               te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*pt(i,j,k) +   &
                            0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                        v(i,j,k)**2+v(i+1,j,k)**2 -  &
                           (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)))
            enddo
            enddo
       else
           do i=is,ie
              te_2d(i,j) = 0.
              phis(i,km+1) = hs(i,j)
           enddo
           do k=km,1,-1
              do i=is,ie
                 phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
              enddo
           enddo

           do k=1,km
              do i=is,ie
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cv_air*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) + &
                                 0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                                 u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                                (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
              enddo
           enddo   ! k-loop
       endif  ! end non-hydro

         do i=is,ie
            te_2d(i,j) = te0_2d(i,j) - te_2d(i,j)
            zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
         enddo
         do k=2,km
            do i=is,ie
               zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
            enddo
         enddo
         if ( hydrostatic ) then
            do i=is,ie
               zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
            enddo
         endif

    enddo   ! j-loop

!$OMP single
         tpe = consv*g_sum(te_2d, is, ie, js, je, ng, gridstruct%area)
      E_Flux = tpe / (grav*pdt*4.*pi*radius**2)    ! unit: W/m**2
                                                   ! Note pdt is "phys" time step
      if ( hydrostatic ) then
           dtmp = tpe / (cp*g_sum(zsum0,  is, ie, js, je, ng, gridstruct%area))
      else
           dtmp = tpe / (cv_air*g_sum(zsum1, is, ie, js, je, ng, gridstruct%area))
      endif
!$OMP end single

  elseif ( consv < -consv_min ) then

!$OMP do
      do j=js,je
         do i=is,ie
            zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
         enddo
         do k=2,km
            do i=is,ie
               zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
            enddo
         enddo
         if ( hydrostatic ) then
            do i=is,ie
               zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
            enddo
         endif
      enddo

      E_Flux = consv
!$OMP single
      if ( hydrostatic ) then
           dtmp = E_flux*(grav*pdt*4.*pi*radius**2) /    &
                 (cp*g_sum(zsum0,  is, ie, js, je, ng, gridstruct%area))
      else
           dtmp = E_flux*(grav*pdt*4.*pi*radius**2) /    &
                 (cv_air*g_sum(zsum1,  is, ie, js, je, ng, gridstruct%area))
      endif
!$OMP end single
  endif        ! end consv check
endif        ! end last_step check

! Note: pt at this stage is T_v

  if ( last_step ) then
       ! Output temperature if last_step
!$OMP do
        do k=1,km
           do j=js,je
                 do i=is,ie
                    pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k)) / (1.+r_vir*q(i,j,k,sphum))
                 enddo
           enddo   ! j-loop
        enddo  ! k-loop
  else  ! not last_step
    if ( kord_tm < 0 ) then
!$OMP do
       do k=1,km
          do j=js,je
             do i=is,ie
                pt(i,j,k) = pt(i,j,k)/pkz(i,j,k)
             enddo
          enddo
       enddo
    endif
  endif
!$OMP end parallel

 end subroutine Lagrangian_to_Eulerian


  subroutine pkez(km, ifirst, ilast, jfirst, jlast, j, &
                  pe, pk, akap, peln, pkz, ptop)

! !INPUT PARAMETERS:
   integer, intent(in):: km, j
   integer, intent(in):: ifirst, ilast        ! Latitude strip
   integer, intent(in):: jfirst, jlast        ! Latitude strip
   real, intent(in):: akap
   real, intent(in):: pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1)
   real, intent(in):: pk(ifirst:ilast,jfirst:jlast,km+1)
   real, intent(IN):: ptop
! !OUTPUT
   real, intent(out):: pkz(ifirst:ilast,jfirst:jlast,km)
   real, intent(inout):: peln(ifirst:ilast, km+1, jfirst:jlast)   ! log (pe)
! Local
   real pk2(ifirst:ilast, km+1)
   real pek
   real lnp
   real ak1
   integer i, k

   ak1 = (akap + 1.) / akap

        pek = pk(ifirst,j,1)
        do i=ifirst, ilast
           pk2(i,1) = pek
        enddo

        do k=2,km+1
           do i=ifirst, ilast
!             peln(i,k,j) =  log(pe(i,k,j))
              pk2(i,k) =  pk(i,j,k)
           enddo
        enddo

!---- GFDL modification
       if( ptop < ptop_min ) then
           do i=ifirst, ilast
               peln(i,1,j) = peln(i,2,j) - ak1
           enddo
       else
           lnp = log( ptop )
           do i=ifirst, ilast
              peln(i,1,j) = lnp
           enddo
       endif
!---- GFDL modification

       do k=1,km
          do i=ifirst, ilast
             pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k) )  /  &
                          (akap*(peln(i,k+1,j) - peln(i,k,j)) )
          enddo
       enddo

 end subroutine pkez



 subroutine map_scalar( km,   pe1,    q1,   qs,           &
                        kn,   pe2,    q2,   i1, i2,       &
                         j,  ibeg, iend, jbeg, jend, iv,  kord, q_min)
! iv=1
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == temp
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real, intent(in) ::   qs(i1:i2)       ! bottom BC
 real, intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real, intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
 real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output
 real, intent(in):: q_min

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real    dp1(i1:i2,km)
   real   q4(4,i1:i2,km)
   real    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map_scalar


 subroutine map1_ppm( km,   pe1,    q1,   qs,           &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord)
 integer, intent(in) :: i1                ! Starting longitude
 integer, intent(in) :: i2                ! Finishing longitude
 integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == ???
                                          !       2 == remap temp with cs scheme
 integer, intent(in) :: kord              ! Method order
 integer, intent(in) :: j                 ! Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                ! Original vertical dimension
 integer, intent(in) :: kn                ! Target vertical dimension
 real, intent(in) ::   qs(i1:i2)       ! bottom BC
 real, intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges
                                       ! (from model top to bottom surface)
                                       ! in the original vertical coordinate
 real, intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges
                                       ! (from model top to bottom surface)
                                       ! in the new vertical coordinate
 real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output

! !DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! !LOCAL VARIABLES:
   real    dp1(i1:i2,km)
   real   q4(4,i1:i2,km)
   real    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm


 subroutine map1_q2(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend, q_min )


! !INPUT PARAMETERS:
      integer, intent(in) :: j
      integer, intent(in) :: i1, i2
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents 1 == ???
      integer, intent(in) :: kord
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real, intent(in) ::  pe1(i1:i2,km+1)     ! pressure at layer edges
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real, intent(in) ::  pe2(i1:i2,kn+1)     ! pressure at layer edges
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real, intent(in) ::  q1(ibeg:iend,jbeg:jend,km) ! Field input
      real, intent(in) ::  dp2(i1:i2,kn)
      real, intent(in) ::  q_min
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: q2(i1:i2,kn) ! Field output
! !LOCAL VARIABLES:
      real   qs(i1:i2)
      real   dp1(i1:i2,km)
      real   q4(4,i1:i2,km)
      real   pl, pr, qsum, dp, esl

      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                    esl = dp / dp1(i,m)
                   qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                       (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                   k0 = m
                   goto 123
                 endif
              enddo
              goto 123
          endif
      endif
100   continue
123   q2(i,k) = qsum / dp2(i,k)
555   continue
1000  continue

 end subroutine map1_q2


 subroutine scalar_profile(qs, a4, delp, km, i1, i2, iv, kord, qmin)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real, intent(in)   ::   qs(i1:i2)
 real, intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real, intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
 real, intent(in):: qmin
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext6
 real  gam(i1:i2,km)
 real    q(i1:i2,km+1)
 real   d4(i1:i2)
 real   bet, a_bot, grat
 real   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo

  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------
!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
               q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord)==16 ) then
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          ext6(i,k) = abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k))
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. a4(1,i,k)<qmin ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              if( a4(1,i,k)<qmin .or. extm(i,k-1) .or. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected; or q is too small -> ehance vertical mixing
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
                   a4(4,i,k) = 0.
              else
! True local extremum
                a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
              endif
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( extm(i,k) ) then
             if ( extm(i,k-1) .and. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
                 a4(4,i,k) = 0.
             else
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
             endif
          else
             a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          endif
       enddo
     elseif ( abs(kord)==14 ) then
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==16 ) then
       do i=i1,i2
          if( ext6(i,k) ) then
             if ( extm(i,k-1) .or. extm(i,k+1) ) then
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
             endif
          endif
       enddo
     else      ! kord = 11, 13
       do i=i1,i2
         if ( extm(i,k) .and. (extm(i,k-1).or.extm(i,k+1).or.a4(1,i,k)<qmin) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine scalar_profile


 subroutine cs_profile(qs, a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: kord
 real, intent(in)   ::   qs(i1:i2)
 real, intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real, intent(inout):: a4(4,i1:i2,km)     ! Interpolated values
!-----------------------------------------------------------------------
 logical:: extm(i1:i2,km)
 real  gam(i1:i2,km)
 real    q(i1:i2,km+1)
 real   d4(i1:i2)
 real   bet, a_bot, grat
 real   pmp_1, lac_1, pmp_2, lac_2
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo

  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif
!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              if( extm(i,k-1) .or. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
                   a4(4,i,k) = 0.
              else
! True local extremum
                a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
              endif
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( extm(i,k) ) then
             if ( extm(i,k-1) .and. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
                 a4(4,i,k) = 0.
             else
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
             endif
          else
             a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          endif
       enddo
     elseif ( abs(kord)==14 ) then
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11
       do i=i1,i2
         if ( extm(i,k) .and. (extm(i,k-1) .or. extm(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile


 subroutine cs_limiters(im, extm, a4, iv)
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real , intent(inout) :: a4(4,im)   ! PPM array
! !LOCAL VARIABLES:
 real  da1, da2, a6da
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters



 subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
                               ! iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               !
 real , intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real , intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
!
! !REVISION HISTORY:
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real    dc(i1:i2,km)
      real    h2(i1:i2,km)
      real  delq(i1:i2,km)
      real   df2(i1:i2,km)
      real    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real  fac
      real  a1, a2, c1, c2, c3, d1, d2
      real  qm, dq, lac, qmp, pmp

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo


! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

#ifdef BOT_MONO
      do i=i1,i2
         a4(4,i,km) = 0
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
         d1 = a4(1,i,km) - a4(2,i,km)
         d2 = a4(3,i,km) - a4(1,i,km)
         if ( d1*d2 < 0. ) then
              a4(2,i,km) = a4(1,i,km)
              a4(3,i,km) = a4(1,i,km)
         else
              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
              a4(2,i,km) = a4(1,i,km) - dq
              a4(3,i,km) = a4(1,i,km) + dq
         endif
      enddo
#else
      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif
#endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile


 subroutine ppm_limiters(dm, a4, itot, lmt)

! !INPUT PARAMETERS:
      real , intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)
! !LOCAL VARIABLES:
      real  qmp
      real  da1, da2, a6da
      real  fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters


 real function g_sum(p, ifirst, ilast, jfirst, jlast, ngc, area)
! Fast version of globalsum
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast, ngc
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real, intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
      integer :: i,j
      real gsum

         gsum = 0.
         do j=jfirst,jlast
            do i=ifirst,ilast
               gsum = gsum + p(i,j)*area(i,j)
            enddo
         enddo

 end function g_sum



end module fv_mapz_mod
