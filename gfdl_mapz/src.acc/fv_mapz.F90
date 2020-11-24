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
    real, dimension(is:ie,js:je):: te_2d, zsum0, zsum1
    real, dimension(is:ie,js:je+1,km)  :: q2, dp2
    real, dimension(is:ie,js:je+1,km+1):: pe1, pe2, pk1, pk2, pn2, phis
    real, dimension(is:ie+1,js:je+1,km+1):: pe0, pe3
    real, dimension(is:ie,js:je+1):: gz
    real rcp, rg, tmp, tpe, rrg, bkh, dtmp, k1k
    logical:: fast_mp_consv
    integer:: i,j,k 
    integer:: nt, iq, n, kmp, kp, k_next

    k1k = rdgas/cv_air   ! akap / (1.-akap) = rg/Cv=0.4
    rg = rdgas
    rcp = 1./ cp
    rrg = -rdgas/grav
!$ACC kernels
!$ACC loop collapse(3)
    do k = 1, km+1
       do j=js,je+1
          do i=is,ie
             pe1(i,j,k) = pe(i,k,j)
          enddo
       enddo
    enddo
!$ACC loop collapse(2)
    do j=js,je+1
       do i=is,ie
          pe2(i,j,   1) = ptop
          pe2(i,j,km+1) = pe(i,km+1,j)
       enddo
    enddo
    if ( kord_tm < 0 ) then
       ! Note: pt at this stage is Theta_v
       if ( hydrostatic ) then
          ! Transform virtual pt to virtual Temp
          do k=1,km
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
                enddo
             enddo
          enddo
       else
          ! Transform "density pt" to "density temp"
          do k=1,km
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                enddo
             enddo
          enddo
       endif         ! hydro test
    elseif ( hydrostatic ) then
       call pkez(km, is, ie, js, je, pe, pk, akap, peln, pkz, ptop)
       ! Compute cp*T + KE
       do k=1,km
          do j=js,je
             do i=is,ie
                te(i,j,k) = 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                     v(i,j,k)**2+v(i+1,j,k)**2 -  &
                     (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))  &
                     + cp_air*pt(i,j,k)*pkz(i,j,k)
             enddo
          enddo
       enddo
    endif
    if ( .not. hydrostatic ) then
       do k=1,km
          do j=js,je
             do i=is,ie
                delz(i,j,k) = -delz(i,j,k) / delp(i,j,k) ! ="specific volume"/grav
             enddo
          enddo
       enddo
    endif

    ! update ps
    do j=js,je
       do i=is,ie
          ps(i,j) = pe1(i,j,km+1)
       enddo
    enddo
    !
    ! Hybrid sigma-P coordinate:
    !
    do k=2,km
       do j=js,je
          do i=is,ie
             pe2(i,j,k) = ak(k) + bk(k)*pe(i,km+1,j)
          enddo
       enddo
    enddo
    do k=1,km
       do j=js,je
          do i=is,ie
             dp2(i,j,k) = pe2(i,j,k+1) - pe2(i,j,k)
             !---update delp
             delp(i,j,k) = dp2(i,j,k)
          enddo
       enddo
    enddo
    !------------------
    ! Compute p**Kappa
    !------------------
    do k=1,km+1
       do j=js,je
          do i=is,ie
             pk1(i,j,k) = pk(i,j,k)
          enddo
       enddo
    enddo
    do j=js,je
       do i=is,ie
          pn2(i,j,   1) = peln(i,   1,j)
          pn2(i,j,km+1) = peln(i,km+1,j)
          pk2(i,j,   1) = pk1(i,j,   1)
          pk2(i,j,km+1) = pk1(i,j,km+1)
       enddo
    enddo
    do k=2,km
       do j=js,je
          do i=is,ie
             pn2(i,j,k) = log(pe2(i,j,k))
             pk2(i,j,k) = exp(akap*pn2(i,j,k))
          enddo
       enddo
    enddo
!$ACC end kernels
    if ( kord_tm<0 ) then
       !----------------------------------
       ! Map t using logp 
       !----------------------------------
          call map_scalar(km,  peln,  pt, gz,   &
               km,  pn2,           pt,              &
               is, ie, js, je, isd, ied, jsd, jed, 1, abs(kord_tm), t_min)
!    else
       ! Map pt using pe
!          call map1_ppm (km,  pe1,  pt,  gz,       &
!               km,  pe2,  pt,                  &
!               is, ie, js, je, isd, ied, jsd, jed, 1, abs(kord_tm))
    endif
    !----------------
    ! Map constituents
    !----------------
    ! Remap one tracer at a time
    do iq=1,nq
          call map1_q2(km, pe1, q(isd,jsd,1,iq),     &
               km, pe2, q2, dp2,             &
               is, ie, 0, kord_tr(iq), js, je, isd, ied, jsd, jed, 0.)
        
       do k=1,km
          do j=js,je
             do i=is,ie
                q(i,j,k,iq) = q2(i,j,k)
             enddo
          enddo
       enddo
    enddo

    if ( .not. hydrostatic ) then
       ! Remap vertical wind:
       call map1_ppm (km,   pe1,  w,  ws,   &
            km,   pe2,  w,              &
            is, ie, js, je, isd, ied, jsd, jed, -2, kord_wz)
       ! Remap delz for hybrid sigma-p coordinate
       call map1_ppm (km,   pe1, delz,  gz,   &
            km,   pe2, delz,              &
            is, ie, js, je, isd,  ied,  jsd,  jed,  1, abs(kord_tm))
       do k=1,km
          do j=js,je
             do i=is,ie
                delz(i,j,k) = -delz(i,j,k)*dp2(i,j,k)
             enddo
          enddo
       enddo
    endif

    !----------
    ! Update pk
    !----------
    do k=1,km+1
       do j=js,je
          do i=is,ie
             pk(i,j,k) = pk2(i,j,k)
          enddo
       enddo
    enddo
    !----------------
    if ( do_omega ) then
       ! Start do_omega
       ! Copy omega field to pe3
       do j=js,je
          do i=is,ie
             pe3(i,j,1) = 0.
          enddo
       enddo
       do k=2,km+1
          do j=js,je
             do i=is,ie
                pe3(i,j,k) = omga(i,j,k-1)
             enddo
          enddo
       enddo
    endif
    do k=1,km+1
       do j=js,je
          do i=is,ie
             pe0(i,j,k)   = peln(i,k,j)
             peln(i,k,j) =  pn2(i,j,k)
          enddo
       enddo
    enddo
    !------------
    ! Compute pkz
    !------------
    if ( hydrostatic ) then
       do k=1,km
          do j=js,je
             do i=is,ie
                pkz(i,j,k) = (pk2(i,j,k+1)-pk2(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
             enddo
          enddo
       enddo
    else
       ! Note: pt at this stage is T_v or T_m
       if ( kord_tm < 0 ) then
          do k=1,km
             do j=js,je
                do i=is,ie
                   pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                   ! Using dry pressure for the definition of the virtual potential temperature
                   !             pkz(i,j,k) = exp(akap*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
                enddo
             enddo
          enddo
       else
          do k=1,km
             do j=js,je
                do i=is,ie
                   pkz(i,j,k) = exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                   ! Using dry pressure for the definition of the virtual potential temperature
                   !             pkz(i,j,k) = exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
                enddo
             enddo
          enddo
          if ( last_step) then
             do k=1,km
                do j=js,je
                   do i=is,ie
                      pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)
                   enddo
                enddo
             enddo
          endif
       endif
    endif

    ! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
    if ( do_omega ) then
       do k=1,km
          do j=js,je
             do i=is,ie
                dp2(i,j,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))
             enddo
          enddo
       enddo
       do j=js,je
          do i=is,ie
             k_next = 1
             do n=1,km
                kp = k_next
                do k=kp,km
                   if( dp2(i,j,n) <= pe0(i,j,k+1) .and. dp2(i,j,n) >= pe0(i,j,k) ) then
                      omga(i,j,n) = pe3(i,j,k)  +  (pe3(i,j,k+1) - pe3(i,j,k)) *    &
                           (dp2(i,j,n)-pe0(i,j,k)) / (pe0(i,j,k+1)-pe0(i,j,k) )
                      k_next = k
                      exit
                   endif
                enddo
             enddo
          enddo
       enddo
    endif     ! end do_omega


    do j=js,je+1
       do i=is,ie+1
          pe0(i,j,1) = pe(i,1,j)
       enddo
    enddo
    !------
    ! map u
    !------
    do k=2,km+1
       do j=js,je+1
          do i=is,ie
             pe0(i,j,k) = 0.5*(pe(i,k,j-1)+pe1(i,j,k))
          enddo
       enddo
    enddo
    do k=1,km+1
       bkh = 0.5*bk(k)
       do j=js,je+1
          do i=is,ie
             pe3(i,j,k) = ak(k) + bkh*(pe(i,km+1,j-1)+pe1(i,j,km+1))
          enddo
       enddo
    enddo
    call map1_ppm( km, pe0,   u,   gz,   &
            km, pe3,   u,               &
            is, ie, js, je+1, isd, ied, jsd, jed+1, -1, kord_mt)
    !------
    ! map v
    !------
    do j=js,je
       do i=is,ie+1
          pe3(i,j,1) = ak(1)
       enddo
    enddo
    
    do k=2,km+1
       bkh = 0.5*bk(k)
       do j=js,je+1
          do i=is,ie+1
             pe0(i,j,k) =         0.5*(pe(i-1,k,   j)+pe(i,k,   j))
             pe3(i,j,k) = ak(k) + bkh*(pe(i-1,km+1,j)+pe(i,km+1,j))
          enddo
       enddo
    enddo
    call map1_ppm (km, pe0,  v, gz,    &
            km, pe3,  v, is, ie+1,    &
            js, je, isd, ied+1, jsd, jed, -1, kord_mt)
    do k=1,km
       do j=js,je+1
          do i=is,ie
             ua(i,j,k) = pe2(i,j,k+1)
          enddo
       enddo
    enddo

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
          
          if ( hydrostatic ) then
             do j=js,je
                do i=is,ie
                   gz(i,j) = hs(i,j)
                enddo
             enddo
             do k=1,km
                do j=js,je
                   do i=is,ie
                      gz(i,j) = gz(i,j) + rg*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
                   enddo
                enddo
             enddo
             do j=js,je
                do i=is,ie
                   te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i,j)
                enddo
             enddo

             do k=1,km
                do j=js,je
                   do i=is,ie
                      te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*pt(i,j,k) +   &
                           0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                           v(i,j,k)**2+v(i+1,j,k)**2 -  &
                           (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)))
                   enddo
                enddo
             enddo
          else
             do j=js,je
                do i=is,ie
                   te_2d(i,j) = 0.
                   phis(i,j,km+1) = hs(i,j)
                enddo
             enddo
             do k=km,1,-1
                do j=js,je
                   do i=is,ie
                      phis(i,j,k) = phis(i,j,k+1) - grav*delz(i,j,k)
                   enddo
                enddo
             enddo

             do k=1,km
                do j=js,je
                   do i=is,ie
                      te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cv_air*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) + &
                           0.5*(phis(i,j,k)+phis(i,j,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                           u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                           (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
                   enddo
                enddo   ! k-loop
             enddo
          endif  ! end non-hydro
          do j=js,je
             do i=is,ie
                te_2d(i,j) = te0_2d(i,j) - te_2d(i,j)
                zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
             enddo
          enddo
          do k=2,km
             do j=js,je
                do i=is,ie
                   zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
                enddo
             enddo
          enddo
          if ( hydrostatic ) then
             do j=js,je
                do i=is,ie
                   zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
                enddo
             enddo
          endif

          tpe = consv*g_sum(te_2d, is, ie, js, je, ng, gridstruct%area)
          E_Flux = tpe / (grav*pdt*4.*pi*radius**2)    ! unit: W/m**2
          ! Note pdt is "phys" time step
          if ( hydrostatic ) then
             dtmp = tpe / (cp*g_sum(zsum0,  is, ie, js, je, ng, gridstruct%area))
          else
             dtmp = tpe / (cv_air*g_sum(zsum1, is, ie, js, je, ng, gridstruct%area))
          endif

       elseif ( consv < -consv_min ) then

          do j=js,je
             do i=is,ie
                zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
             enddo
          enddo
          do k=2,km
             do j=js,je
                do i=is,ie
                   zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
                enddo
             enddo
          enddo
          if ( hydrostatic ) then
             do j=js,je
                do i=is,ie
                   zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
                enddo
             enddo
          endif

          E_Flux = consv
          if ( hydrostatic ) then
             dtmp = E_flux*(grav*pdt*4.*pi*radius**2) /    &
                  (cp*g_sum(zsum0,  is, ie, js, je, ng, gridstruct%area))
          else
             dtmp = E_flux*(grav*pdt*4.*pi*radius**2) /    &
                  (cv_air*g_sum(zsum1,  is, ie, js, je, ng, gridstruct%area))
          endif
       endif        ! end consv check
    endif        ! end last_step check

    ! Note: pt at this stage is T_v

    if ( last_step ) then
       ! Output temperature if last_step
       do k=1,km
          do j=js,je
             do i=is,ie
                pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k)) / (1.+r_vir*q(i,j,k,sphum))
             enddo
          enddo   ! j-loop
       enddo  ! k-loop
    else  ! not last_step
       if ( kord_tm < 0 ) then
          do k=1,km
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = pt(i,j,k)/pkz(i,j,k)
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine Lagrangian_to_Eulerian


  subroutine pkez(km, ifirst, ilast, jfirst, jlast, &
       pe, pk, akap, peln, pkz, ptop)
!$ACC routine seq
    ! !INPUT PARAMETERS:
    integer, intent(in):: km
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
    integer i, j, k

    ak1 = (akap + 1.) / akap
    do j=jfirst,jlast
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
    enddo

  end subroutine pkez



  subroutine map_scalar( km,   pe1,    q1,   qs,           &
       kn,   pe2,    q2,   i1, i2, j1, j2,      &
       ibeg, iend, jbeg, jend, iv,  kord, q_min)

    ! iv=1
    integer, intent(in) :: i1                ! Starting longitude
    integer, intent(in) :: i2                ! Finishing longitude
    integer, intent(in) :: j1                ! Starting longitude
    integer, intent(in) :: j2                ! Finishing longitude
    integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == temp
    !       2 == remap temp with cs scheme
    integer, intent(in) :: kord              ! Method order
    integer, intent(in) :: ibeg, iend, jbeg, jend
    integer, intent(in) :: km                ! Original vertical dimension
    integer, intent(in) :: kn                ! Target vertical dimension
    real, intent(in) ::   qs(i1:i2,j1:j2+1)       ! bottom BC
    real, intent(in) ::  pe1(i1:i2,km+1,j1:j2)  ! pressure at layer edges 
    ! (from model top to bottom surface)
    ! in the original vertical coordinate
    real, intent(in) ::  pe2(i1:i2,j1:j2+1,kn+1)  ! pressure at layer edges 
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
    real    dp1(i1:i2,j1:j2,km)
    real   q4(4,i1:i2,j1:j2,km)
    real    pl, pr, qsum, dp, esl
    integer i, j, k, l, m, k0
    logical keep_going

!$ACC kernels
    do k=1,km
       do j=j1, j2
          do i=i1,i2
             dp1(i,j,k) = pe1(i,k+1,j) - pe1(i,k,j)
             q4(1,i,j,k) = q1(i,j,k)
          enddo
       enddo
    enddo
!$ACC end kernels
    ! Compute vertical subgrid distribution
    if ( kord >7 ) then
       call scalar_profile( qs, q4, dp1, km, i1, i2, j1, j2, iv, kord, q_min )
    else
       call ppm_profile( q4, dp1, km, i1, i2, j1, j2, iv, kord )
    endif
!$ACC kernels
    do j=j1,j2
       do i=i1,i2
          k0 = 1
          !$acc loop seq
          do k=1,kn
             keep_going = .true.
             do l=k0,km
                ! locate the top edge: pe2(i,j,k)
                if( keep_going .and. pe2(i,j,k) >= pe1(i,l,j) .and. pe2(i,j,k) <= pe1(i,l+1,j) ) then
                   pl = (pe2(i,j,k)-pe1(i,l,j)) / dp1(i,j,l)
                   if( pe2(i,j,k+1) <= pe1(i,l+1,j) ) then
                      ! entire new grid is within the original grid
                      pr = (pe2(i,j,k+1)-pe1(i,l,j)) / dp1(i,j,l)
                      q2(i,j,k) = q4(2,i,j,l) + 0.5*(q4(4,i,j,l)+q4(3,i,j,l)-q4(2,i,j,l))  &
                           *(pr+pl)-q4(4,i,j,l)*r3*(pr*(pr+pl)+pl**2)
                      k0 = l
                   else
                      ! Fractional area...
                      qsum = (pe1(i,l+1,j)-pe2(i,j,k))*(q4(2,i,j,l)+0.5*(q4(4,i,j,l)+   &
                           q4(3,i,j,l)-q4(2,i,j,l))*(1.+pl)-q4(4,i,j,l)*           &
                           (r3*(1.+pl*(1.+pl))))
                      do m=l+1,km
                         ! locate the bottom edge: pe2(i,k+1)
                         if( pe2(i,j,k+1) > pe1(i,m+1,j) ) then
                            ! Whole layer
                            qsum = qsum + dp1(i,j,m)*q4(1,i,j,m)
                         else
                            dp = pe2(i,j,k+1)-pe1(i,m,j)
                            esl = dp / dp1(i,j,m)
                            qsum = qsum + dp*(q4(2,i,j,m)+0.5*esl*               &
                                 (q4(3,i,j,m)-q4(2,i,j,m)+q4(4,i,j,m)*(1.-r23*esl)))
                            k0 = m
                            exit
                         endif
                      enddo
                      q2(i,j,k) = qsum / ( pe2(i,j,k+1) - pe2(i,j,k) )
                   endif
                   keep_going = .false.
                endif
             enddo
          enddo
       enddo
    enddo
!$ACC end kernels

  end subroutine map_scalar


  subroutine map1_ppm( km,   pe1,    q1,   qs,           &
       kn,   pe2,    q2,   i1, i2,       &
       j1, j2, ibeg, iend, jbeg, jend, iv,  kord)
    integer, intent(in) :: i1                ! Starting longitude
    integer, intent(in) :: i2                ! Finishing longitude
    integer, intent(in) :: j1                ! Starting latitude
    integer, intent(in) :: j2                ! Finishing latitude
    integer, intent(in) :: iv                ! Mode: 0 == constituents  1 == ???
    !       2 == remap temp with cs scheme
    integer, intent(in) :: kord              ! Method order
    integer, intent(in) :: ibeg, iend, jbeg, jend
    integer, intent(in) :: km                ! Original vertical dimension
    integer, intent(in) :: kn                ! Target vertical dimension
    real, intent(in) ::   qs(i1:,j1:)       ! bottom BC
    real, intent(in) ::  pe1(i1:,j1:,:)  ! pressure at layer edges 
    ! (from model top to bottom surface)
    ! in the original vertical coordinate
    real, intent(in) ::  pe2(i1:,j1:,:)  ! pressure at layer edges 
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
    real    dp1(i1:i2,j1:j2,km)
    real   q4(4,i1:i2,j1:j2,km)
    real    pl, pr, qsum, dp, esl
    integer i, j, k, l, m, k0
    logical :: keep_going 

    do k=1,km
       do j=j1,j2
          do i=i1,i2
             dp1(i,j,k) = pe1(i,j,k+1) - pe1(i,j,k)
             q4(1,i,j,k) = q1(i,j,k)
          enddo
       enddo
    enddo

    ! Compute vertical subgrid distribution
    if ( kord >7 ) then
       call  cs_profile( qs, q4, dp1, km, i1, i2, j1, j2, iv, kord )
    else
       call ppm_profile( q4, dp1, km, i1, i2, j1, j2, iv, kord )
    endif

    do j=j1,j2
       do i=i1,i2
          k0 = 1
          do k=1,kn
             keep_going = .true.
             do l=k0,km
                ! locate the top edge: pe2(i,j,k)
                if( pe2(i,j,k) >= pe1(i,j,l) .and. pe2(i,j,k) <= pe1(i,j,l+1) ) then
                   pl = (pe2(i,j,k)-pe1(i,j,l)) / dp1(i,j,l)
                   if( pe2(i,j,k+1) <= pe1(i,j,l+1) ) then
                      ! entire new grid is within the original grid
                      pr = (pe2(i,j,k+1)-pe1(i,j,l)) / dp1(i,j,l)
                      q2(i,j,k) = q4(2,i,j,l) + 0.5*(q4(4,i,j,l)+q4(3,i,j,l)-q4(2,i,j,l))  &
                           *(pr+pl)-q4(4,i,j,l)*r3*(pr*(pr+pl)+pl**2)
                      k0 = l
                   else
                      ! Fractional area...
                      qsum = (pe1(i,j,l+1)-pe2(i,j,k))*(q4(2,i,j,l)+0.5*(q4(4,i,j,l)+   &
                           q4(3,i,j,l)-q4(2,i,j,l))*(1.+pl)-q4(4,i,j,l)*           &
                           (r3*(1.+pl*(1.+pl))))
                      do m=l+1,km
                         ! locate the bottom edge: pe2(i,j,k+1)
                         if( pe2(i,j,k+1) > pe1(i,j,m+1) ) then
                            ! Whole layer
                            qsum = qsum + dp1(i,j,m)*q4(1,i,j,m)
                         else
                            dp = pe2(i,j,k+1)-pe1(i,j,m)
                            esl = dp / dp1(i,j,m)
                            qsum = qsum + dp*(q4(2,i,j,m)+0.5*esl*               &
                                 (q4(3,i,j,m)-q4(2,i,j,m)+q4(4,i,j,m)*(1.-r23*esl)))
                            k0 = m
                            exit
                         endif
                      enddo
                      q2(i,j,k) = qsum / ( pe2(i,j,k+1) - pe2(i,j,k) )
                   endif
                   keep_going = .false.
                endif
             enddo
          enddo
       enddo
    enddo

  end subroutine map1_ppm


  subroutine map1_q2(km,   pe1,   q1,            &
       kn,   pe2,   q2,   dp2,     &
       i1,   i2,    iv,   kord, j1, j2, &
       ibeg, iend, jbeg, jend, q_min )

    ! !INPUT PARAMETERS:
    integer, intent(in) :: j1, j2
    integer, intent(in) :: i1, i2
    integer, intent(in) :: ibeg, iend, jbeg, jend
    integer, intent(in) :: iv                ! Mode: 0 ==  constituents 1 == ???
    integer, intent(in) :: kord
    integer, intent(in) :: km                ! Original vertical dimension
    integer, intent(in) :: kn                ! Target vertical dimension

    real, intent(in) ::  pe1(i1:i2,j1:j2+1, km+1)     ! pressure at layer edges 
    ! (from model top to bottom surface)
    ! in the original vertical coordinate
    real, intent(in) ::  pe2(i1:i2,j1:j2+1,kn+1)     ! pressure at layer edges 
    ! (from model top to bottom surface)
    ! in the new vertical coordinate
    real, intent(in) ::  q1(ibeg:iend,jbeg:jend,km) ! Field input
    real, intent(in) ::  dp2(i1:i2,j1:j2+1,kn)
    real, intent(in) ::  q_min
    ! !INPUT/OUTPUT PARAMETERS:
    real, intent(inout):: q2(i1:i2,j1:j2+1,kn) ! Field output
    ! !LOCAL VARIABLES:
    real   qs(i1:i2,j1:j2)
    real   dp1(i1:i2,j1:j2,km)
    real   q4(4,i1:i2,j1:j2,km)
    real   pl, pr, qsum, dp, esl

    integer i, j, k, l, m, k0
    logical :: keep_going

    do j=j1,j2
       do k=1,km
          do i=i1,i2
             dp1(i,j,k) = pe1(i,j,k+1) - pe1(i,j,k)
             q4(1,i,j,k) = q1(i,j,k)
          enddo
       enddo
    enddo
    ! Compute vertical subgrid distribution
    if ( kord >7 ) then
       call  scalar_profile( qs, q4, dp1, km, i1, i2, j1, j2, iv, kord, q_min )
    else
       call ppm_profile( q4, dp1, km, i1, i2, j1, j2, iv, kord )
    endif

    ! Mapping
    do j=j1,j2
       do i=i1,i2
          k0 = 1
          do k=1,kn
             keep_going = .true.
             do l=k0,km
                ! locate the top edge: pe2(i,j,k)
                if(keep_going .and. pe2(i,j,k) >= pe1(i,j,l) .and. pe2(i,j,k) <= pe1(i,j,l+1)) then
                   pl = (pe2(i,j,k)-pe1(i,j,l)) / dp1(i,j,l)
                   if(pe2(i,j,k+1) <= pe1(i,j,l+1)) then
                      ! entire new grid is within the original grid
                      pr = (pe2(i,j,k+1)-pe1(i,j,l)) / dp1(i,j,l)
                      q2(i,j,k) = q4(2,i,j,l) + 0.5*(q4(4,i,j,l)+q4(3,i,j,l)-q4(2,i,j,l))  &
                           *(pr+pl)-q4(4,i,j,l)*r3*(pr*(pr+pl)+pl**2)
                      k0 = l
                   else
                      ! Fractional area...
                      qsum = (pe1(i,j,l+1)-pe2(i,j,k))*(q4(2,i,j,l)+0.5*(q4(4,i,j,l)+   &
                           q4(3,i,j,l)-q4(2,i,j,l))*(1.+pl)-q4(4,i,j,l)*           &
                           (r3*(1.+pl*(1.+pl))))
                      do m=l+1,km
                         ! locate the bottom edge: pe2(i,j,k+1)
                         if(pe2(i,j,k+1) > pe1(i,j,m+1) ) then
                            ! Whole layer..
                            qsum = qsum + dp1(i,j,m)*q4(1,i,j,m)
                         else
                            dp = pe2(i,j,k+1)-pe1(i,j,m)
                            esl = dp / dp1(i,j,m)
                            qsum = qsum + dp*(q4(2,i,j,m)+0.5*esl*               &
                                 (q4(3,i,j,m)-q4(2,i,j,m)+q4(4,i,j,m)*(1.-r23*esl)))
                            k0 = m
                            !                   q2(i,j,k) = qsum / dp2(i,j,k)
                            exit
                         endif
                      enddo
                      q2(i,j,k) = qsum / dp2(i,j,k)
                   endif
                   keep_going = .false.
                endif
             enddo
          enddo
       enddo
    enddo

  end subroutine map1_q2


 subroutine scalar_profile(qs, a4, delp, km, i1, i2, j1, j2, iv, kord, qmin)
   ! Optimized vertical profile reconstruction:
   ! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
   integer, intent(in):: i1, i2, j1, j2
   integer, intent(in):: km      ! vertical dimension
   integer, intent(in):: iv      ! iv =-1: winds
   ! iv = 0: positive definite scalars
   ! iv = 1: others
   integer, intent(in):: kord
   real, intent(in)   ::   qs(i1:,j1:)
   real, intent(in)   :: delp(i1:i2,j1:j2,km)     ! layer pressure thickness
   real, intent(inout):: a4(4,i1:i2,j1:j2,km)     ! Interpolated values
   real, intent(in):: qmin
   !-----------------------------------------------------------------------
   logical, dimension(i1:i2,j1:j2,km):: extm, ext6
   real  gam(i1:i2,j1:j2,km)
   real    q(i1:i2,j1:j2,km+1)
   real   d4(i1:i2,j1:j2)
   real   bet, a_bot, grat 
   real   pmp_1, lac_1, pmp_2, lac_2
   integer i, j, k, im

!$acc kernels
   if ( iv .eq. -2 ) then
      do j=j1,j2
         do i=i1,i2
            gam(i,j,2) = 0.5
            q(i,j,1) = 1.5*a4(1,i,j,1)
         enddo
      enddo
      do k=2,km-1
         do j=j1,j2
            do i=i1, i2
               grat = delp(i,j,k-1) / delp(i,j,k)
               bet =  2. + grat + grat - gam(i,j,k)
               q(i,j,k) = (3.*(a4(1,i,j,k-1)+a4(1,i,j,k)) - q(i,j,k-1))/bet
               gam(i,j,k+1) = grat / bet
            enddo
         enddo
      enddo
      do j=j1,j2
         do i=i1,i2
            grat = delp(i,j,km-1) / delp(i,j,km) 
            q(i,j,km) = (3.*(a4(1,i,j,km-1)+a4(1,i,j,km)) - grat*qs(i,j) - q(i,j,km-1)) /  &
                 (2. + grat + grat - gam(i,j,km))
            q(i,j,km+1) = qs(i,j)
         enddo
      enddo
      do k=km-1,1,-1
         do j=j1,j2
            do i=i1,i2
               q(i,j,k) = q(i,j,k) - gam(i,j,k+1)*q(i,j,k+1)
            enddo
         enddo
      enddo
   else
      do j=j1,j2
         do i=i1,i2
            grat = delp(i,j,2) / delp(i,j,1)   ! grid ratio
            bet = grat*(grat+0.5)
            q(i,j,1) = ( (grat+grat)*(grat+1.)*a4(1,i,j,1) + a4(1,i,j,2) ) / bet
            gam(i,j,1) = ( 1. + grat*(grat+1.5) ) / bet
         enddo
      enddo
      do k=2,km
         do j=j1,j2
            do i=i1,i2
               d4(i,j) = delp(i,j,k-1) / delp(i,j,k)
               bet =  2. + d4(i,j) + d4(i,j) - gam(i,j,k-1)
               q(i,j,k) = ( 3.*(a4(1,i,j,k-1)+d4(i,j)*a4(1,i,j,k)) - q(i,j,k-1) )/bet
               gam(i,j,k) = d4(i,j) / bet
            enddo
         enddo
      enddo
      do j=j1,j2
         do i=i1,i2
            a_bot = 1. + d4(i,j)*(d4(i,j)+1.5)
            q(i,j,km+1) = (2.*d4(i,j)*(d4(i,j)+1.)*a4(1,i,j,km)+a4(1,i,j,km-1)-a_bot*q(i,j,km))  &
                 / ( d4(i,j)*(d4(i,j)+0.5) - a_bot*gam(i,j,km) )
         enddo
      enddo
      do k=km,1,-1
         do j=j1,j2
            do i=i1,i2
               q(i,j,k) = q(i,j,k) - gam(i,j,k)*q(i,j,k+1)
            enddo
         enddo
      enddo
   endif
!$ACC end kernels
   !----- Perfectly linear scheme --------------------------------
   if ( abs(kord) > 16 ) then
!$ACC kernels
      do k=1,km
         do j=j1,j2
            do i=i1,i2
               a4(2,i,j,k) = q(i,j,k  )
               a4(3,i,j,k) = q(i,j,k+1)
               a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
            enddo
         enddo
      enddo
!$ACC end kernels
      return
   endif
!$ACC kernels

   !----- Perfectly linear scheme --------------------------------
   !------------------
   ! Apply constraints
   !------------------
   im = i2 - i1 + 1

   ! Apply *large-scale* constraints
   do j=j1,j2 
      do i=i1,i2
         q(i,j,2) = min( q(i,j,2), max(a4(1,i,j,1), a4(1,i,j,2)) )
         q(i,j,2) = max( q(i,j,2), min(a4(1,i,j,1), a4(1,i,j,2)) )
      enddo
   enddo

   do k=2,km
      do j=j1,j2
         do i=i1,i2
            gam(i,j,k) = a4(1,i,j,k) - a4(1,i,j,k-1)
         enddo
      enddo
   enddo

      ! Interior:
   do k=3,km-1
      do j=j1,j2
         do i=i1,i2
            if ( gam(i,j,k-1)*gam(i,j,k+1)>0. ) then
               ! Apply large-scale constraint to ALL fields if not local max/min
               q(i,j,k) = min( q(i,j,k), max(a4(1,i,j,k-1),a4(1,i,j,k)) )
               q(i,j,k) = max( q(i,j,k), min(a4(1,i,j,k-1),a4(1,i,j,k)) )
            else
               if ( gam(i,j,k-1) > 0. ) then
                  ! There exists a local max
                  q(i,j,k) = max(q(i,j,k), min(a4(1,i,j,k-1),a4(1,i,j,k)))
               else
                  ! There exists a local min
                  q(i,j,k) = min(q(i,j,k), max(a4(1,i,j,k-1),a4(1,i,j,k)))
                  if ( iv==0 ) q(i,j,k) = max(0., q(i,j,k))
               endif
            endif
         enddo
      enddo
   enddo

      ! Bottom:
   do j=j1,j2
      do i=i1,i2
         q(i,j,km) = min( q(i,j,km), max(a4(1,i,j,km-1), a4(1,i,j,km)) )
         q(i,j,km) = max( q(i,j,km), min(a4(1,i,j,km-1), a4(1,i,j,km)) )
      enddo
   enddo
   do k=1,km
      do j=j1,j2
         do i=i1,i2
            a4(2,i,j,k) = q(i,j,k  )
            a4(3,i,j,k) = q(i,j,k+1)
         enddo
      enddo

      do j=j1,j2
         if ( k==1 .or. k==km ) then
            do i=i1,i2
               extm(i,j,k) = (a4(2,i,j,k)-a4(1,i,j,k)) * (a4(3,i,j,k)-a4(1,i,j,k)) > 0.
            enddo
         else
            do i=i1,i2
               extm(i,j,k) = gam(i,j,k)*gam(i,j,k+1) < 0.
            enddo
         endif
         if ( abs(kord)==16 ) then
            do i=i1,i2
               a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
               ext6(i,j,k) = abs(a4(4,i,j,k)) > abs(a4(2,i,j,k)-a4(3,i,j,k))
            enddo
         endif
      enddo
   enddo
   !---------------------------
   ! Apply subgrid constraints:
   !---------------------------
   ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
   ! Top 2 and bottom 2 layers always use monotonic mapping

   if ( iv==0 ) then
      do j=j1,j2
         do i=i1,i2
            a4(2,i,j,1) = max(0., a4(2,i,j,1))
         enddo
      enddo
   elseif ( iv==-1 ) then
      do j=j1,j2 
         do i=i1,i2
            if ( a4(2,i,j,1)*a4(1,i,j,1) <= 0. ) a4(2,i,j,1) = 0.
         enddo
      enddo
   elseif ( iv==2 ) then
      do j=j1,j2
         do i=i1,i2
            a4(2,i,j,1) = a4(1,i,j,1)
            a4(3,i,j,1) = a4(1,i,j,1)
            a4(4,i,j,1) = 0.
         enddo
      enddo
   endif
!$ACC end kernels
   if ( iv/=2 ) then
!$ACC kernels
      do j=j1,j2
         do i=i1,i2
            a4(4,i,j,1) = 3.*(2.*a4(1,i,j,1) - (a4(2,i,j,1)+a4(3,i,j,1)))
         enddo
      enddo
!$ACC end kernels
      call cs_limiters(i1,i2,j1,j2,1,1,extm(i1,j1,1), a4(1,i1,j1,1), 1)
   endif
   ! k=2
!$ACC kernels
   do j=j1,j2
      do i=i1,i2
         a4(4,i,j,2) = 3.*(2.*a4(1,i,j,2) - (a4(2,i,j,2)+a4(3,i,j,2)))
      enddo
   enddo
!$ACC end kernels
   call cs_limiters(i1,i2,j1,j2,2,2,extm(i1,j1,2), a4(1,i1,j1,2), 2)
   !-------------------------------------
   ! Huynh's 2nd constraint for interior:
   !-------------------------------------
!$ACC kernels
   if ( abs(kord)<9 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               ! Left  edges
               pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
               lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
               a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),   &
                    max(a4(1,i,j,k), pmp_1, lac_1) )
               ! Right edges
               pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
               lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
               a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),    &
                    max(a4(1,i,j,k), pmp_2, lac_2) )

               a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
            enddo
         enddo
      enddo
   elseif ( abs(kord)==9 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if ( extm(i,j,k) .and. extm(i,j,k-1) ) then
                  ! grid-scale 2-delta-z wave detected
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else if ( extm(i,j,k) .and. extm(i,j,k+1) ) then
                  ! grid-scale 2-delta-z wave detected
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else if ( extm(i,j,k) .and. a4(1,i,j,k)<qmin ) then
                  ! grid-scale 2-delta-z wave detected
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else
                  a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
                  ! Check within the smooth region if subgrid profile is non-monotonic
                  if( abs(a4(4,i,j,k)) > abs(a4(2,i,j,k)-a4(3,i,j,k)) ) then
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),  &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),  &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
                  endif
               endif
            enddo
         enddo
      enddo
   elseif ( abs(kord)==10 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if( extm(i,j,k) ) then
                  if( a4(1,i,j,k)<qmin .or. extm(i,j,k-1) .or. extm(i,j,k+1) ) then
                     ! grid-scale 2-delta-z wave detected; or q is too small -> ehance vertical mixing
                     a4(2,i,j,k) = a4(1,i,j,k)
                     a4(3,i,j,k) = a4(1,i,j,k)
                     a4(4,i,j,k) = 0.
                  else
                     ! True local extremum
                     a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  endif
               else        ! not a local extremum
                  a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  ! Check within the smooth region if subgrid profile is non-monotonic
                  if( abs(a4(4,i,j,k)) > abs(a4(2,i,j,k)-a4(3,i,j,k)) ) then
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),  &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),  &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  endif
               endif
            enddo
         enddo
      enddo
   elseif ( abs(kord)==12 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if( extm(i,j,k) ) then
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else        ! not a local extremum
                  a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  ! Check within the smooth region if subgrid profile is non-monotonic
                  if( abs(a4(4,i,j,k)) > abs(a4(2,i,j,k)-a4(3,i,j,k)) ) then
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),  &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),  &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  endif
               endif
            enddo
         enddo
      enddo
   elseif ( abs(kord)==13 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if( extm(i,j,k) ) then
                  if ( extm(i,j,k-1) .and. extm(i,j,k+1) ) then
                     ! grid-scale 2-delta-z wave detected
                     a4(2,i,j,k) = a4(1,i,j,k)
                     a4(3,i,j,k) = a4(1,i,j,k)
                     a4(4,i,j,k) = 0.
                  else
                     ! Left  edges
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),   &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     ! Right edges
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),    &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
                  endif
               else
                  a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
               endif
            enddo
         enddo
      enddo
   elseif ( abs(kord)==14 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
            enddo
         enddo
      enddo
   elseif ( abs(kord)==16 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if( ext6(i,j,k) ) then
                  if ( extm(i,j,k-1) .or. extm(i,j,k+1) ) then
                     ! Left  edges
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),   &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     ! Right edges
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),    &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
                  endif
               endif
            enddo
         enddo
      enddo
   else      ! kord = 11, 13
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if ( extm(i,j,k) .and. (extm(i,j,k-1).or.extm(i,j,k+1).or.a4(1,i,j,k)<qmin) ) then
                  ! Noisy region:
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else
                  a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
               endif
            enddo
         enddo
      enddo
   endif
!$acc end kernels
   ! Additional constraint to ensure positivity
   if ( iv==0 ) call cs_limiters(i1,i2,j1,j2,3,km-2, extm(i1,j1,3), a4(1,i1,j1,3), 0)
   !----------------------------------
   ! Bottom layer subgrid constraints:
   !----------------------------------
!$ACC kernels
   if ( iv==0 ) then
      do j=j1,j2
         do i=i1,i2
            a4(3,i,j,km) = max(0., a4(3,i,j,km))
         enddo
      enddo
   elseif ( iv .eq. -1 ) then 
      do j=j1,j2
         do i=i1,i2
            if ( a4(3,i,j,km)*a4(1,i,j,km) <= 0. )  a4(3,i,j,km) = 0.
         enddo
      enddo
   endif

   do k=km-1,km
      do j=j1,j2
         do i=i1,i2
            a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
         enddo
      enddo
   enddo
!$ACC end kernels
   call cs_limiters(i1,i2,j1,j2,km-1,km-1, extm(i1,j1,km-1), a4(1,i1,j1,km-1), 2)
   call cs_limiters(i1,i2,j1,j2,km,km, extm(i1,j1,km), a4(1,i1,j1,km), 1)

 end subroutine scalar_profile


 subroutine cs_profile(qs, a4, delp, km, i1, i2, j1, j2, iv, kord)
   ! Optimized vertical profile reconstruction:
   ! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
   integer, intent(in):: i1, i2, j1, j2
   integer, intent(in):: km      ! vertical dimension
   integer, intent(in):: iv      ! iv =-1: winds
   ! iv = 0: positive definite scalars
   ! iv = 1: others
   integer, intent(in):: kord
   real, intent(in)   ::   qs(i1:,j1:)
   real, intent(in)   :: delp(i1:i2,j1:j2,km)     ! layer pressure thickness
   real, intent(inout):: a4(4,i1:i2,j1:j2,km)     ! Interpolated values
   !-----------------------------------------------------------------------
   logical:: extm(i1:i2,j1:j2,km) 
   real  gam(i1:i2,j1:j2,km)
   real    q(i1:i2,j1:j2,km+1)
   real   d4(i1:i2,j1:j2)
   real   bet, a_bot, grat 
   real   pmp_1, lac_1, pmp_2, lac_2
   integer i, j, k, im

   if ( iv .eq. -2 ) then
      do j=j1,j2
         do i=i1,i2
            gam(i,j,2) = 0.5
            q(i,j,1) = 1.5*a4(1,i,j,1)
         enddo
      enddo
      do k=2,km-1
         do j=j1,j2
            do i=i1, i2
               grat = delp(i,j,k-1) / delp(i,j,k)
               bet =  2. + grat + grat - gam(i,j,k)
               q(i,j,k) = (3.*(a4(1,i,j,k-1)+a4(1,i,j,k)) - q(i,j,k-1))/bet
               gam(i,j,k+1) = grat / bet
            enddo
         enddo
      enddo
      do j=j1,j2
         do i=i1,i2
            grat = delp(i,j,km-1) / delp(i,j,km) 
            q(i,j,km) = (3.*(a4(1,i,j,km-1)+a4(1,i,j,km)) - grat*qs(i,j) - q(i,j,km-1)) /  &
                 (2. + grat + grat - gam(i,j,km))
            q(i,j,km+1) = qs(i,j)
         enddo
      enddo
      do k=km-1,1,-1
         do j=j1,j2
            do i=i1,i2
               q(i,j,k) = q(i,j,k) - gam(i,j,k+1)*q(i,j,k+1)
            enddo
         enddo
      enddo
   else
      do j=j1,j2
         do i=i1,i2
            grat = delp(i,j,2) / delp(i,j,1)   ! grid ratio
            bet = grat*(grat+0.5)
            q(i,j,1) = ( (grat+grat)*(grat+1.)*a4(1,i,j,1) + a4(1,i,j,2) ) / bet
            gam(i,j,1) = ( 1. + grat*(grat+1.5) ) / bet
         enddo
      enddo
      do k=2,km
         do j=j1,j2
            do i=i1,i2
               d4(i,j) = delp(i,j,k-1) / delp(i,j,k)
               bet =  2. + d4(i,j) + d4(i,j) - gam(i,j,k-1)
               q(i,j,k) = ( 3.*(a4(1,i,j,k-1)+d4(i,j)*a4(1,i,j,k)) - q(i,j,k-1) )/bet
               gam(i,j,k) = d4(i,j) / bet
            enddo
         enddo
      enddo
      do j=j1,j2
         do i=i1,i2
            a_bot = 1. + d4(i,j)*(d4(i,j)+1.5)
            q(i,j,km+1) = (2.*d4(i,j)*(d4(i,j)+1.)*a4(1,i,j,km)+a4(1,i,j,km-1)-a_bot*q(i,j,km))  &
                 / ( d4(i,j)*(d4(i,j)+0.5) - a_bot*gam(i,j,km) )
         enddo
      enddo
      do k=km,1,-1
         do j=j1,j2
            do i=i1,i2
               q(i,j,k) = q(i,j,k) - gam(i,j,k)*q(i,j,k+1)
            enddo
         enddo
      enddo
   endif
   !----- Perfectly linear scheme --------------------------------
   if ( abs(kord) > 16 ) then
      do k=1,km
         do j=j1,j2
            do i=i1,i2
               a4(2,i,j,k) = q(i,j,k  )
               a4(3,i,j,k) = q(i,j,k+1)
               a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
            enddo
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
   do j=j1,j2
      do i=i1,i2
         q(i,j,2) = min( q(i,j,2), max(a4(1,i,j,1), a4(1,i,j,2)) )
         q(i,j,2) = max( q(i,j,2), min(a4(1,i,j,1), a4(1,i,j,2)) )
      enddo
   enddo

   do k=2,km
      do j=j1,j2
         do i=i1,i2
            gam(i,j,k) = a4(1,i,j,k) - a4(1,i,j,k-1)
         enddo
      enddo
   enddo

      ! Interior:
   do k=3,km-1
      do j=j1,j2
         do i=i1,i2
            if ( gam(i,j,k-1)*gam(i,j,k+1)>0. ) then
               ! Apply large-scale constraint to ALL fields if not local max/min
               q(i,j,k) = min( q(i,j,k), max(a4(1,i,j,k-1),a4(1,i,j,k)) )
               q(i,j,k) = max( q(i,j,k), min(a4(1,i,j,k-1),a4(1,i,j,k)) )
            else
               if ( gam(i,j,k-1) > 0. ) then
                  ! There exists a local max
                  q(i,j,k) = max(q(i,j,k), min(a4(1,i,j,k-1),a4(1,i,j,k)))
               else
                  ! There exists a local min
                  q(i,j,k) = min(q(i,j,k), max(a4(1,i,j,k-1),a4(1,i,j,k)))
                  if ( iv==0 ) q(i,j,k) = max(0., q(i,j,k))
               endif
            endif
         enddo
      enddo
   enddo
      ! Bottom:
   do j=j1,j2
      do i=i1,i2
         q(i,j,km) = min( q(i,j,km), max(a4(1,i,j,km-1), a4(1,i,j,km)) )
         q(i,j,km) = max( q(i,j,km), min(a4(1,i,j,km-1), a4(1,i,j,km)) )
      enddo
   enddo

   do k=1,km
      do j=j1,j2
         do i=i1,i2
            a4(2,i,j,k) = q(i,j,k  )
            a4(3,i,j,k) = q(i,j,k+1)
         enddo
      enddo

      if ( k==1 .or. k==km ) then
         do j=j1,j2
            do i=i1,i2
               extm(i,j,k) = (a4(2,i,j,k)-a4(1,i,j,k)) * (a4(3,i,j,k)-a4(1,i,j,k)) > 0.
            enddo
         enddo
      else
         do j=j1,j2
            do i=i1,i2
               extm(i,j,k) = gam(i,j,k)*gam(i,j,k+1) < 0.
            enddo
         enddo
      endif
   enddo
   !---------------------------
   ! Apply subgrid constraints:
   !---------------------------
   ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
   ! Top 2 and bottom 2 layers always use monotonic mapping
   if ( iv==0 ) then
      do j=j1,j2
         do i=i1,i2
            a4(2,i,j,1) = max(0., a4(2,i,j,1))
         enddo
      enddo
   elseif ( iv==-1 ) then 
      do j=j1,j2
         do i=i1,i2
            if ( a4(2,i,j,1)*a4(1,i,j,1) <= 0. ) a4(2,i,j,1) = 0.
         enddo
      enddo
   elseif ( iv==2 ) then
      do j=j1,j2
         do i=i1,i2
            a4(2,i,j,1) = a4(1,i,j,1)
            a4(3,i,j,1) = a4(1,i,j,1)
            a4(4,i,j,1) = 0.
         enddo
      enddo
   endif

   if ( iv/=2 ) then
      do j=j1,j2
         do i=i1,i2
            a4(4,i,j,1) = 3.*(2.*a4(1,i,j,1) - (a4(2,i,j,1)+a4(3,i,j,1)))
            a4(4,i,j,2) = 3.*(2.*a4(1,i,j,2) - (a4(2,i,j,2)+a4(3,i,j,2)))
         enddo
      enddo
      call cs_limiters(i1,i2,j1,j2,1,1, extm(i1,j1,1), a4(1,i1,j1,1), 1)
   else
      do j=j1,j2
         do i=i1,i2
            a4(4,i,j,2) = 3.*(2.*a4(1,i,j,2) - (a4(2,i,j,2)+a4(3,i,j,2)))
         enddo
      enddo
   endif
   call cs_limiters(i1,i2,j1,j2,2,2, extm(i1,j1,2), a4(1,i1,j1,2), 2)

   !-------------------------------------
   ! Huynh's 2nd constraint for interior:
   !-------------------------------------
   if ( abs(kord)<9 ) then
      do k=3,km-2
         do j=j1,j2      
            do i=i1,i2
               ! Left  edges
               pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
               lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
               a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),   &
                    max(a4(1,i,j,k), pmp_1, lac_1) )
               ! Right edges
               pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
               lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
               a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),    &
                    max(a4(1,i,j,k), pmp_2, lac_2) )

               a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
            enddo
         enddo
      enddo
   elseif ( abs(kord)==9 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if ( extm(i,j,k) .and. extm(i,j,k-1) ) then  ! c90_mp122
                  ! grid-scale 2-delta-z wave detected
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else if ( extm(i,j,k) .and. extm(i,j,k+1) ) then  ! c90_mp122
                  ! grid-scale 2-delta-z wave detected
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else
                  a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  ! Check within the smooth region if subgrid profile is non-monotonic
                  if( abs(a4(4,i,j,k)) > abs(a4(2,i,j,k)-a4(3,i,j,k)) ) then
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),  &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),  &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  endif
               endif
            enddo
         enddo
      enddo
   elseif ( abs(kord)==10 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if( extm(i,j,k) ) then
                  if( extm(i,j,k-1) .or. extm(i,j,k+1) ) then
                     ! grid-scale 2-delta-z wave detected
                     a4(2,i,j,k) = a4(1,i,j,k)
                     a4(3,i,j,k) = a4(1,i,j,k)
                     a4(4,i,j,k) = 0.
                  else
                     ! True local extremum
                     a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  endif
               else        ! not a local extremum
                  a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  ! Check within the smooth region if subgrid profile is non-monotonic
                  if( abs(a4(4,i,j,k)) > abs(a4(2,i,j,k)-a4(3,i,j,k)) ) then
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),  &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),  &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  endif
               endif
            enddo
         enddo
      enddo
   elseif ( abs(kord)==12 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if( extm(i,j,k) ) then
                  ! grid-scale 2-delta-z wave detected
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else        ! not a local extremum
                  a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  ! Check within the smooth region if subgrid profile is non-monotonic
                  if( abs(a4(4,i,j,k)) > abs(a4(2,i,j,k)-a4(3,i,j,k)) ) then
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),  &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),  &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 6.*a4(1,i,j,k) - 3.*(a4(2,i,j,k)+a4(3,i,j,k))
                  endif
               endif
            enddo
         enddo
      enddo
   elseif ( abs(kord)==13 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if( extm(i,j,k) ) then
                  if ( extm(i,j,k-1) .and. extm(i,j,k+1) ) then
                     ! grid-scale 2-delta-z wave detected
                     a4(2,i,j,k) = a4(1,i,j,k)
                     a4(3,i,j,k) = a4(1,i,j,k)
                     a4(4,i,j,k) = 0.
                  else
                     ! Left  edges
                     pmp_1 = a4(1,i,j,k) - 2.*gam(i,j,k+1)
                     lac_1 = pmp_1 + 1.5*gam(i,j,k+2)
                     a4(2,i,j,k) = min(max(a4(2,i,j,k), min(a4(1,i,j,k), pmp_1, lac_1)),   &
                          max(a4(1,i,j,k), pmp_1, lac_1) )
                     ! Right edges
                     pmp_2 = a4(1,i,j,k) + 2.*gam(i,j,k)
                     lac_2 = pmp_2 - 1.5*gam(i,j,k-1)
                     a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), pmp_2, lac_2)),    &
                          max(a4(1,i,j,k), pmp_2, lac_2) )
                     a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
                  endif
               else
                  a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
               endif
            enddo
         enddo
      enddo
   elseif ( abs(kord)==14 ) then
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
            enddo
         enddo
      enddo
   else      ! kord = 11
      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               if ( extm(i,j,k) .and. (extm(i,j,k-1) .or. extm(i,j,k+1)) ) then
                  ! Noisy region:
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else
                  a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
               endif
            enddo
         enddo
      enddo
   endif

   ! Additional constraint to ensure positivity
   if ( iv==0 ) call cs_limiters(i1,i2,j1,j2,3,km-2, extm(i1,j1,3), a4(1,i1,j1,3), 0)
   !----------------------------------
   ! Bottom layer subgrid constraints:
   !----------------------------------
   if ( iv==0 ) then
      do j=j1,j2
         do i=i1,i2
            a4(3,i,j,km) = max(0., a4(3,i,j,km))
         enddo
      enddo
   elseif ( iv .eq. -1 ) then 
      do j=j1,j2
         do i=i1,i2
            if ( a4(3,i,j,km)*a4(1,i,j,km) <= 0. )  a4(3,i,j,km) = 0.
         enddo
      enddo
   endif

   do k=km-1,km
       do j=j1,j2
         do i=i1,i2
            a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
         enddo
      enddo
   enddo
   call cs_limiters(i1,i2,j1,j2,km-1,km-1, extm(i1,j1,km-1), a4(1,i1,j1,km-1), 2)
   call cs_limiters(i1,i2,j1,j2,km,km, extm(i1,j1,km), a4(1,i1,j1,km), 1)

 end subroutine cs_profile


 subroutine cs_limiters(i1, i2, j1, j2, k1, k2, extm, a4, iv)
   integer, intent(in) :: i1, i2, j1, j2, k1, k2
   integer, intent(in) :: iv
   logical, intent(in) :: extm(i1:i2,j1:j2,k1:k2)
   real , intent(inout) :: a4(4,i1:i2,j1:j2,k1:k2)   ! PPM array
   ! !LOCAL VARIABLES:
   real  da1, da2, a6da
   integer i,j,k
!$acc kernels
   if ( iv==0 ) then
      ! Positive definite constraint
      do j=j1,j2
         do k=k1,k2
            do i=i1,i2
               if( a4(1,i,j,k)<=0.) then
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else
                  if( abs(a4(3,i,j,k)-a4(2,i,j,k)) < -a4(4,i,j,k) ) then
                     if( (a4(1,i,j,k)+0.25*(a4(3,i,j,k)-a4(2,i,j,k))**2/a4(4,i,j,k)+a4(4,i,j,k)*r12) < 0. ) then
                        ! local minimum is negative
                        if( a4(1,i,j,k)<a4(3,i,j,k) .and. a4(1,i,j,k)<a4(2,i,j,k) ) then
                           a4(3,i,j,k) = a4(1,i,j,k)
                           a4(2,i,j,k) = a4(1,i,j,k)
                           a4(4,i,j,k) = 0.
                        elseif( a4(3,i,j,k) > a4(2,i,j,k) ) then
                           a4(4,i,j,k) = 3.*(a4(2,i,j,k)-a4(1,i,j,k))
                           a4(3,i,j,k) = a4(2,i,j,k) - a4(4,i,j,k)
                        else
                           a4(4,i,j,k) = 3.*(a4(3,i,j,k)-a4(1,i,j,k))
                           a4(2,i,j,k) = a4(3,i,j,k) - a4(4,i,j,k)
                        endif
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
   elseif ( iv==1 ) then
      do j=j1,j2
         do k=k1,k2
            do i=i1,i2
               if( (a4(1,i,j,k)-a4(2,i,j,k))*(a4(1,i,j,k)-a4(3,i,j,k))>=0. ) then
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else
                  da1  = a4(3,i,j,k) - a4(2,i,j,k)
                  da2  = da1**2
                  a6da = a4(4,i,j,k)*da1
                  if(a6da < -da2) then
                     a4(4,i,j,k) = 3.*(a4(2,i,j,k)-a4(1,i,j,k))
                     a4(3,i,j,k) = a4(2,i,j,k) - a4(4,i,j,k)
                  elseif(a6da > da2) then
                     a4(4,i,j,k) = 3.*(a4(3,i,j,k)-a4(1,i,j,k))
                     a4(2,i,j,k) = a4(3,i,j,k) - a4(4,i,j,k)
                  endif
               endif
            enddo
         enddo
      enddo
   else
      ! Standard PPM constraint
      do j=j1,j2
         do k=k1,k2
            do i=i1,i2
               if( extm(i,j,k) ) then
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else
                  da1  = a4(3,i,j,k) - a4(2,i,j,k)
                  da2  = da1**2
                  a6da = a4(4,i,j,k)*da1
                  if(a6da < -da2) then
                     a4(4,i,j,k) = 3.*(a4(2,i,j,k)-a4(1,i,j,k))
                     a4(3,i,j,k) = a4(2,i,j,k) - a4(4,i,j,k)
                  elseif(a6da > da2) then
                     a4(4,i,j,k) = 3.*(a4(3,i,j,k)-a4(1,i,j,k))
                     a4(2,i,j,k) = a4(3,i,j,k) - a4(4,i,j,k)
                  endif
               endif
            enddo
         enddo
      enddo
   endif
!$ACC end kernels
 end subroutine cs_limiters



 subroutine ppm_profile(a4, delp, km, i1, i2, j1, j2, iv, kord)
!$ACC routine seq
   ! !INPUT PARAMETERS:
   integer, intent(in):: iv      ! iv =-1: winds
   ! iv = 0: positive definite scalars
   ! iv = 1: others
   ! iv = 2: temp (if remap_t) and w (iv=-2)
   integer, intent(in):: i1      ! Starting longitude
   integer, intent(in):: i2      ! Finishing longitude
   integer, intent(in):: j1      ! Starting latitude
   integer, intent(in):: j2      ! Finishing latitude
   integer, intent(in):: km      ! vertical dimension
   integer, intent(in):: kord    ! Order (or more accurately method no.):
   ! 
   real , intent(in):: delp(i1:i2,j1:j2,km)     ! layer pressure thickness

   ! !INPUT/OUTPUT PARAMETERS:
   real , intent(inout):: a4(4,i1:i2,j1:j2,km)  ! Interpolated values

   ! DESCRIPTION:
   !
   !   Perform the piecewise parabolic reconstruction
   ! 
   ! !REVISION HISTORY: 
   ! S.-J. Lin   revised at GFDL 2007
   !-----------------------------------------------------------------------
   ! local arrays:
   real    dc(i1:i2,j1:j2,km)
   real    h2(i1:i2,j1:j2,km)
   real  delq(i1:i2,j1:j2,km)
   real   df2(i1:i2,j1:j2,km)
   real    d4(i1:i2,j1:j2,km)

   ! local scalars:
   integer i, j, k, km1, lmt, it
   real  fac
   real  a1, a2, c1, c2, c3, d1, d2
   real  qm, dq, lac, qmp, pmp

   km1 = km - 1
   it = i2 - i1 + 1
   do k=2,km
      do j=j1,j2
         do i=i1,i2
            delq(i,j,k-1) =   a4(1,i,j,k) - a4(1,i,j,k-1)
            d4(i,j,k  ) = delp(i,j,k-1) + delp(i,j,k)
         enddo
      enddo
   enddo

   do k=2,km1
      do j=j1,j2
         do i=i1,i2
            c1  = (delp(i,j,k-1)+0.5*delp(i,j,k))/d4(i,j,k+1)
            c2  = (delp(i,j,k+1)+0.5*delp(i,j,k))/d4(i,j,k)
            df2(i,j,k) = delp(i,j,k)*(c1*delq(i,j,k) + c2*delq(i,j,k-1)) /      &
                 (d4(i,j,k)+delp(i,j,k+1))
            dc(i,j,k) = sign( min(abs(df2(i,j,k)),              &
                 max(a4(1,i,j,k-1),a4(1,i,j,k),a4(1,i,j,k+1))-a4(1,i,j,k),  &
                 a4(1,i,j,k)-min(a4(1,i,j,k-1),a4(1,i,j,k),a4(1,i,j,k+1))), df2(i,j,k) )
         enddo
      enddo
   enddo
      !-----------------------------------------------------------
      ! 4th order interpolation of the provisional cell edge value
      !-----------------------------------------------------------

   do k=3,km1
      do j=j1,j2
         do i=i1,i2
            c1 = delq(i,j,k-1)*delp(i,j,k-1) / d4(i,j,k)
            a1 = d4(i,j,k-1) / (d4(i,j,k) + delp(i,j,k-1))
            a2 = d4(i,j,k+1) / (d4(i,j,k) + delp(i,j,k))
            a4(2,i,j,k) = a4(1,i,j,k-1) + c1 + 2./(d4(i,j,k-1)+d4(i,j,k+1)) *    &
                 ( delp(i,j,k)*(c1*(a1 - a2)+a2*dc(i,j,k-1)) -          &
                 delp(i,j,k-1)*a1*dc(i,j,k  ) )
         enddo
      enddo
   enddo

      ! Area preserving cubic with 2nd deriv. = 0 at the boundaries
      ! Top
   do j=j1,j2
      do i=i1,i2
         d1 = delp(i,j,1)
         d2 = delp(i,j,2)
         qm = (d2*a4(1,i,j,1)+d1*a4(1,i,j,2)) / (d1+d2)
         dq = 2.*(a4(1,i,j,2)-a4(1,i,j,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,j,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,j,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
         ! Top edge:
         !-------------------------------------------------------
         a4(2,i,j,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,j,2)
         !-------------------------------------------------------
         !        a4(2,i,j,1) = (12./7.)*a4(1,i,j,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
         !-------------------------------------------------------
         ! No over- and undershoot condition
         a4(2,i,j,2) = max( a4(2,i,j,2), min(a4(1,i,j,1), a4(1,i,j,2)) )
         a4(2,i,j,2) = min( a4(2,i,j,2), max(a4(1,i,j,1), a4(1,i,j,2)) )
         dc(i,j,1) =  0.5*(a4(2,i,j,2) - a4(1,i,j,1))
      enddo
   enddo
      ! Enforce monotonicity  within the top layer

   if( iv==0 ) then
      do j=j1,j2
         do i=i1,i2
            a4(2,i,j,1) = max(0., a4(2,i,j,1))
            a4(2,i,j,2) = max(0., a4(2,i,j,2))
         enddo
      enddo
   elseif( iv==-1 ) then
      do j=j1,j2
         do i=i1,i2
            if ( a4(2,i,j,1)*a4(1,i,j,1) <= 0. ) a4(2,i,j,1) = 0.
         enddo
      enddo
   elseif( abs(iv)==2 ) then
      do j=j1,j2
         do i=i1,i2
            a4(2,i,j,1) = a4(1,i,j,1)
            a4(3,i,j,1) = a4(1,i,j,1)
         enddo
      enddo
   endif

      ! Bottom
      ! Area preserving cubic with 2nd deriv. = 0 at the surface
   do j=j1,j2
      do i=i1,i2
         d1 = delp(i,j,km)
         d2 = delp(i,j,km1)
         qm = (d2*a4(1,i,j,km)+d1*a4(1,i,j,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,j,km1)-a4(1,i,j,km)) / (d1+d2)
         c1 = (a4(2,i,j,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,j,km) = qm - c1*d1*d2*(d2+3.*d1)
         ! Bottom edge:
         !-----------------------------------------------------
         a4(3,i,j,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,j,km)
         !        dc(i,j,km) = 0.5*(a4(3,i,j,km) - a4(1,i,j,km))
         !-----------------------------------------------------
         !        a4(3,i,j,km) = (12./7.)*a4(1,i,j,km)-(13./14.)*a4(1,i,j,km-1)+(3./14.)*a4(1,i,j,km-2)
         ! No over- and under-shoot condition
         a4(2,i,j,km) = max( a4(2,i,j,km), min(a4(1,i,j,km), a4(1,i,j,km1)) )
         a4(2,i,j,km) = min( a4(2,i,j,km), max(a4(1,i,j,km), a4(1,i,j,km1)) )
         dc(i,j,km) = 0.5*(a4(1,i,j,km) - a4(2,i,j,km))
      enddo
   enddo

      ! Enforce constraint on the "slope" at the surface

#ifdef BOT_MONO
   do j=j1,j2
      do i=i1,i2
         a4(4,i,j,km) = 0
         if( a4(3,i,j,km) * a4(1,i,j,km) <= 0. ) a4(3,i,j,km) = 0.
         d1 = a4(1,i,j,km) - a4(2,i,j,km)
         d2 = a4(3,i,j,km) - a4(1,i,j,km)
         if ( d1*d2 < 0. ) then
            a4(2,i,j,km) = a4(1,i,j,km)
            a4(3,i,j,km) = a4(1,i,j,km)
         else
            dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,j,km-1))), d1)
            a4(2,i,j,km) = a4(1,i,j,km) - dq
            a4(3,i,j,km) = a4(1,i,j,km) + dq
         endif
      enddo
   enddo
#else
   if( iv==0 ) then
       do j=j1,j2
         do i=i1,i2
            a4(2,i,j,km) = max(0.,a4(2,i,j,km))
            a4(3,i,j,km) = max(0.,a4(3,i,j,km))
         enddo
      enddo
   elseif( iv<0 ) then
      do j=j1,j2
         do i=i1,i2
            if( a4(1,i,j,km)*a4(3,i,j,km) <= 0. )  a4(3,i,j,km) = 0.
         enddo
      enddo
   endif
#endif

   do k=1,km1
      do j=j1,j2
         do i=i1,i2
            a4(3,i,j,k) = a4(2,i,j,k+1)
         enddo
      enddo
   enddo
   !-----------------------------------------------------------
   ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
   !-----------------------------------------------------------
   ! Top 2 and bottom 2 layers always use monotonic mapping
   do k=1,2
      do j=j1,j2
         do i=i1,i2
            a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
         enddo
      enddo
   enddo
   call ppm_limiters(i1,i2,j1,j2,1,2,dc(i1,j1,1), a4(1,i1,j1,1), 0)

   if(kord >= 7) then
      !-----------------------
      ! Huynh's 2nd constraint
      !-----------------------
      do k=2,km1
         do j=j1,j2
            do i=i1,i2
               ! Method#1
               !           h2(i,j,k) = delq(i,j,k) - delq(i,j,k-1)
               ! Method#2 - better
               h2(i,j,k) = 2.*(dc(i,j,k+1)/delp(i,j,k+1) - dc(i,j,k-1)/delp(i,j,k-1))  &
                    / ( delp(i,j,k)+0.5*(delp(i,j,k-1)+delp(i,j,k+1)) )        &
                    * delp(i,j,k)**2 
               ! Method#3
!!!            h2(i,j,k) = dc(i,j,k+1) - dc(i,j,k-1)
            enddo
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
         do j=j1,j2
            do i=i1,i2
               ! Right edges
               !        qmp   = a4(1,i,j,k) + 2.0*delq(i,j,k-1)
               !        lac   = a4(1,i,j,k) + fac*h2(i,j,k-1) + 0.5*delq(i,j,k-1)
               !
               pmp   = 2.*dc(i,j,k)
               qmp   = a4(1,i,j,k) + pmp
               lac   = a4(1,i,j,k) + fac*h2(i,j,k-1) + dc(i,j,k)
               a4(3,i,j,k) = min(max(a4(3,i,j,k), min(a4(1,i,j,k), qmp, lac)),    &
                    max(a4(1,i,j,k), qmp, lac) )
               ! Left  edges
               !        qmp   = a4(1,i,j,k) - 2.0*delq(i,j,k)
               !        lac   = a4(1,i,j,k) + fac*h2(i,j,k+1) - 0.5*delq(i,j,k)
               !
               qmp   = a4(1,i,j,k) - pmp
               lac   = a4(1,i,j,k) + fac*h2(i,j,k+1) - dc(i,j,k)
               a4(2,i,j,k) = min(max(a4(2,i,j,k),  min(a4(1,i,j,k), qmp, lac)),   &
                    max(a4(1,i,j,k), qmp, lac))
               !-------------
               ! Recompute A6
               !-------------
               a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
            enddo
         enddo
      enddo
      ! Additional constraint to ensure positivity when kord=7
      if (iv == 0 .and. kord >= 6 )  call ppm_limiters(i1,i2,j1,j2,3,km-2,dc(i1,j1,3), a4(1,i1,j1,3), 2)

   else

      lmt = kord - 3
      lmt = max(0, lmt)
      if (iv == 0) lmt = min(2, lmt)
      if( kord /= 4) then
         do k=3,km-2
            do j=j1,j2
               do i=i1,i2
                  a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
               enddo
            enddo
         enddo
      endif
      if(kord/=6) call ppm_limiters(i1,i2,j1,j2,3,km-2,dc(i1,j1,3), a4(1,i1,j1,3), lmt)
   endif

   do k=km1,km
      do j=j1,j2
         do i=i1,i2
            a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)))
         enddo
      enddo
   enddo

   call ppm_limiters(i1,i2,j1,j2,km1,km,dc(i1,j1,km1), a4(1,i1,j1,km1), 0)

 end subroutine ppm_profile


 subroutine ppm_limiters(i1,i2,j1,j2,k1,k2,dm, a4, lmt)
!$ACC routine seq
   ! !INPUT PARAMETERS:
   integer, intent(in) :: i1,i2,j1,j2,k1,k2
   real , intent(in):: dm(i1:i2,j1:j2,k1:k2)     ! the linear slope
   integer, intent(in) :: lmt       ! 0: Standard PPM constraint
   ! 1: Improved full monotonicity constraint (Lin)
   ! 2: Positive definite constraint
   ! 3: do nothing (return immediately)
   ! !INPUT/OUTPUT PARAMETERS:
   real , intent(inout) :: a4(4,i1:i2,j1:j2,k1:k2)   ! PPM array
   ! AA <-- a4(1,i)
   ! AL <-- a4(2,i)
   ! AR <-- a4(3,i)
   ! A6 <-- a4(4,i)
   ! !LOCAL VARIABLES:
   real  qmp
   real  da1, da2, a6da
   real  fmin
   integer i, j, k

   ! Developer: S.-J. Lin

   if ( lmt == 3 ) return

   if(lmt == 0) then
      ! Standard PPM constraint
      do k = k1,k2
         do j=j1,j2
            do i=i1,i2
               if(dm(i,j,k) == 0.) then
                  a4(2,i,j,k) = a4(1,i,j,k)
                  a4(3,i,j,k) = a4(1,i,j,k)
                  a4(4,i,j,k) = 0.
               else
                  da1  = a4(3,i,j,k) - a4(2,i,j,k)
                  da2  = da1**2
                  a6da = a4(4,i,j,k)*da1
                  if(a6da < -da2) then
                     a4(4,i,j,k) = 3.*(a4(2,i,j,k)-a4(1,i,j,k))
                     a4(3,i,j,k) = a4(2,i,j,k) - a4(4,i,j,k)
                  elseif(a6da > da2) then
                     a4(4,i,j,k) = 3.*(a4(3,i,j,k)-a4(1,i,j,k))
                     a4(2,i,j,k) = a4(3,i,j,k) - a4(4,i,j,k)
                  endif
               endif
            enddo
         enddo
      enddo

   elseif (lmt == 1) then

      ! Improved full monotonicity constraint (Lin 2004)
      ! Note: no need to provide first guess of A6 <-- a4(4,i,j,k)
      do k = k1,k2
         do j=j1,j2
            do i=i1,i2
               qmp = 2.*dm(i,j,k)
               a4(2,i,j,k) = a4(1,i,j,k)-sign(min(abs(qmp),abs(a4(2,i,j,k)-a4(1,i,j,k))), qmp)
               a4(3,i,j,k) = a4(1,i,j,k)+sign(min(abs(qmp),abs(a4(3,i,j,k)-a4(1,i,j,k))), qmp)
               a4(4,i,j,k) = 3.*( 2.*a4(1,i,j,k) - (a4(2,i,j,k)+a4(3,i,j,k)) )
            enddo
         enddo
      enddo
   elseif (lmt == 2) then

      ! Positive definite constraint
      do k = k1,k2
         do j=j1,j2
            do i=i1,i2
               if( abs(a4(3,i,j,k)-a4(2,i,j,k)) < -a4(4,i,j,k) ) then
                  fmin = a4(1,i,j,k)+0.25*(a4(3,i,j,k)-a4(2,i,j,k))**2/a4(4,i,j,k)+a4(4,i,j,k)*r12
                  if( fmin < 0. ) then
                     if(a4(1,i,j,k)<a4(3,i,j,k) .and. a4(1,i,j,k)<a4(2,i,j,k)) then
                        a4(3,i,j,k) = a4(1,i,j,k)
                        a4(2,i,j,k) = a4(1,i,j,k)
                        a4(4,i,j,k) = 0.
                     elseif(a4(3,i,j,k) > a4(2,i,j,k)) then
                        a4(4,i,j,k) = 3.*(a4(2,i,j,k)-a4(1,i,j,k))
                        a4(3,i,j,k) = a4(2,i,j,k) - a4(4,i,j,k)
                     else
                        a4(4,i,j,k) = 3.*(a4(3,i,j,k)-a4(1,i,j,k))
                        a4(2,i,j,k) = a4(3,i,j,k) - a4(4,i,j,k)
                     endif
                  endif
               endif
            enddo
         enddo
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
