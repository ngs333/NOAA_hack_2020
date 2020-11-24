!***********************************************************************
!*                   GNU Lesser General Public License                 
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it 
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or 
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be 
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty 
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.  
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'nh_utils' peforms non-hydrostatic computations.
!@author S. J. Lin, NOAA/GFDL
!>@todo Include moisture effect in pt

module nh_core_mod

   implicit none
   private

   real, parameter:: GRAV   = 9.8
   real, parameter:: RDGAS  = 287.04


   public Riem_Solver_c

   real, parameter:: dz_min = 2.
   real, parameter:: r3 = 1./3.
   integer :: call_count_to_write = 48

CONTAINS 

  subroutine Riem_Solver_c(ms,   dt,  is,   ie,   js, je, km,   ng,  &
                           akap, cappa, cp,  ptop, hs, w3,  pt, q_con, &
                           delp, gz,  pef,  ws, p_fac, a_imp, scale_m)

   integer, intent(in):: is, ie, js, je, ng, km
   integer, intent(in):: ms
   real, intent(in):: dt,  akap, cp, ptop, p_fac, a_imp, scale_m
   real, intent(in):: ws(is-ng:ie+ng,js-ng:je+ng)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, delp
   real, intent(in), dimension(is-ng:,js-ng:,1:):: q_con, cappa
   real, intent(in)::   hs(is-ng:ie+ng,js-ng:je+ng)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w3
! OUTPUT PARAMETERS 
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz
   real, intent(  out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: pef
! Local:
  real, dimension(is-1:ie+1,km  ):: dm, dz2, w2, pm2, gm2, cp2
  real, dimension(is-1:ie+1,km+1):: pem, pe2, peg
  real gama, rgrav
  integer i, j, k
  integer is1, ie1


    gama = 1./(1.-akap)
   rgrav = 1./grav

   is1 = is - 1
   ie1 = ie + 1

!$OMP parallel do default(none) shared(js,je,is1,ie1,km,delp,pef,ptop,gz,rgrav,w3,pt, &
!$OMP                                  a_imp,dt,gama,akap,ws,p_fac,scale_m,ms,hs,q_con,cappa) &
!$OMP                          private(cp2,gm2, dm, dz2, w2, pm2, pe2, pem, peg)
   do 2000 j=js-1, je+1

      do k=1,km
         do i=is1, ie1
            dm(i,k) = delp(i,j,k)
         enddo
      enddo

      do i=is1, ie1
         pef(i,j,1) = ptop                     ! full pressure at top
         pem(i,1) = ptop
#ifdef USE_COND
         peg(i,1) = ptop
#endif
      enddo

      do k=2,km+1
         do i=is1, ie1
            pem(i,k) = pem(i,k-1) + dm(i,k-1)
#ifdef USE_COND
! Excluding contribution from condensates:
            peg(i,k) = peg(i,k-1) + dm(i,k-1)*(1.-q_con(i,j,k-1))
#endif
         enddo
      enddo

      do k=1,km
         do i=is1, ie1
            dz2(i,k) = gz(i,j,k+1) - gz(i,j,k)
#ifdef USE_COND
            pm2(i,k) = (peg(i,k+1)-peg(i,k))/log(peg(i,k+1)/peg(i,k))

#ifdef MOIST_CAPPA
            cp2(i,k) = cappa(i,j,k)
            gm2(i,k) = 1. / (1.-cp2(i,k))
#endif

#else
            pm2(i,k) = dm(i,k)/log(pem(i,k+1)/pem(i,k))
#endif
             dm(i,k) = dm(i,k) * rgrav
             w2(i,k) = w3(i,j,k)
         enddo
      enddo

      if ( a_imp > 0.5 ) then
           call SIM1_solver(dt, is1, ie1, km, rdgas, gama, gm2, cp2, akap, pe2,  &
                            dm, pm2, pem, w2, dz2, pt(is1:ie1,j,1:km), ws(is1,j), p_fac)
      endif

      do k=2,km+1
         do i=is1, ie1
            pef(i,j,k) = pe2(i,k) + pem(i,k)  ! add hydrostatic full-component
         enddo
      enddo
! Compute Height * grav (for p-gradient computation)
      do i=is1, ie1
         gz(i,j,km+1) = hs(i,j)
      enddo

      do k=km,1,-1
         do i=is1, ie1
            gz(i,j,k) = gz(i,j,k+1) - dz2(i,k)*grav
         enddo
      enddo

2000  continue



  end subroutine Riem_Solver_c


 subroutine SIM1_solver(dt,  is,  ie, km, rgas, gama, gm2, cp2, kappa, pe, dm2,   &
                        pm2, pem, w2, dz2, pt2, ws, p_fac)
   integer, intent(in):: is, ie, km
   real,    intent(in):: dt, rgas, gama, kappa, p_fac
   real, intent(in), dimension(is:ie,km):: dm2, pt2, pm2, gm2, cp2
   real, intent(in )::  ws(is:ie)
   real, intent(in ), dimension(is:ie,km+1):: pem
   real, intent(out)::  pe(is:ie,km+1)
   real, intent(inout), dimension(is:ie,km):: dz2, w2
! Local
   real, dimension(is:ie,km  ):: aa, bb, dd, w1, g_rat, gam
   real, dimension(is:ie,km+1):: pp
   real, dimension(is:ie):: p1, bet
   real t1g, rdt, capa1
   integer i, k


#ifdef MOIST_CAPPA
      t1g = 2.*dt*dt
#else
      t1g = gama * 2.*dt*dt
#endif
      rdt = 1. / dt
    capa1 = kappa - 1.

    do k=1,km
       do i=is, ie
#ifdef MOIST_CAPPA
          pe(i,k) = exp(gm2(i,k)*log(-dm2(i,k)/dz2(i,k)*rgas*pt2(i,k))) - pm2(i,k)
#else
          pe(i,k) = exp(gama*log(-dm2(i,k)/dz2(i,k)*rgas*pt2(i,k))) - pm2(i,k)
#endif
          w1(i,k) = w2(i,k)
       enddo
    enddo

    do k=1,km-1
       do i=is, ie
          g_rat(i,k) = dm2(i,k)/dm2(i,k+1)
             bb(i,k) = 2.*(1.+g_rat(i,k))
             dd(i,k) = 3.*(pe(i,k) + g_rat(i,k)*pe(i,k+1))
       enddo
    enddo

    do i=is, ie
         bet(i) = bb(i,1)
        pp(i,1) = 0.
        pp(i,2) = dd(i,1) / bet(i)
       bb(i,km) = 2.
       dd(i,km) = 3.*pe(i,km)
    enddo

    do k=2,km
      do i=is, ie
          gam(i,k) =  g_rat(i,k-1) / bet(i)
            bet(i) =  bb(i,k) - gam(i,k)
         pp(i,k+1) = (dd(i,k) - pp(i,k) ) / bet(i)
      enddo
    enddo

    do k=km, 2, -1
       do i=is, ie
          pp(i,k) = pp(i,k) - gam(i,k)*pp(i,k+1)
       enddo
    enddo

! Start the w-solver
    do k=2, km
       do i=is, ie
#ifdef MOIST_CAPPA
          aa(i,k) = t1g*0.5*(gm2(i,k-1)+gm2(i,k))/(dz2(i,k-1)+dz2(i,k)) * (pem(i,k)+pp(i,k))
#else
          aa(i,k) = t1g/(dz2(i,k-1)+dz2(i,k)) * (pem(i,k)+pp(i,k))
#endif
       enddo
    enddo
    do i=is, ie
       bet(i)  = dm2(i,1) - aa(i,2)
       w2(i,1) = (dm2(i,1)*w1(i,1) + dt*pp(i,2)) / bet(i)
    enddo
    do k=2,km-1
       do i=is, ie
          gam(i,k) = aa(i,k) / bet(i)
            bet(i) =  dm2(i,k) - (aa(i,k) + aa(i,k+1) + aa(i,k)*gam(i,k))
           w2(i,k) = (dm2(i,k)*w1(i,k)+dt*(pp(i,k+1)-pp(i,k))-aa(i,k)*w2(i,k-1)) / bet(i)
       enddo
    enddo
    do i=is, ie
#ifdef MOIST_CAPPA
           p1(i) = t1g*gm2(i,km)/dz2(i,km)*(pem(i,km+1)+pp(i,km+1))
#else
           p1(i) = t1g/dz2(i,km)*(pem(i,km+1)+pp(i,km+1))
#endif
       gam(i,km) = aa(i,km) / bet(i)
          bet(i) =  dm2(i,km) - (aa(i,km)+p1(i) + aa(i,km)*gam(i,km))
        w2(i,km) = (dm2(i,km)*w1(i,km)+dt*(pp(i,km+1)-pp(i,km))-p1(i)*ws(i)-aa(i,km)*w2(i,km-1))/bet(i)
    enddo
    do k=km-1, 1, -1
       do i=is, ie
          w2(i,k) = w2(i,k) - gam(i,k+1)*w2(i,k+1)
       enddo
    enddo

    do i=is, ie
       pe(i,1) = 0.
    enddo
    do k=1,km
       do i=is, ie
          pe(i,k+1) = pe(i,k) + dm2(i,k)*(w2(i,k)-w1(i,k))*rdt
       enddo
    enddo

    do i=is, ie
           p1(i) = ( pe(i,km) + 2.*pe(i,km+1) )*r3
#ifdef MOIST_CAPPA
       dz2(i,km) = -dm2(i,km)*rgas*pt2(i,km)*exp((cp2(i,km)-1.)*log(max(p_fac*pm2(i,km),p1(i)+pm2(i,km))))
#else
       dz2(i,km) = -dm2(i,km)*rgas*pt2(i,km)*exp(capa1*log(max(p_fac*pm2(i,km),p1(i)+pm2(i,km))))
#endif
    enddo

    do k=km-1, 1, -1
       do i=is, ie
          p1(i) = (pe(i,k) + bb(i,k)*pe(i,k+1) + g_rat(i,k)*pe(i,k+2))*r3 - g_rat(i,k)*p1(i)
#ifdef MOIST_CAPPA
          dz2(i,k) = -dm2(i,k)*rgas*pt2(i,k)*exp((cp2(i,k)-1.)*log(max(p_fac*pm2(i,k),p1(i)+pm2(i,k))))
#else
          dz2(i,k) = -dm2(i,k)*rgas*pt2(i,k)*exp(capa1*log(max(p_fac*pm2(i,k),p1(i)+pm2(i,k))))
#endif
       enddo
    enddo

 end subroutine SIM1_solver



end module nh_core_mod
