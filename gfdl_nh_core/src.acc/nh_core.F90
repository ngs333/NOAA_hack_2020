!>@brief The module 'nh_utils' peforms non-hydrostatic computations.
!@author S. J. Lin, NOAA/GFDL
!>@todo Include moisture effect in pt

#ifdef _OPENACC
#define _DEVICE_ ,device
#else
#define _DEVICE_
#endif

module nh_core_mod

   implicit none
   private

   real, parameter:: GRAV   = 9.8
   real, parameter:: RDGAS  = 287.04


   public Riem_Solver_c, allocate_data

   real, parameter:: dz_min = 2.
   real, parameter:: r3 = 1./3.

   real, allocatable, dimension(:,:,:) :: dm, dz2, w2, pm2, gm2, cp2
   real, allocatable, dimension(:,:,:) :: pem, pe2, peg
   real, allocatable, dimension(:,:,:) :: aa, bb, dd, w1, g_rat, gam, pp
   real, allocatable, dimension(:,:  ) :: p1, bet

!$ACC declare create(dm, dz2, w2, pm2, gm2, cp2,pem, pe2, peg, &
!$ACC                aa, bb, dd, w1, g_rat, gam,pp,p1,bet)

CONTAINS

   subroutine allocate_data(is,ie,js,je,km)
    integer, intent(in) :: is,ie,js,je,km

    allocate(dm(is-1:ie+1,js-1:je+1,km) )
    allocate(dz2(is-1:ie+1,js-1:je+1,km) )
    allocate(w2(is-1:ie+1,js-1:je+1,km) )
    allocate(pm2(is-1:ie+1,js-1:je+1,km) )
    allocate(gm2(is-1:ie+1,js-1:je+1,km) )
    allocate(cp2(is-1:ie+1,js-1:je+1,km) )
    allocate(pem(is-1:ie+1,js-1:je+1,km+1) )
    allocate(pe2(is-1:ie+1,js-1:je+1,km+1) )
    allocate(peg(is-1:ie+1,js-1:je+1,km+1) )
    allocate(aa (is-1:ie+1,js-1:je+1,km)   )
    allocate(bb (is-1:ie+1,js-1:je+1,km)   )
    allocate(dd (is-1:ie+1,js-1:je+1,km)   )
    allocate(w1 (is-1:ie+1,js-1:je+1,km)   )
    allocate(g_rat (is-1:ie+1,js-1:je+1,km))
    allocate(gam (is-1:ie+1,js-1:je+1,km)  )
    allocate(pp (is-1:ie+1,js-1:je+1,km+1) )
    allocate(p1 (is-1:ie+1,js-1:je+1) )
    allocate(bet (is-1:ie+1,js-1:je+1) )

  end subroutine allocate_data

  subroutine Riem_Solver_c(ms,   dt,  is,   ie,   js, je, km,   ng,  &
                           akap, cappa, cp,  ptop, hs, w3,  pt, q_con, &
                           delp, gz,  pef,  ws, p_fac, a_imp, scale_m)

   integer, intent(in):: is, ie, js, je, ng, km
   integer, intent(in):: ms
   real, intent(in):: dt,  akap, cp, ptop, p_fac, a_imp, scale_m
   real _DEVICE_, intent(in):: ws(is-ng:ie+ng,js-ng:je+ng)
   real _DEVICE_, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, delp
   real _DEVICE_, intent(in), dimension(is-ng:,js-ng:,1:):: q_con, cappa
   real _DEVICE_, intent(in)::   hs(is-ng:ie+ng,js-ng:je+ng)
   real _DEVICE_, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w3
! OUTPUT PARAMETERS
   real _DEVICE_, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz
   real _DEVICE_, intent(  out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: pef
! Local:
  real gama, rgrav
  integer i, j, k
  integer is1, ie1, js1, je1, isd, ied, jsd, jed


    gama = 1./(1.-akap)
   rgrav = 1./grav

   is1 = is - 1
   ie1 = ie + 1
   js1 = js - 1
   je1 = je + 1
   isd = is - ng
   ied = ie + ng
   jsd = js - ng
   jed = je + ng

!$ACC kernels present(ms,dt,is,ie,js,je,km,ng,akap,cappa, &
!$ACC                 cp,ptop,hs,w3,pt,q_con,delp,gz,pef,ws,       &
!$ACC                 p_fac,a_imp,scale_m,dm, dz2, w2, pm2, gm2, cp2,pem, pe2,peg)
!$acc loop collapse(2)

   do k=1,km
      do j=js1,je1
         do i=is1, ie1
            dm(i,j,k) = delp(i,j,k)
         enddo
      enddo
   enddo

   do j=js1,je1
      do i=is1, ie1
         pef(i,j,1) = ptop                     ! full pressure at top
         pem(i,j,1) = ptop
#ifdef USE_COND
         peg(i,j,1) = ptop
#endif
      enddo
   enddo

   do k=2,km+1
      do j=js1,je1
         do i=is1, ie1
            pem(i,j,k) = pem(i,j,k-1) + dm(i,j,k-1)
#ifdef USE_COND
! Excluding contribution from condensates:
            peg(i,j,k) = peg(i,j,k-1) + dm(i,j,k-1)*(1.-q_con(i,j,k-1))
#endif
         enddo
      enddo
   enddo

   do k=1,km
      do j=js1,je1
         do i=is1, ie1
            dz2(i,j,k) = gz(i,j,k+1) - gz(i,j,k)
#ifdef USE_COND
            pm2(i,j,k) = (peg(i,j,k+1)-peg(i,j,k))/log(peg(i,j,k+1)/peg(i,j,k))

#ifdef MOIST_CAPPA
            cp2(i,j,k) = cappa(i,j,k)
            gm2(i,j,k) = 1. / (1.-cp2(i,j,k))
#endif

#else
            pm2(i,j,k) = dm(i,j,k)/log(pem(i,j,k+1)/pem(i,j,k))
#endif
             dm(i,j,k) = dm(i,j,k) * rgrav
             w2(i,j,k) = w3(i,j,k)
         enddo
      enddo
   enddo
!$ACC end kernels
      if ( a_imp > 0.5 ) then
           call SIM1_solver(dt, is1, ie1, js1, je1, isd, ied, jsd, jed, km, rdgas, gama, gm2, cp2, akap, pe2,  &
                            dm, pm2, pem, w2, dz2, pt, ws, p_fac)
      endif
!$ACC kernels present(ms,dt,km,hs,gz,pef,pe2,pem,dz2)
    do k=2,km+1
      do j=js1,je1
         do i=is1, ie1
            pef(i,j,k) = pe2(i,j,k) + pem(i,j,k)  ! add hydrostatic full-component
         enddo
      enddo
    enddo
! Compute Height * grav (for p-gradient computation)
    do j=js1,je1
      do i=is1, ie1
         gz(i,j,km+1) = hs(i,j)
      enddo
    enddo

    do k=km,1,-1
      do j=js1,je1
         do i=is1, ie1
            gz(i,j,k) = gz(i,j,k+1) - dz2(i,j,k)*grav
         enddo
      enddo
    enddo
!$ACC end kernels



  end subroutine Riem_Solver_c


 subroutine SIM1_solver(dt,  is,  ie, js, je, isd, ied, jsd, jed, km, rgas, gama, gm2, cp2, kappa, pe, dm2,   &
                        pm2, pem, w2, dz2, pt2, ws, p_fac)
   integer, intent(in):: is, ie, js, je, isd, ied, jsd, jed, km
   real,    intent(in):: dt, rgas, gama, kappa, p_fac
   real _DEVICE_,  intent(in), dimension(is:ie,js:je,km):: dm2, pm2, gm2, cp2
   real _DEVICE_, intent(in), dimension(isd:ied,jsd:jed,km) :: pt2
   real _DEVICE_,  intent(in )::  ws(isd:ied,jsd:jed)
   real _DEVICE_, intent(in ), dimension(is:ie,js:je,km+1):: pem
   real _DEVICE_, intent(out)::  pe(is:ie,js:je,km+1)
   real _DEVICE_, intent(inout), dimension(is:ie,js:je,km):: dz2, w2
! Local
!   real _DEVICE_, dimension(is:ie,js:je,km  ):: aa, bb, dd, w1, g_rat, gam
!   real _DEVICE_, dimension(is:ie,js:je,km+1):: pp
!   real _DEVICE_, dimension(is:ie,js:je):: p1, bet
   real t1g, rdt, capa1
   integer i, j, k


!$ACC kernels
#ifdef MOIST_CAPPA
      t1g = 2.*dt*dt
#else
      t1g = gama * 2.*dt*dt
#endif
      rdt = 1. / dt
    capa1 = kappa - 1.

  do k=1,km
    do j=js,je
       do i=is, ie
#ifdef MOIST_CAPPA
          pe(i,j,k) = exp(gm2(i,j,k)*log(-dm2(i,j,k)/dz2(i,j,k)*rgas*pt2(i,j,k))) - pm2(i,j,k)
#else
          pe(i,j,k) = exp(gama*log(-dm2(i,j,k)/dz2(i,j,k)*rgas*pt2(i,j,k))) - pm2(i,j,k)
#endif
          w1(i,j,k) = w2(i,j,k)
       enddo
    enddo
  enddo

  do k=1,km-1
    do j=js,je
       do i=is, ie
          g_rat(i,j,k) = dm2(i,j,k)/dm2(i,j,k+1)
             bb(i,j,k) = 2.*(1.+g_rat(i,j,k))
             dd(i,j,k) = 3.*(pe(i,j,k) + g_rat(i,j,k)*pe(i,j,k+1))
       enddo
    enddo
  enddo
!$acc loop collapse(2)
  do j=js,je
    do i=is, ie
         bet(i,j) = bb(i,j,1)
        pp(i,j,1) = 0.
        pp(i,j,2) = dd(i,j,1) / bet(i,j)
       bb(i,j,km) = 2.
       dd(i,j,km) = 3.*pe(i,j,km)
    enddo
  enddo

  do k=2,km
    do j=js,je
      do i=is, ie
          gam(i,j,k) =  g_rat(i,j,k-1) / bet(i,j)
            bet(i,j) =  bb(i,j,k) - gam(i,j,k)
         pp(i,j,k+1) = (dd(i,j,k) - pp(i,j,k) ) / bet(i,j)
      enddo
    enddo
  enddo

  do k=km, 2, -1
    do j=js,je
       do i=is, ie
          pp(i,j,k) = pp(i,j,k) - gam(i,j,k)*pp(i,j,k+1)
       enddo
    enddo
  enddo

! Start the w-solver
  do k=2, km
    do j=js,je
       do i=is, ie
#ifdef MOIST_CAPPA
          aa(i,j,k) = t1g*0.5*(gm2(i,j,k-1)+gm2(i,j,k))/(dz2(i,j,k-1)+dz2(i,j,k)) * (pem(i,j,k)+pp(i,j,k))
#else
          aa(i,j,k) = t1g/(dz2(i,j,k-1)+dz2(i,j,k)) * (pem(i,j,k)+pp(i,j,k))
#endif
       enddo
    enddo
  enddo
  do j=js,je
    do i=is, ie
       bet(i,j)  = dm2(i,j,1) - aa(i,j,2)
       w2(i,j,1) = (dm2(i,j,1)*w1(i,j,1) + dt*pp(i,j,2)) / bet(i,j)
    enddo
  enddo

  do k=2,km-1
    do j=js,je
       do i=is, ie
          gam(i,j,k) = aa(i,j,k) / bet(i,j)
            bet(i,j) =  dm2(i,j,k) - (aa(i,j,k) + aa(i,j,k+1) + aa(i,j,k)*gam(i,j,k))
           w2(i,j,k) = (dm2(i,j,k)*w1(i,j,k)+dt*(pp(i,j,k+1)-pp(i,j,k))-aa(i,j,k)*w2(i,j,k-1)) / bet(i,j)
       enddo
    enddo
  enddo

  do j=js,je
    do i=is, ie
#ifdef MOIST_CAPPA
           p1(i,j) = t1g*gm2(i,j,km)/dz2(i,j,km)*(pem(i,j,km+1)+pp(i,j,km+1))
#else
           p1(i,j) = t1g/dz2(i,j,km)*(pem(i,j,km+1)+pp(i,j,km+1))
#endif
       gam(i,j,km) = aa(i,j,km) / bet(i,j)
          bet(i,j) =  dm2(i,j,km) - (aa(i,j,km)+p1(i,j) + aa(i,j,km)*gam(i,j,km))
        w2(i,j,km) = (dm2(i,j,km)*w1(i,j,km)+dt*(pp(i,j,km+1)-pp(i,j,km))-p1(i,j)*ws(i,j)-aa(i,j,km)*w2(i,j,km-1))/bet(i,j)
    enddo
  enddo

  do k=km-1, 1, -1
    do j=js,je
       do i=is, ie
          w2(i,j,k) = w2(i,j,k) - gam(i,j,k+1)*w2(i,j,k+1)
       enddo
    enddo
  enddo
  do j=js,je
    do i=is, ie
       pe(i,j,1) = 0.
    enddo
  enddo

  do k=1,km
    do j=js,je
       do i=is, ie
          pe(i,j,k+1) = pe(i,j,k) + dm2(i,j,k)*(w2(i,j,k)-w1(i,j,k))*rdt
       enddo
    enddo
  enddo

  do j=js,je
    do i=is, ie
           p1(i,j) = ( pe(i,j,km) + 2.*pe(i,j,km+1) )*r3
#ifdef MOIST_CAPPA
       dz2(i,j,km) = -dm2(i,j,km)*rgas*pt2(i,j,km)*exp((cp2(i,j,km)-1.)*log(max(p_fac*pm2(i,j,km),p1(i,j)+pm2(i,j,km))))
#else
       dz2(i,j,km) = -dm2(i,j,km)*rgas*pt2(i,j,km)*exp(capa1*log(max(p_fac*pm2(i,j,km),p1(i,j)+pm2(i,j,km))))
#endif
    enddo
  enddo

  do k=km-1, 1, -1
    do j=js,je
       do i=is, ie
          p1(i,j) = (pe(i,j,k) + bb(i,j,k)*pe(i,j,k+1) + g_rat(i,j,k)*pe(i,j,k+2))*r3 - g_rat(i,j,k)*p1(i,j)
#ifdef MOIST_CAPPA
          dz2(i,j,k) = -dm2(i,j,k)*rgas*pt2(i,j,k)*exp((cp2(i,j,k)-1.)*log(max(p_fac*pm2(i,j,k),p1(i,j)+pm2(i,j,k))))
#else
          dz2(i,j,k) = -dm2(i,j,k)*rgas*pt2(i,j,k)*exp(capa1*log(max(p_fac*pm2(i,j,k),p1(i,j)+pm2(i,j,k))))
#endif
       enddo
    enddo
  enddo
!$ACC end kernels

 end subroutine SIM1_solver



end module nh_core_mod
