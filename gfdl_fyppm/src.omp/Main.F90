!
! Argonne Leadership Computing Facility
! GFDL application kernel for Held Suarez benchmark
! Vitali Morozov (morozov@anl.gov)
! Aug 11, 2010
!
program Main

use fyppm_mod

implicit none

include 'mpif.h'

character(132) linebuf
integer   :: i, j, k, info, myrank, mythread, omp_get_thread_num, iter

integer :: ifirst, ilast, isd, ied
integer :: jfirst, jlast, jsd, jed, js, je

integer, parameter :: km = 91
integer, parameter :: nx = 96
integer, parameter :: ny = 96

integer :: jord, npy, ng

real, allocatable, dimension (:,:,:) :: c
real, allocatable, dimension (:,:,:) :: q
real, allocatable, dimension (:,:,:) :: wk
real, allocatable, dimension (:,:,:) :: area
real, allocatable, dimension (:,:) :: dya

real :: flux(nx, ny+1)
real :: dm(nx,-1:ny+2)

real :: a, b

real(8) :: time1, time2, elapsed

real :: ppm_limiter = 2.0

print*, "begin the run"

! call mpi_init(info)
! call mpi_comm_rank(mpi_comm_world, myrank, info)


  isd    = 1; ied    = nx
  jsd    = 1; jed    = ny
  ifirst = 1; ilast  = nx
  jfirst = 1; jlast  = ny

  js = 1; je = ny; ng = 3; npy = je + 1

  allocate( c   (isd:ied,js:je+1,1:km) )
  allocate( dya (isd:ied,js:je  ) )

  allocate( q   (ifirst:ilast,jfirst-ng:jlast+ng,1:km) )
  allocate( area(ifirst:ilast,jfirst   :jlast+1, 1:km ) )
  allocate(   wk(ifirst:ilast,jfirst   :jlast+1, 1:km ) )

  do k = 1, km
    do j=jfirst-ng,jlast+ng
      do i = ifirst,ilast
        q(i,j,k) = i + j + k
      enddo
    enddo
  enddo

  do k = 1, km
     do j = js, je+1
        do i = isd, ied
           c(i,j,k) = 1.0+i*1.e-2+j*2.e-4+k*3*1.e-6
        enddo
     enddo
  enddo

  do k = 1, km
     do j = jfirst, jlast+1
        do i = ifirst, ilast
           wk(i,j,k) = 1.0+k*1.e-2+j*2.e-4+i*3*1.e-6
           area(i,j,k) = 12.0+k*1.e-2+j*2.e-4+i*3*1.e-6
        enddo
     enddo
  enddo

  do j = js, je
     do i = isd, ied
        dya(i,j) = 2.0+i*1.e-2+j*2.e-4
     enddo
  enddo

    call CPU_TIME(time1)
!$ACC DATA copy(wk) copyin(area,c,q,dya,dm, flux)
    do iter = 1, 10000
      if (mod(iter,10) .lt. 6) then
        jord = 9
      else
        jord = 12
      end if
!$omp parallel do private(flux, dm)
!!$acc parallel loop gang private(flux,dm)
!$acc kernels
      do k=1,km
        flux(:,:) = 1.0
        dm  (:,:) = 1.0
        call fyppm(c(isd,js,k),  q(ifirst,jfirst-ng,k),  flux, jord, jfirst, jlast, npy, dm, ppm_limiter, &
                   dya, isd, ied, jsd, jed, js, je, ng )

        do j = jfirst, jlast+1
          do i = ifirst, ilast
            wk(i,j,k) = wk(i,j,k) + flux(i,j)/area(i,j,k)
          end do
        end do
      end do
!$acc end kernels
!!$acc end parallel
    end do
!$ACC END DATA
    call CPU_TIME(time2)


 print*, "end of the run"

!     write( linebuf,'(A)' ) 'expect   :  6.2082139432E+10'
     write( linebuf,'(A)' ) 'expect   :  9.6326547367E+10'
     write( 0, *) trim(linebuf)
     write( linebuf,'(A,ES18.10)' )  'got      :', SUM( wk(:,:,:) )
     write( 0, *) trim(linebuf)
     print *, ' elapsed time (secs) = ', time2 - time1
     write( 0,*) trim(linebuf)

    deallocate( c )
    deallocate( q )
    deallocate( area )
    deallocate( wk )

!    call mpi_finalize(info)

end program Main

