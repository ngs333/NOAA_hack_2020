!
! Argonne Leadership Computing Facility
! GFDL application kernel for Held Suarez benchmark
! Vitali Morozov (morozov@anl.gov)
! Aug 11, 2010
!
program Main
#ifdef _OPENACC
use openacc
#endif
use fyppm_mod

implicit none

include 'mpif.h'

character(132) linebuf
integer   :: i, j, k, info, myrank, mythread, omp_get_thread_num, iter

integer :: ifirst, ilast, isd, ied
integer :: jfirst, jlast, jsd, jed, js, je
integer :: layout(2) = (/1,1/)

integer :: km = 91
integer :: nx = 96
integer :: ny = 96

integer :: jord, npy, ng

real, allocatable, dimension (:,:,:) :: c
real, allocatable, dimension (:,:,:) :: q
real, allocatable, dimension (:,:,:) :: wk
real, allocatable, dimension (:,:,:) :: area, rarea
real, allocatable, dimension (:,:,:) :: al, bl, br, dq
real, allocatable, dimension (:,:) :: dya
real, allocatable, dimension (:,:,:) :: flux, dm, flx1

real :: a, b

real(8) :: time1, time2, elapsed

real :: ppm_limiter = 2.0
integer :: ndev
integer :: npes, nxc, nyc, ipos, jpos
integer :: istart, iend, jstart, jend
real    :: gsum
integer(KIND=8) :: gchksum

 call mpi_init(info)
 call mpi_comm_rank(mpi_comm_world, myrank, info)
 call MPI_COMM_SIZE(mpi_comm_world, npes, info )
ndev=0
#ifdef _OPENACC
 ndev = acc_get_num_devices(acc_device_nvidia)

 if( ndev > 0 .and. npes==ndev) then
    if(myrank==0) print*, "calling acc_set_device_num, npes, ndev=", npes, ndev
    call acc_set_device_num(myrank, acc_device_nvidia)
 endif
#endif
if(myrank==0) print*, "ndev = ", ndev, ", npes = ", npes
 call define_layout(nx,ny,npes,layout)
 if(myrank == 0) print*, "layout=", layout
  !make sure nx and ny is divisibor by layout
  if(mod(nx,layout(1)) .NE. 0) call error_handler("main: nx is not divisible by layout(1)")
  if(mod(ny,layout(2)) .NE. 0) call error_handler("main: nx is not divisible by layout(2)")
  nxc = nx/layout(1)
  nyc = ny/layout(2)

  ipos = mod(myrank,layout(1))
  jpos = myrank/layout(1)

  isd = ipos*nxc+1; ied = isd + nxc - 1
  jsd = jpos*nyc+1; jed = jsd + nyc - 1
  ifirst = isd; ilast  = ied
  jfirst = jsd; jlast  = jed
  js = jsd; je = jed; ng = 3; npy = ny + 1

  allocate( c   (isd:ied,js:je+1,1:km) )
  allocate( dya (isd:ied,js:je  ) )

  allocate( q   (ifirst:ilast,jfirst-ng:jlast+ng,1:km) )
  allocate( area(ifirst:ilast,jfirst   :jlast+1, 1:km ) )
  allocate( rarea(ifirst:ilast,jfirst   :jlast+1, 1:km ) )
  allocate(   wk(ifirst:ilast,jfirst   :jlast+1, 1:km ) )

  allocate( flux (isd:ied,js:je+1,1:km) )
flux(:,:,:) = 0.
  allocate( flx1 (isd:ied,js:je+1,1:km) )
  allocate( dm   (isd:ied,js-2:je+2,1:km) )

  allocate(al(ifirst:ilast,jfirst-1:jlast+2,km))
  allocate(bl(ifirst:ilast,jfirst-1:jlast+1,km))
  allocate(br(ifirst:ilast,jfirst-1:jlast+1,km))
  allocate(dq(ifirst:ilast,jfirst-3:jlast+2,km))

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

 rarea(:,:,:) = 0.
  do k = 1, km
     do j = jfirst, jlast+1
        do i = ifirst, ilast
           wk(i,j,k) = 1.0+k*1.e-2+j*2.e-4+i*3*1.e-6
           area(i,j,k) = 12.0+k*1.e-2+j*2.e-4+i*3*1.e-6
           rarea(i,j,k) = 1./area(i,j,k)
        enddo
     enddo
  enddo

  do j = js, je
     do i = isd, ied
        dya(i,j) = 2.0+i*1.e-2+j*2.e-4
     enddo
  enddo

    call CPU_TIME(time1)
!!$ACC enter data copyin(wk,area,c,q,dya,dm,flux)
!$ACC DATA copy(wk) copyin(area,rarea,c,q,dya,dm,flux,flx1,ppm_limiter,jord,jfirst,jlast,npy) copyout(al,bl,br,dq)
    do iter = 1, 10000
      if (mod(iter,10) .lt. 6) then
        jord = 9
      else
        jord = 12
      end if
!$omp parallel do
!$acc kernels async
      do k=1,km
        flux(:,:,k) = 1.0
        dm  (:,:,k) = 1.0
      enddo
!$acc end kernels
      call fyppm(c,  q,  flux, jord, ifirst, ilast, jfirst, jlast, npy, dm, ppm_limiter, &
                   dya, isd, ied, jsd, jed, js, je, ng, km, al, bl, br, dq )
!!$acc kernels loop collapse(3) present(wk,flux,flx1,rarea,area) async
!      do k=1,km
!        do j = jfirst, jlast+1
!          do i = ifirst, ilast
!            flx1(i,j,k) = flux(i,j,k)
!          end do
!        end do
!      end do
!!$acc end kernels

!$omp parallel do
!$acc kernels loop collapse(3) present(wk,flux,flx1,rarea,area) async
      do k=1,km
        do j = jfirst, jlast+1
          do i = ifirst, ilast
            wk(i,j,k) = wk(i,j,k) + flux(i,j,k)*rarea(i,j,k)
          end do
        end do
      end do
!$acc end kernels
!$acc wait
    end do
!$ACC END DATA
!!!$ACC exit data copyout(wk)
    call CPU_TIME(time2)

  istart = ifirst; iend = ilast
  jstart = jfirst; jend = jlast
  if(jlast==ny) jend=jend+1
  gsum = global_sum(wk(istart:iend,jstart:jend,:))
  gchksum = global_chksum(wk(istart:iend,jstart:jend,:))

  if(myrank == 0) then
!     write( linebuf,'(A)' ) 'expect   :  6.2082139432E+10'
     write( linebuf,'(A)' ) 'expect   :  9.6326547367E+10'

     write( 0, *) trim(linebuf)
     write( linebuf,'(A,ES18.10)' )  'got      :', gsum
     write( 0, *) trim(linebuf)
     write( linebuf,* )  'chksum      :', gchksum
     write( 0, *) trim(linebuf)
     print *, ' elapsed time (secs) = ', time2 - time1
  endif

    deallocate( c )
    deallocate( q )
    deallocate( area )
    deallocate( rarea )
    deallocate( wk )

    call mpi_finalize(info)

contains

subroutine error_handler(errmsg)
  character(len=*), intent(in) :: errmsg


  write(*,*) "ERROR: ", trim(errmsg)
  call MPI_ABORT()

end subroutine error_handler

subroutine define_layout(isz,jsz,ndivs,layout)
   integer, intent(in) :: isz, jsz
   integer, intent(in) :: ndivs
   integer, intent(out) :: layout(:)
   integer :: idiv, jdiv

   idiv = nint( sqrt(float(ndivs*isz)/jsz) )
   idiv = max(idiv,1) !for isz=1 line above can give 0
   do while( mod(ndivs,idiv).NE.0 )
      idiv = idiv - 1
   end do                 !will terminate at idiv=1 if not before
idiv=1
   jdiv = ndivs/idiv

   layout = (/ idiv, jdiv /)

end subroutine define_layout

function global_chksum(data)
   real,             intent(in) :: data(:,:,:)
   integer(kind=8) :: idata(size(data))
   integer(kind=8) :: mold(1), global_chksum, sum_data
   integer :: error

   idata = transfer(data, mold)
   sum_data = sum(idata)
   call MPI_ALLREDUCE(sum_data, global_chksum, 1, MPI_INTEGER8, MPI_SUM, mpi_comm_world, error)

end function global_chksum

function global_sum(data)
   real,             intent(in) :: data(:,:,:)
   real :: sum_data, global_sum
   integer :: error

   sum_data = sum(data)
   call MPI_ALLREDUCE(sum_data, global_sum, 1, MPI_REAL8, MPI_SUM, mpi_comm_world, error)

end function global_sum



end program Main

