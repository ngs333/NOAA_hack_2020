program Main

  use netcdf
  use fv_mapz_mod, only : CP_AIR, RVGAS, RDGAS
  use fv_mapz_mod, only : Lagrangian_to_Eulerian, fv_grid_type

  implicit none
include 'mpif.h'


  !namelist
  integer :: nx = 48
  integer :: ny = 48
  integer :: nz = 64
  integer :: num_iter = 100
  character(len=128) :: input_file = "fv_mapz_data.nc"

  namelist /main_nml/nx, ny, nz, num_iter, input_file


  integer :: ng = 3
  real    :: consv = 0.
  real    :: pdt = 1800.0
  integer :: nq = 10
  integer :: sphum=1
  real    :: kappa = RDGAS/CP_AIR
  real    :: zvir = rvgas/rdgas - 1.
  integer :: kord_mt = 10
  integer :: kord_wz = 10
  integer, allocatable :: kord_tr(:)
  integer :: kord_tm = 10
  real    :: ptop = 1.0
  logical :: do_omega = .false.

  integer   :: i, j, k, n, iter
  integer   :: isc, iec, jsc, jec
  integer   :: isd, ied, jsd, jed
  integer   :: npx, npy

  real, allocatable,       dimension(:) :: ak, bk
  real, allocatable,     dimension(:,:) :: ps, hs, te0_2d, ws
  real, allocatable,   dimension(:,:,:) :: pe, delp, pkz, pk, u, v, w, pt
  real, allocatable,   dimension(:,:,:) :: peln, ua, va, omga, te, delz
  real, allocatable,   dimension(:,:,:) :: pkz_out, pt_out
  real, allocatable, dimension(:,:,:,:) :: q, q_old

  real, allocatable,   dimension(:,:,:) :: u_old, v_old, pt_old, delz_old, delp_old
  real, allocatable,   dimension(:,:,:) :: pe_old, te_old, w_old, ua_old, va_old
  real, allocatable,   dimension(:,:,:) :: pk_old, omga_old, peln_old
  real, allocatable,     dimension(:,:) :: ps_old

  type(fv_grid_type) :: gridstruct

  logical :: fexist, opened, last_step
  integer :: unit, io_status
  integer :: info, myrank
  integer :: start_time, end_time, cr, cm
  real    :: rate, total_time
  integer :: ncid, did_axis, vid
  integer :: nxd, nyd, nzq, kq, iq
  integer :: nwat
  real    :: r_vir, cp, akap
  logical :: hydrostatic
  ! read namelist

! call mpi_init(info)
! call mpi_comm_rank(mpi_comm_world, myrank, info)

  CALL system_clock(count_rate=cr)
  CALL system_clock(count_max=cm)
  rate = REAL(cr)


  inquire(file='input.nml', exist = fexist)
  if(fexist) then
     do unit = 103, 512
        inquire(unit, opened=opened)
        if(.not. opened) exit
     enddo
     open(unit, file='input.nml', iostat = io_status)
     if( io_status > 0) call error_handler("Main program: failed at open input.nml")
     read(unit, main_nml, iostat = io_status)
     if( io_status > 0) call error_handler("Main program: failed at read main_nml")
     close(unit)
  endif


  !--- check if input file exist. If exist, read grid size from the file.
  inquire(file=trim(input_file), exist=fexist)
  if(.not. fexist) call error_handler("file "//trim(input_file)// " does not exist")

  print*, "NOTE from Main program: will read data from file "//trim(input_file)

  call handle_error(nf90_open(trim(input_file),nf90_nowrite,ncid))
  call get_data_i0d(ncid, 'npx', npx)
  call get_data_i0d(ncid, 'npy', npy)
  call get_data_i0d(ncid, 'km', nz)
  call get_data_i0d(ncid, 'nq', nq)
  nx = npx-1
  ny = npy-1
  print*, "nx,ny,nz,nq=", nx, ny, nz, nq

  isc = 1; iec = nx
  jsc = 1; jec = ny
  isd = isc-ng; ied = iec+ng
  jsd = jsc-ng; jed = jec+ng
 
  allocate(kord_tr(nq))
  allocate(ps(isd:ied,jsd:jed))
  allocate(hs(isd:ied,jsd:jed))
  allocate(pe(isc-1:iec+1,nz+1,jsc-1:jec+1))
  allocate(delp(isd:ied,jsd:jed,nz))
  allocate(pk(isc:iec,jsc:jec,nz+1))
  allocate(u(isd:ied,jsd:jed+1,nz))
  allocate(v(isd:ied+1,jsd:jed,nz))
  allocate(w(isd:ied,jsd:jed,nz))
  allocate(delz(isd:ied,jsd:jed,nz))
  allocate(pt(isd:ied,jsd:jed,nz))
  allocate(q(isd:ied,jsd:jed,nz,nq))
  allocate(peln(isc:iec,nz+1,jsc:jec))
  allocate(te0_2d(isc:iec,jsc:jec))
  allocate(omga(isd:ied,jsd:jed,nz))
  allocate(ws(isc:iec,jsc:jec))
  allocate(ak(nz+1), bk(nz+1))
  allocate(gridstruct%rsin2(isd:ied,jsd:jed))
  allocate(gridstruct%area(isd:ied,jsd:jed))
  allocate(gridstruct%cosa_s(isd:ied,jsd:jed))

  allocate(pkz(isc:iec,jsc:jec,nz))
  allocate(te(isd:ied,jsd:jed,nz))
  allocate(ua(isd:ied,jsd:jed,nz))
  allocate(va(isd:ied,jsd:jed,nz))
  allocate(pkz_out(isc:iec,jsc:jec,nz))
  allocate(pt_out(isd:ied,jsd:jed,nz))

  allocate(pe_old(isc-1:iec+1,nz+1,jsc-1:jec+1))
  allocate(delp_old(isd:ied,jsd:jed,nz))
  allocate(u_old(isd:ied,jsd:jed+1,nz))
  allocate(v_old(isd:ied+1,jsd:jed,nz))
  allocate(w_old(isd:ied,jsd:jed,nz))
  allocate(delz_old(isd:ied,jsd:jed,nz))
  allocate(pt_old(isd:ied,jsd:jed,nz))
  allocate(q_old(isd:ied,jsd:jed,nz,nq))
  allocate(te_old(isd:ied,jsd:jed,nz))
  allocate(ua_old(isd:ied,jsd:jed,nz))
  allocate(va_old(isd:ied,jsd:jed,nz))
  allocate(ps_old(isd:ied,jsd:jed))
  allocate(pk_old(isc:iec,jsc:jec,nz+1))
  allocate(omga_old(isd:ied,jsd:jed,nz))
  allocate(peln_old(isc:iec,nz+1,jsc:jec))

  call get_data_l0d(ncid, 'last_step', last_step)
  call get_data_r0d(ncid, 'consv', consv)
  call get_data_r2d(ncid, 'ps', ps)
  call get_data_r3d(ncid, 'pe', pe)
  call get_data_r3d(ncid, 'delp', delp)
  call get_data_r3d(ncid, 'pk', pk)
  call get_data_r0d(ncid, 'pdt', pdt)
  call get_data_i0d(ncid, 'nwat', nwat)
  call get_data_i0d(ncid, 'sphum', sphum)
  call get_data_r3d(ncid, 'u', u)
  call get_data_r3d(ncid, 'v', v)
  call get_data_r3d(ncid, 'w', w)
  call get_data_r3d(ncid, 'delz', delz)
  call get_data_r3d(ncid, 'pt', pt)
  call get_data_r3d(ncid, 'q', q(:,:,:,1))
  call get_data_r2d(ncid, 'hs', hs)
  call get_data_r0d(ncid, 'r_vir', r_vir)
  call get_data_r0d(ncid, 'cp', cp)
  call get_data_r0d(ncid, 'akap', akap)
  call get_data_i0d(ncid, 'kord_mt', kord_mt)
  call get_data_i0d(ncid, 'kord_wz', kord_wz)
  call get_data_i1d(ncid, 'kord_tr', kord_tr)
  call get_data_i0d(ncid, 'kord_tm', kord_tm)
  call get_data_r3d(ncid, 'peln', peln)
  call get_data_r2d(ncid, 'te0_2d', te0_2d)
  call get_data_r3d(ncid, 'omga', omga)
  call get_data_r2d(ncid, 'ws', ws)
  call get_data_r0d(ncid, 'ptop', ptop)
  call get_data_r1d(ncid, 'ak', ak)
  call get_data_r1d(ncid, 'bk', bk)
  call get_data_r2d(ncid, 'rsin2', gridstruct%rsin2)
  call get_data_r2d(ncid, 'cosa_s', gridstruct%cosa_s)
  call get_data_r2d(ncid, 'area', gridstruct%area)
  call get_data_l0d(ncid, 'hydrostatic', hydrostatic)
  call get_data_l0d(ncid, 'do_omega', do_omega)
  call get_data_r3d(ncid, 'pkz', pkz_out)
  call get_data_r3d(ncid, 'pt_out', pt_out)


  !--- output data
  pkz = 0
  te = 0

  u_old = u
  v_old = v
  pt_old = pt
  delz_old = delz
  delp_old = delp
  pe_old = pe
  te_old = te
  q_old = q
  w_old = w
  ua_old = ua
  va_old = va
  pk_old = pk
  ps_old = ps
  omga_old = omga
  peln_old = peln
  call system_clock(start_time)

  do iter = 1, num_iter
     pt = pt_old
     delz = delz_old
     delp = delp_old
     pe = pe_old
     q = q_old
     peln = peln_old
     call Lagrangian_to_Eulerian(last_step, consv, ps, pe, delp, pkz, pk,   &
                                   pdt, nz, isc,iec,jsc,jec, isd,ied,jsd,jed,       &
                      nq, nwat, sphum, u, v, w, delz, pt, q, hs, r_vir, cp,  &
                      akap, kord_mt, kord_wz, kord_tr, kord_tm,  peln, te0_2d,        &
                      ng, ua, va, omga, te, ws,   &
                      ptop, ak, bk, gridstruct, &
                      hydrostatic, do_omega)
  end do

call system_clock(end_time)
total_time = (end_time-start_time)/rate

  call compare_data("pkz", "pkz_out", pkz, pkz_out)
  call compare_data("pt", "pt_out", pt, pt_out)
  call print_chksum(pkz, "pkz")
  call print_chksum(pt, "pt")

  write(*,*) ' elapsed time (secs) = ', total_time


contains

  subroutine get_data_r0d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     real, intent(out) :: data

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_r0d

  subroutine get_data_r2d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     real, dimension(:,:), intent(out) :: data

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_r2d

  subroutine get_data_r3d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     real, dimension(:,:,:), intent(out) :: data

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_r3d

  subroutine get_data_i0d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     integer, intent(out) :: data

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_i0d

  subroutine get_data_i1d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     integer, intent(out) :: data(:)

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_i1d

  subroutine get_data_r1d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     real, intent(out) :: data(:)

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_r1d



  subroutine get_data_l0d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     logical, intent(out) :: data
     integer :: idata

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,idata))

     if(idata == 0) then
        data = .false.
     else if(idata==1) then
        data = .true.
     else
        call error_handler("field "//trim(field)//" value is neither 0 nor 1")
     endif

  end subroutine get_data_l0d


  subroutine compare_data(str1, str2, data1, data2)
     real, dimension(:,:,:), intent(in) :: data1, data2
     character(len=*), intent(in)       :: str1, str2

     write(*,'(A,ES18.10)' ) "Sum of "//trim(str1)//" = ", sum(data1)
     write(*,'(A,ES18.10)' ) "Sum of "//trim(str2)//" = ", sum(data2)
     write(*,'(A,ES18.10)' ) "max diff of "//trim(str1)//" = ", maxval(data1-data2)


  end subroutine compare_data


  subroutine error_handler(errmsg)
    character(len=*), intent(in) :: errmsg


    write(*,*) "ERROR: ", trim(errmsg)
    call ABORT()

  end subroutine error_handler

  subroutine handle_error(error_flag,err_string)
    ! Simple error handle for NetCDF
    integer,intent(in) :: error_flag
    character(*), intent(in),optional :: err_string
    if ( error_flag  /= nf90_noerr ) then
       write(*,*) 'FATAL ERROR:',nf90_strerror(error_flag)
       if (present(err_string)) write(*,*) trim(err_string)
       call ABORT()
    endif
  end subroutine handle_error

  subroutine print_chksum(data, str)
     real,             intent(in) :: data(:,:,:)
     character(len=*), intent(in) :: str
     integer(kind=8) :: idata(size(data))
     integer(kind=8) :: mold(1), chksum_data, gchksum
     real :: sum_data, gsum_data
     integer :: isc, iec, jsc, jec, error

     idata = transfer(data, mold)
     chksum_data = sum(idata)

     print*, trim(str)//" chksum is ", chksum_data

  end subroutine print_chksum


end program Main

