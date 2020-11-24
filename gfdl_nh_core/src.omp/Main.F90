program Main

  use netcdf
  use nh_core_mod, only :  Riem_Solver_c

  implicit none

  !namelist
  integer :: num_iter = 100
  character(len=128) :: riem_solver_c_file = "Riem_Solver_c_data.nc"

  namelist /main_nml/num_iter

  integer :: nx=0, ny=0, nz=0
  integer :: ng = 3
  integer :: ms
  real    :: dt, akap, cp, ptop, scale_m, p_fac, a_imp

  integer   :: npx, npy
  integer   :: isc, iec, jsc, jec
  integer   :: isd, ied, jsd, jed


  integer   :: i, j, k, iter


  real, allocatable,     dimension(:,:) :: hs, ws
  real, allocatable,   dimension(:,:,:) :: cappa, w3, pt, q_con, delp, gz, pef, gz_out, pef_out

  real, allocatable,     dimension(:,:) :: hs_old, ws_old
  real, allocatable,   dimension(:,:,:) :: cappa_old, w3_old, pt_old, q_con_old, delp_old, gz_old



  logical :: fexist, opened
  integer :: unit, io_status
  integer :: info, myrank
  integer :: start_time, end_time, cr, cm
  real    :: rate, total_time
  integer :: ncid, did_axis, vid

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
  inquire(file=trim(riem_solver_c_file), exist=fexist)
  if(.not. fexist) call error_handler("file "//trim(riem_solver_c_file)// " does not exist")

  print*, "NOTE from Main program: read riem_solver_c data from file "//trim(riem_solver_c_file)

  call handle_error(nf90_open(trim(riem_solver_c_file),nf90_nowrite,ncid))
  call handle_error(nf90_inq_varid(ncid,'npx',vid))
  call handle_error(nf90_get_var(ncid,vid,npx))
  call handle_error(nf90_inq_varid(ncid,'npy',vid))
  call handle_error(nf90_get_var(ncid,vid,npy))
  call handle_error(nf90_inq_varid(ncid,'km',vid))
  call handle_error(nf90_get_var(ncid,vid,nz))
  call handle_error(nf90_inq_varid(ncid,'ng',vid))
  call handle_error(nf90_get_var(ncid,vid,ng))

  nx = npx-1; ny = npy-1
  print*, "nx,ny,nz", nx, ny, nz, ng


  isc = 1; iec = nx
  jsc = 1; jec = ny
  isd = isc-ng; ied = iec+ng
  jsd = jsc-ng; jed = jec+ng

  allocate(cappa(isd:ied,jsd:jed,nz))
  allocate(hs(isd:ied,jsd:jed))
  allocate(w3(isd:ied,jsd:jed,nz))
  allocate(pt(isd:ied,jsd:jed,nz))
  allocate(q_con(isd:ied,jsd:jed,nz))
  allocate(delp(isd:ied,jsd:jed,nz))
  allocate(gz(isd:ied,jsd:jed,nz+1))
  allocate(ws(isd:ied,jsd:jed))
  allocate(gz_out(isd:ied,jsd:jed,nz+1))
  allocate(pef(isd:ied,jsd:jed,nz+1))
  allocate(pef_out(isd:ied,jsd:jed,nz+1))

  allocate(cappa_old(isd:ied,jsd:jed,nz))
  allocate(hs_old(isd:ied,jsd:jed))
  allocate(w3_old(isd:ied,jsd:jed,nz))
  allocate(pt_old(isd:ied,jsd:jed,nz))
  allocate(q_con_old(isd:ied,jsd:jed,nz))
  allocate(delp_old(isd:ied,jsd:jed,nz))
  allocate(gz_old(isd:ied,jsd:jed,nz+1))
  allocate(ws_old(isd:ied,jsd:jed))

  pef = 1.e8
  call handle_error(nf90_inq_varid(ncid,'ms',vid))
  call handle_error(nf90_get_var(ncid,vid,ms))
  call handle_error(nf90_inq_varid(ncid,'dt',vid))
  call handle_error(nf90_get_var(ncid,vid,dt))
  call handle_error(nf90_inq_varid(ncid,'akap',vid))
  call handle_error(nf90_get_var(ncid,vid,akap))
!  call handle_error(nf90_inq_varid(ncid,'cappa',vid))
!  call handle_error(nf90_get_var(ncid,vid,cappa))
  call handle_error(nf90_inq_varid(ncid,'cp',vid))
  call handle_error(nf90_get_var(ncid,vid,cp))
  call handle_error(nf90_inq_varid(ncid,'ptop',vid))
  call handle_error(nf90_get_var(ncid,vid,ptop))
  call handle_error(nf90_inq_varid(ncid,'hs',vid))
  call handle_error(nf90_get_var(ncid,vid,hs))
  call handle_error(nf90_inq_varid(ncid,'w3',vid))
  call handle_error(nf90_get_var(ncid,vid,w3))
  call handle_error(nf90_inq_varid(ncid,'pt',vid))
  call handle_error(nf90_get_var(ncid,vid,pt))
!  call handle_error(nf90_inq_varid(ncid,'q_con',vid))
!  call handle_error(nf90_get_var(ncid,vid,q_con))
  call handle_error(nf90_inq_varid(ncid,'delp',vid))
  call handle_error(nf90_get_var(ncid,vid,delp))
  call handle_error(nf90_inq_varid(ncid,'gz',vid))
  call handle_error(nf90_get_var(ncid,vid,gz))
  call handle_error(nf90_inq_varid(ncid,'ws',vid))
  call handle_error(nf90_get_var(ncid,vid,ws))
  call handle_error(nf90_inq_varid(ncid,'p_fac',vid))
  call handle_error(nf90_get_var(ncid,vid,p_fac))
  call handle_error(nf90_inq_varid(ncid,'a_imp',vid))
  call handle_error(nf90_get_var(ncid,vid,a_imp))
  call handle_error(nf90_inq_varid(ncid,'scale_m',vid))
  call handle_error(nf90_get_var(ncid,vid,scale_m))
  call handle_error(nf90_inq_varid(ncid,'gz_out',vid))
  call handle_error(nf90_get_var(ncid,vid,gz_out))
  call handle_error(nf90_inq_varid(ncid,'pef',vid))
  call handle_error(nf90_get_var(ncid,vid,pef_out))


!  cappa_old = cappa
  hs_old = hs
  w3_old = w3
  pt_old = pt
!  q_con_old = q_con
  delp_old = delp
  gz_old = gz
  ws_old = ws

  call system_clock(start_time)

  do iter = 1, num_iter
!     cappa = cappa_old
     hs = hs_old
     w3 = w3_old
     pt = pt_old
!     q_con = q_con_old
     delp = delp_old
     gz = gz_old
     ws = ws_old

     call Riem_Solver_c(ms, dt, isc, iec, jsc, jec, nz, ng, akap,   &
                          cappa, cp, ptop, hs, w3,  pt, q_con, &
                          delp, gz, pef, ws, p_fac, a_imp, scale_m)

  end do


  call system_clock(end_time)
  total_time = (end_time-start_time)/rate

  !--- check if data match the output from model
  write(*,'(A,ES18.10)' ) "Sum of gz_out  = ", sum(gz_out)
  write(*,'(A,ES18.10)' ) "Sum of gz      = ", sum(gz)
  write(*,'(A,ES18.10)' ) "Sum of pef_out = ", sum(pef_out)
  write(*,'(A,ES18.10)' ) "Sum of pef     = ", sum(pef)
  write(*,'(A,ES18.10)' ) "max diff of gz = ", maxval(abs(gz-gz_out))
  write(*,'(A,ES18.10)' ) "max relative diff of pef= ", maxval(abs((pef(isc-1:iec+1,jsc-1:jec+1,:)- &
                                 pef_out(isc-1:iec+1,jsc-1:jec+1,:))/pef(isc-1:iec+1,jsc-1:jec+1,:)))

  write(*,*) ' elapsed time (secs) = ', total_time


contains

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



end program Main

