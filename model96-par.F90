module printout

  use mpi
  interface print_dbg
    module procedure print_dbg_int, print_dbg_real
  end interface

  contains

  subroutine print_dbg_int( str, val )
    character(len=*) :: str
    integer :: val
    integer :: mpi_rank
    call mpi_comm_rank(mpi_comm_world, mpi_rank, ierr)
    if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 4x, a, I5)") 'DEBUG', &
      str, val
    if(mpi_rank .eq. 0) call flush(6)
    call mpi_barrier(mpi_comm_world,ierr)
  end subroutine
  
  subroutine print_dbg_real( str, val )
    character(len=*) :: str
    real :: val
    integer :: mpi_rank
    call mpi_comm_rank(mpi_comm_world, mpi_rank, ierr)
    if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 4x, a, es12.4)") 'DEBUG', &
      str, val
    if(mpi_rank .eq. 0) call flush(6)
    call mpi_barrier(mpi_comm_world,ierr)
  end subroutine

end module

program lorenz96_seq

  use iso_c_binding
  use timer
  use printout
  use mpi
  USE parser, &           ! Parser function
    ONLY: parse
  implicit none

  INCLUDE 'melissa_api.f90'

  !include 'mpif.h'
  ! for distributing along the ranks if q = NG/mpi_size not integer
  integer, parameter                              :: MPI_MIN_BLK = 1

  integer, parameter                              :: NG = 1024
  integer, parameter                              :: NT = 10
  real, parameter                     :: F  = 0.2
  real, parameter                     :: dt = 0.01
  real, allocatable, dimension(:)     :: x
  real(kind=c_double), allocatable, dimension(:)     :: x_d
  real, allocatable, dimension(:)     :: x_old
  real, allocatable, dimension(:)     :: ki
  real, allocatable, dimension(:)     :: kj
  integer                                         :: i,j
  integer                                         :: ierr

  integer                                         :: mpi_rank
  integer                                         :: mpi_size
  integer                                         :: mpi_left, mpi_right
  
  integer                                         :: nl, nlt
  integer                                         :: nl_mod
  integer, allocatable, dimension(:)              :: nl_all
  integer(8)                                      :: nl_off
  integer                                         :: state_min_p
  integer                                         :: state_max_p

  integer                                         :: obs_block = 4
  integer                                         :: obs_last = 1
  integer                                         :: obs_share = 20
  real                                            :: obs_percent
  integer                                         :: dbg_var_int
  integer, parameter                              :: rnd_blk_size =32
  integer                                         :: member = 1
  CHARACTER(len=5)                                :: ensstr
  CHARACTER(len=5)                                :: memstr
  CHARACTER(len=5)                                :: mpestr
  CHARACTER(len=5)                                :: epostr
  integer(4)                                      :: epoch = 0
  real           :: mean = 0.0d0
  real           :: stdv = 0.2d0
  integer(4), parameter :: seed_init = 310780
  integer(4)        :: seed
  CHARACTER(len=32) :: handle  ! handle for command line parser
  logical :: generate_obs = .false.
  character(len=512) :: obs_file
  character(len=512) :: dataset_path

  integer :: timer_all = 0
  integer :: timer_iter = 1

  integer :: nsteps
  
  call timeit(timer_all, 'ini')
  call timeit(timer_iter, 'ini')

  call init_parallel()
  call mpi_barrier(mpi_comm_world, ierr)
  call timeit(timer_all, 'new')

  !   parse commandline args
  handle = 'member'             ! Control application of model error
  CALL parse(handle, member)
  handle = 'obs_gen'             ! Control application of model error
  CALL parse(handle, generate_obs)
  handle = 'seed'             ! Control application of model error
  CALL parse(handle, seed)
  handle = 'epoch'             ! Control application of model error
  CALL parse(handle, epoch)
  handle = 'obs_block'             ! Control application of model error
  CALL parse(handle, obs_block)
  handle = 'obs_share'             ! Control application of model error
  CALL parse(handle, obs_share)
  handle = 'obs_last'             ! Control application of model error
  CALL parse(handle, obs_last)
  obs_percent = real(obs_share)/100
  
  if (.not. generate_obs) then
    nsteps = 1
  else
    nsteps = obs_last
  end if

  if(mpi_rank .eq. 0) WRITE(*,"(a)") '-----------------------------------------------------'
  if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 4x, a)") 'INFO', &
    '==> model simulation started'
  if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 8x, a, I12)") 'INFO', &
    'dimension of state:      ', NG
  if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 8x, a, I12)") 'INFO', &
    'ensemble member:         ', member 
  if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 8x, a, I12)") 'INFO', &
    'number of observations:  ', int(obs_percent*NG) 
  if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 8x, a, L12)") 'INFO', &
    'generate observations:   ', generate_obs 
  if(mpi_rank .eq. 0) call flush(6)
  call mpi_barrier(mpi_comm_world,ierr)

  allocate( x(nlt) )   
  allocate( x_d(nl) )   
  allocate( x_old(nlt) )   
  allocate( ki(nlt) )
  allocate( kj(nlt) )

! create initial state with some noise and a 
! small perturbation
  if ( epoch .eq. 0 ) then
    x = 0.0
    if ( generate_obs ) then
      if (mpi_rank .eq. 0) then
        x(5) = 1.5
      endif
    else 
      if (mpi_rank .eq. 0) then
        x(5+member) = 1.0
      end if
      seed = seed_init
      call add_noise( x(3:3+nl-1), nl, 0.02, seed )
    end if
  else
    if ( generate_obs ) then
      call read_ens('true_state.txt')
    end if
  end if
  
  if ( .NOT. generate_obs ) then
    CALL MELISSA_INIT_F('lorenz_field', nl, 0, MPI_COMM_WORLD)
  end if 
  
  call print_dbg("rank ", mpi_rank)
  call print_dbg("size ", mpi_size)

  stepping: DO WHILE (nsteps > 0)
    
    do i = 1, nt
      ! use runga kutte RK4 method to solve lorenz96
      ! https://en.wikipedia.org/wiki/Runge-Kutta_methods 
      x_old   = x
      ki      = x
      call d96(ki, kj, F)
      x       = x + dt * kj/6.0
      ki      = x_old + dt * kj/2.0
      call exchange(ki)
      call d96(ki, kj, F)
      x       = x + dt * kj/3.0
      ki      = x_old + dt * kj/2.0 
      call exchange(ki)
      call d96(ki, kj, F)
      x       = x + dt * kj/3.0
      ki      = x_old + dt * kj 
      call exchange(ki)
      call d96(ki, kj, F)
      x       = x + dt * kj/6.0 
      call exchange(x)
    end do
    
    if ( .NOT. generate_obs ) then
      do i=1,nl
        x_d(i) = x(i)
      end do
      nsteps = melissa_expose('lorenz_field', x_d)
      do i=1,nl
        x(i) = x_d(i)
      end do
    else
      nsteps = nsteps - 1
    end if

    if (generate_obs) then
      call write_ens('true_state.txt')
!     simulate a small error in time direction
!     only if number of timesteps is large enough
      do i = 1,int(0.01*nt)
        ! use runga kutte RK4 method to solve lorenz96
        ! https://en.wikipedia.org/wiki/Runge-Kutta_methods 
        x_old   = x
        ki      = x
        call d96(ki, kj, F)
        x       = x + dt * kj/6.0
        ki      = x_old + dt * kj/2.0
        call exchange(ki)
        call d96(ki, kj, F)
        x       = x + dt * kj/3.0
        ki      = x_old + dt * kj/2.0 
        call exchange(ki)
        call d96(ki, kj, F)
        x       = x + dt * kj/3.0
        ki      = x_old + dt * kj 
        call exchange(ki)
        call d96(ki, kj, F)
        x       = x + dt * kj/6.0 
        call exchange(x)
      end do
      
      seed = seed_init
      call add_noise( x(3:3+nl-1), nl, 0.05, seed )
      
      write(ensstr, "(I5.5)") obs_last - nsteps
      
      call get_environment_variable( 'DATASET_PATH', dataset_path )
  
      obs_file = TRIM(dataset_path)//'/obs-'//TRIM(ensstr)//'.dat'
      
      call write_obs(obs_file, x(3:3+nl-1), obs_percent, obs_block)
    end if

  END DO stepping

  deallocate(nl_all)
  deallocate(x)
  deallocate(x_old)
  deallocate(ki)
  deallocate(kj)

  call mpi_barrier(mpi_comm_world, ierr)
  call timeit(timer_all, 'old')
 
  if(mpi_rank .eq. 0) WRITE(*,"(a)") '-----------------------------------------------------'
  if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 4x, a)") 'INFO', &
    '==> model simulation finished'
  if(mpi_rank .eq. 0) WRITE(*,"(2x, a, 8x, a, es12.4)") 'INFO', &
    'total time:              ', time_tot(timer_all) 
  if(mpi_rank .eq. 0) WRITE(*,"(a)") '-----------------------------------------------------'
  if(mpi_rank .eq. 0) call flush(6)
  call mpi_barrier(mpi_comm_world,ierr)

  !call timeit(timer_all, 'fin')
  
  500 call mpi_finalize(ierr)
  

contains
  subroutine d96(x, d, F)
    real, dimension(:), intent(IN)      :: x
    real, dimension(:), intent(OUT)     :: d
    real, intent(IN)                    :: F
    integer                                         :: N
    integer                                         :: i

    N = size(x)
    do i = 3,N-1
      d(i) = ( x(i+1) - x(i-2) ) * x(i-1) - x(i)
    end do
    d = d + F
  end subroutine

  subroutine add_noise( x, n, sigma, seed )
    real, intent(inout), dimension(:)   :: x
    integer, intent(in)                             :: n
    real, intent(in)                             :: sigma
    integer, intent(in)                             :: seed
    real                                         :: mean = 0.0
    integer(4), parameter                           :: buf_size = 1024
    real, dimension(buf_size)                    :: buf
    integer                                         :: i,j,k

    do i=0,mpi_size-1
      if ( i .eq. mpi_rank ) then
        if ( mpi_rank .gt. 0 ) then
          call mpi_recv( seed, 1, MPI_INTEGER, mpi_rank-1, 0, &
            mpi_comm_world, mpi_status_ignore, ierr ) 
        endif

        do j=1,n
          if(modulo(j, buf_size) .eq. 1) then
            call r8vec_normal_ab(buf_size, mean, sigma, seed, buf )
          endif
          k = modulo(j-1, buf_size) + 1
          x(j) = x(j) + buf(k)
        end do

        if ( i .lt. (mpi_size-1) ) then
          call mpi_send( seed, 1, MPI_INTEGER, mpi_rank+1, 0, &
            mpi_comm_world, ierr ) 
        endif
      end if
    end do


  end subroutine

  subroutine init_parallel()
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, mpi_size, ierr)

    allocate( nl_all(mpi_size) )

    ! middle ranks
    mpi_left    = mpi_rank - 1
    mpi_right   = mpi_rank + 1

    ! first and last rank
    if (mpi_rank == 0) then 
      mpi_left = mpi_size - 1
    elseif (mpi_rank == mpi_size-1) then
      mpi_right = 0
    endif

    nl_all = NG / mpi_size
    nl_mod = modulo( NG, mpi_size )
    do while (nl_mod .gt. 0)
      do i = 1, mpi_size
        if (nl_mod .gt. MPI_MIN_BLK) then
          nl_all(i) = nl_all(i) + MPI_MIN_BLK
          nl_mod = nl_mod - MPI_MIN_BLK
        else
          nl_all(i) = nl_all(i) + nl_mod
          nl_mod = 0
          exit
        endif
      end do
    end do

    nl = nl_all(mpi_rank+1)
    nlt = nl + 3
    
    nl_off = 0
    do i=1,mpi_rank
      nl_off = nl_off + nl_all(i)
    end do
    state_min_p = nl_off + 1
    state_max_p = nl_off + nl

  end subroutine

  subroutine exchange(x)
    real, dimension(:), intent(INOUT)      :: x

    if (modulo(mpi_rank,2) .eq. 0) then
      call mpi_send(x(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
      call mpi_recv(x(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
        MPI_STATUS_IGNORE, ierr)
      call mpi_send(x(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
      call mpi_recv(x(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
        MPI_STATUS_IGNORE, ierr)
    else
      call mpi_recv(x(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
        MPI_STATUS_IGNORE, ierr)
      call mpi_send(x(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
      call mpi_recv(x(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
        MPI_STATUS_IGNORE, ierr)
      call mpi_send(x(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
    endif
  end subroutine

  subroutine write_ens( file_name )
    integer     :: thefile
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer     :: ierr
    CHARACTER(len=*) :: file_name

    call mpi_file_open(MPI_COMM_WORLD, file_name, & 
      MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
      MPI_INFO_NULL, thefile, ierr)  

    disp = 0

    do i = 1, mpi_rank
      disp = disp + nl_all(i) * sizeof( x(3) )
    end do

    call mpi_file_set_view(thefile, disp, MPI_DOUBLE_PRECISION, & 
      MPI_DOUBLE_PRECISION, 'native', & 
      MPI_INFO_NULL, ierr)

    call mpi_file_write(thefile, x(3), nl, MPI_DOUBLE_PRECISION, & 
      MPI_STATUS_IGNORE, ierr)

    call mpi_file_close(thefile, ierr)
  end subroutine
  
  subroutine read_ens( file_name )
    integer     :: thefile
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer     :: ierr
    CHARACTER(len=*) :: file_name

    call mpi_file_open(MPI_COMM_WORLD, file_name, & 
      MPI_MODE_RDONLY, & 
      MPI_INFO_NULL, thefile, ierr)  

    disp = 0

    do i = 1, mpi_rank
      disp = disp + nl_all(i) * sizeof( x(3) )
    end do

    call mpi_file_set_view(thefile, disp, MPI_DOUBLE_PRECISION, & 
      MPI_DOUBLE_PRECISION, 'native', & 
      MPI_INFO_NULL, ierr)

    call mpi_file_read(thefile, x(3), nl, MPI_DOUBLE_PRECISION, & 
      MPI_STATUS_IGNORE, ierr)

    call mpi_file_close(thefile, ierr)
  end subroutine

  subroutine write_obs( file_name, x, share, blk_size )
    CHARACTER(len=*)  :: file_name
    real, dimension(:), intent(in) :: x
    real, intent(in)                :: share
    integer, intent(in)             :: blk_size
    integer                         :: dim_obs
    integer                         :: num_reg
    integer                         :: stride
    integer                         :: dim_obs_p
    real, allocatable, dimension(:) :: obs_p
    real, allocatable, dimension(:) :: obs_full_p
    integer                         :: offset
    integer                         :: index_tmp
    integer                         :: cnt_obs_p
    integer                         :: cnt_obs
    integer                         :: dim_obs_all(mpi_size)
    
    integer     :: thefile
    integer(kind=MPI_OFFSET_KIND) :: disp

    ! compute total number of observations
    dim_obs = share * NG
    if (dim_obs .eq. 0) then
      dim_obs = 1
    end if

    ! compute number of regions
    num_reg = dim_obs / blk_size
    if ( MODULO(dim_obs, blk_size) .ne. 0 ) then
      num_reg = num_reg + 1
    end if

    ! compute stride for regions
    stride = NG / num_reg

    ! determine number of obs in pe
    dim_obs_p = 0
    cnt_obs = 0
    do i=1,num_reg
      offset = (i-1) * stride
      do j=1,blk_size
        index_tmp = offset + j
        if ( (index_tmp .ge. state_min_p) .and. (index_tmp .le. state_max_p) ) then
          dim_obs_p = dim_obs_p + 1
        end if
        cnt_obs = cnt_obs + 1
        if ( cnt_obs .eq. dim_obs ) exit
        if ( index_tmp .eq. state_max_p ) exit
      end do
      if ( cnt_obs .eq. dim_obs ) exit
      if ( index_tmp .eq. state_max_p ) exit
    end do

    ALLOCATE( obs_p(dim_obs_p) )
    ALLOCATE( obs_full_p(nl) )

    obs_full_p = -999.0

    ! assign indices to index array
    cnt_obs_p = 0
    do i=1,num_reg
      offset = (i-1) * stride
      do j=1,blk_size
        index_tmp = offset + j
        
        if ( (index_tmp .ge. state_min_p) .and. (index_tmp .le. state_max_p) ) then
          cnt_obs_p = cnt_obs_p + 1
          obs_p(cnt_obs_p) = x(index_tmp - (state_min_p - 1))
          obs_full_p(index_tmp - (state_min_p - 1)) = obs_p(cnt_obs_p)
        end if
        if ( cnt_obs_p .eq. dim_obs_p ) exit
      end do
      if ( cnt_obs_p .eq. dim_obs_p ) exit
    end do
    
    ! write observations
    call mpi_allgather( dim_obs_p, 1, MPI_INTEGER, dim_obs_all, &
      1, MPI_INTEGER, mpi_comm_world, ierr )

    disp = 0

    do i = 1, mpi_rank
      disp = disp + dim_obs_all(i) * sizeof( obs_p(1) )
    end do

    call mpi_file_open(MPI_COMM_WORLD, file_name, & 
      MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
      MPI_INFO_NULL, thefile, ierr)  

    call mpi_file_set_view(thefile, disp, MPI_DOUBLE_PRECISION, & 
      MPI_DOUBLE_PRECISION, 'native', & 
      MPI_INFO_NULL, ierr)

    call mpi_file_write(thefile, obs_p, dim_obs_p, MPI_DOUBLE_PRECISION, & 
      MPI_STATUS_IGNORE, ierr)

    call mpi_file_close(thefile, ierr)
    
    DEALLOCATE( obs_p )
    DEALLOCATE( obs_full_p )
    
  end subroutine

end program lorenz96_seq
