  subroutine add_noise( x, n, sigma, seed )
	
    USE mod_parallel_pdaf, &
      ONLY: COMM_filter, mype_filter, npes_filter, &
        MPI_INTEGER, MPI_STATUS_IGNORE, MPIerr
    
    implicit none 

    integer, intent(in)                             :: n
		real, intent(inout)   :: x(n)
    real, intent(in)                             :: sigma
    integer, intent(in)                             :: seed
    real                                         :: mean = 0.0
    integer(4), parameter                           :: buf_size = 1024
    real, dimension(buf_size)                    :: buf
    integer                                         :: i,j,k
    
    do i=0,npes_filter-1
      if ( i .eq. mype_filter ) then
        if ( mype_filter .gt. 0 ) then
          call mpi_recv( seed, 1, MPI_INTEGER, mype_filter-1, 0, &
            COMM_filter, mpi_status_ignore, MPIerr ) 
        endif

        do j=1,n
          if(modulo(j, buf_size) .eq. 1) then
            call r8vec_normal_ab(buf_size, mean, sigma, seed, buf )
          endif
          k = modulo(j-1, buf_size) + 1
          x(j) = x(j) + buf(k)
        end do

        if ( i .lt. (npes_filter-1) ) then
          call mpi_send( seed, 1, MPI_INTEGER, mype_filter+1, 0, &
            COMM_filter, MPIerr ) 
        endif
      end if
    end do


  end subroutine
