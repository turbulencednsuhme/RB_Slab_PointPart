!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ReadInitCond.F90                               !
!    CONTAINS: subroutine ReadInitCond                    !
!                                                         ! 
!    PURPOSE: Read the initial conditions (temperature,   !
!     velocities)                                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ReadInitCond
      use mpih
      use mpi_param, only: kstart,kend,countk
      use local_arrays, only: temp,vy,vz,vx
      use param
      IMPLICIT NONE
      real ::  aaa
      character*70 :: filcnw2
      integer :: ihist
      integer :: n1o,n2o,n3o,n3om
      integer :: istro3,j,i,kstarto,kendo
      integer (kind=MPI_ADDRESS_KIND) :: extent,lb
      real :: stro3
      real :: intinfo(1:4)
      real, allocatable, dimension(:,:,:) :: tempold,vyold,vzold,vxold
      integer, parameter :: ghosts = 50 ! min(ghosts) = 1
      integer           :: kstartog,kendog
      integer, allocatable, dimension(:) :: countko

      
!EP   Reading old grid information by rank0
      if (myid .eq. 0) then
      filcnw2 = 'continua_grid.dat'
      open(13,file=filcnw2,status='unknown')
      rewind(13)                                                      
      read(13,*) n1o,n2o,n3o
      read(13,*) aaa,time
      read(13,*) istro3,stro3
      close(13)
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!EP   Bcasting old grid information and time
      call MPI_BCAST(n1o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n2o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n3o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(istro3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stro3,1,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(time,1,MDP,0,MPI_COMM_WORLD,ierr)
      
!EP   Check whether grid specifications have been updated
      if(n2o.ne.n2.or.n3o.ne.n3.or.n1o.ne.n1 &
       .or.istro3.ne.istr3.or.stro3.ne.str3) then
      if(myid.eq.0) write(*,*) "Interpolating new grid"
      if(n1.gt.n1o*2.or.n2.gt.n2o*2.or.n3.gt.n3o*2) then
      if(myid.eq.0) write(*,*) "New grid resolution can not be more",&
       " than twice the old resolution"
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif

      n3om = n3o - 1
      
      intinfo(1) = istro3
      intinfo(2) = stro3

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      allocate(countko(0:numtasks-1))

      call block(n3o-1, numtasks, myid, kstarto, kendo, countko) 

      kstartog = max(kstarto-ghosts,1)
      kendog   = min(kendo+ghosts,n3om)
      
      lb=0
      call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,extent,ierr)
!EP   temp
      allocate(tempold(1:n1o,1:n2o,kstartog-1:kendog+1))

      call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,4, &
       tempold(1:n1o,1:n2o,kstartog-1:kendog+1))

      if(myid.eq.numtasks-1) then
      do j=1,n2o
      do i=1,n1o
      tempold(i,j,n3o) = 0.0
      enddo
      enddo
      endif

      call interp(tempold,temp(1:n1,1:n2,kstart:kend),n1o,n2o,n3o,istro3,stro3,4,kstartog,kendog)

      deallocate(tempold)

!EP   vx
      allocate(vxold(1:n1o,1:n2o,kstartog-1:kendog+1))

      call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,1,vxold(1:n1o,1:n2o,kstartog-1:kendog+1))

      if(myid.eq.numtasks-1) then
      do j=1,n2o
      do i=1,n1o
      vxold(i,j,n3o) = 0.0
      enddo
      enddo
      endif

      call interp(vxold,vx(1:n1,1:n2,kstart:kend),n1o,n2o,n3o,istro3,stro3,1,kstartog,kendog)

      deallocate(vxold)

!EP   vy
      allocate(vyold(1:n1o,1:n2o,kstartog-1:kendog+1))

      call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,2,vyold(1:n1o,1:n2o,kstartog-1:kendog+1))

      if(myid.eq.numtasks-1) then
      do j=1,n2o
      do i=1,n1o
      vyold(i,j,n3o) = 0.0
      enddo
      enddo
      endif

      call interp(vyold,vy(1:n1,1:n2,kstart:kend),n1o,n2o,n3o,istro3,stro3,2,kstartog,kendog)

      deallocate(vyold)

!EP   vz
      allocate(vzold(1:n1o,1:n2o,kstartog-1:kendog+1))

      call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,3,vzold(1:n1o,1:n2o,kstartog-1:kendog+1))

      if(myid.eq.numtasks-1) then
      do j=1,n2o
      do i=1,n1o
      vzold(i,j,n3o) = 0.0
      enddo
      enddo
      endif

      call interp(vzold,vz(1:n1,1:n2,kstart:kend),n1o,n2o,n3o,istro3,stro3,3,kstartog,kendog)

      deallocate(vzold)

      else

!EP   One to one HDF read
      call mpi_read_continua(n1,n2,n3,kstart,kend,4,temp)
      call mpi_read_continua(n1,n2,n3,kstart,kend,1,vx)
      call mpi_read_continua(n1,n2,n3,kstart,kend,2,vy)
      call mpi_read_continua(n1,n2,n3,kstart,kend,3,vz)

      endif

      if(myid.eq.numtasks-1) then
      do j=1,n2
      do i=1,n1
      temp(i,j,n3) = 0.0
      vx(i,j,n3) = 0.0
      vy(i,j,n3) = 0.0
      vz(i,j,n3) = 0.0
      enddo
      enddo
      endif

      return                                                            
      end                                                               
                                                                       
