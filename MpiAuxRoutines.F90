!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: HdfRoutines.F90                                !
!    CONTAINS: subroutines MPI*                           !
!                                                         ! 
!    PURPOSE: Wrappers for MPI Routines                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine MpiBcastInt(n)
      use mpih
      implicit none
      integer, intent(inout) :: n
      
      call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      return
      end subroutine MpiBcastInt

!==============================================================================

      subroutine MpiBcastReal(n)
      use mpih
      implicit none
      real, intent(inout) :: n
      
      call MPI_BCAST(n,1,MDP,0,MPI_COMM_WORLD,ierr)

      return
      end subroutine MpiBcastReal
!==============================================================================

      subroutine MpiBcastReal1DArray(var,sx)
      use mpih
      implicit none
      integer, intent(in) :: sx
      real, dimension(1:sx), intent(inout) :: var
      
      call MPI_BCAST(var,sx,MDP,0,MPI_COMM_WORLD,ierr)

      return
      end subroutine MpiBcastReal1DArray
!==============================================================================

      subroutine MpiBcastReal2DArray(var,sx,sy)
      use mpih
      implicit none
      integer, intent(in) :: sx,sy
      real, dimension(1:sx,1:sy), intent(inout) :: var
      
      call MPI_BCAST(var,sx*sy,MDP,0,MPI_COMM_WORLD,ierr)

      return
      end subroutine MpiBcastReal2DArray
!==============================================================================

      subroutine MpiBarrier
      use mpih
      implicit none
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end subroutine MpiBarrier

!==============================================================================

      subroutine MpiSumRealScalar(var)
      use mpih
      implicit none
      real, intent(inout) :: var
      real :: buf
      
       call MPI_REDUCE(var,buf,1, &
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       var = buf

      return
      end subroutine MpiSumRealScalar
!==============================================================================

      subroutine MpiAllSumRealScalar(var)
      use mpih
      implicit none
      real, intent(inout) :: var
      real :: buf
      
       call MPI_ALLREDUCE(var,buf,1, &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

       var = buf

      return
      end subroutine MpiAllSumRealScalar
!==============================================================================

      subroutine MpiMaxRealScalar(var)
      use mpih
      implicit none
      real, intent(inout) :: var
      real :: buf
      
       call MPI_REDUCE(var,buf,1, &
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiMaxRealScalar
!==============================================================================

      subroutine MpiAllMaxRealScalar(var)
      use mpih
      implicit none
      real, intent(inout) :: var
      real :: buf
      
       call MPI_ALLREDUCE(var,buf,1, &
        MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiAllMaxRealScalar
!==============================================================================

      subroutine MpiMinRealScalar(var)
      use mpih
      implicit none
      real, intent(inout) :: var
      real :: buf
      
       call MPI_REDUCE(var,buf,1, &
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiMinRealScalar
!==============================================================================

      subroutine MpiSumReal1D(var,sz)
      use mpih
      implicit none
      integer, intent(in) :: sz
      real, intent(inout), dimension(1:sz) :: var
      real, dimension(1:sz) :: buf
      
       call MPI_REDUCE(var,buf,sz, &
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiSumReal1D
!==============================================================================
      subroutine MpiAllSumReal1D(var,sz)
      use mpih
      implicit none
      integer, intent(in) :: sz
      real, intent(inout), dimension(1:sz) :: var
      real, dimension(1:sz) :: buf
      
       call MPI_ALLREDUCE(var,buf,sz, &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiAllSumReal1D
!==============================================================================
      subroutine MpiAllSumReal2D(var,sx,sy)
      use mpih
      implicit none
      integer, intent(in) :: sx,sy
      real, intent(inout), dimension(1:sx,1:sy) :: var
      real, dimension(1:sx,1:sy) :: buf
      
       call MPI_ALLREDUCE(var,buf,sx*sy, &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiAllSumReal2D
!==============================================================================

      subroutine MpiAbort
      use mpih
      implicit none
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

      return
      end subroutine MpiAbort
!==============================================================================

      subroutine InitializeMPI
      use mpih
      use param, only: ismaster
      implicit none

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)

      if(myid.eq.0) ismaster = .true.

      return
      end subroutine InitializeMPI

!==============================================================================

      subroutine FinalizeMPI
      use mpih
      implicit none
      call MPI_FINALIZE(ierr)

      return
      end subroutine FinalizeMPI
!==============================================================================

      subroutine MpiAddGhosts(q,n1,n2)
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout) :: q(n1,n2,kstart-1:kend+1)
      real :: buf(n1,n2,1), buf2(n1,n2,1)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc

      mydata= n1*n2

      my_down= myid-1

      my_up= myid+1

      buf=0.0d0
      buf2=0.0d0

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL

      tag=1

      call MPI_ISEND(q(1,1,kstart-1),mydata,MDP, &
       my_down, tag, MPI_COMM_WORLD, req(1), ierr)

      call MPI_IRECV(buf(1,1,1), mydata, MDP, &
       my_up,tag, MPI_COMM_WORLD, req(2), ierr)

      call MPI_Waitall(2,req,status,ierr)

      call MPI_ISEND(q(1,1,kend),mydata,MDP, &
       my_down, tag, MPI_COMM_WORLD, req(1), ierr)

      call MPI_IRECV(buf2(1,1,1), mydata, MDP, &
       my_up,tag, MPI_COMM_WORLD, req(2), ierr)

      call MPI_Waitall(2,req,status,ierr)

      q(:,:,kend) = q(:,:,kend) + buf(:,:,1)

      q(:,:,kstart-1) = q(:,:,kstart-1) + buf2(:,:,1)

      end subroutine MpiAddGhosts
