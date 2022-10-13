      program papero
      use mpih
      use mpi_param
      use param
      implicit none
      character(len=4) :: dummy
      integer :: n,ns,nt
      real :: ts,te

      call InitializeMPI

      nt = 0
!$OMP  PARALLEL
!$OMP& REDUCTION(+:nt)
      nt = nt + 1
!$OMP  END PARALLEL

      if (ismaster) then 
       write(6,*) 'MPI tasks=', numtasks
       write(6,*) 'OMP threads per task=', nt
       write(6,*) 'No. of processors=', nt*numtasks
      end if

      call ReadInputFile
      call ReadPPartInput
      call MpiBarrier

      if(ismaster) call OpenLogs

      if(ismaster) then

      write(6,112)rext/alx3
  112 format(//,20x,'R A Y B E N ',//,10x, &
       '2D Cell with aspect-ratio:  D/H = ',f5.2)

      write(6,142) 
  142 format(//,8x,'Periodic lateral wall boundary condition')

      write(6,202) ray,pra
  202 format(/,5x,'Parameters: ',' Ra=',e10.3,' Pr= ',e10.3)

       if(idtv.eq.1) then
         write(6,204) cflmax
  204 format(/,5x,'Variable dt and fixed cfl= ', &
       e11.4,/ )            
       else 
         write(6,205) dtmax,cfllim
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=', &
        e11.4,/ )            
       endif

      endif

      call InitTimeMarchScheme

      call mpi_workdistribution

      call InitArrays

!     the solution of the problem starts

      ts=MPI_WTIME()

      call gcurv
      
      te=MPI_WTIME()
      
      if(myid.eq.0) then
        open(27,file="Total_time.out")
        write(27,*)"Total simulation time in sec.: ", te-ts
        close(27)
      endif

      call CloseLogs

      call DeallocateArrays
      
      call MPI_FINALIZE(ierr)

      stop                                                              
      end                                                               
