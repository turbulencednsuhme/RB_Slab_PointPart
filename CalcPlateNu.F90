!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcPlateNu.F90                                !
!    CONTAINS: subroutine CalcPlateNu                     !
!                                                         ! 
!    PURPOSE: Calculate the Nusselt number at the top     !
!     and bottom plates and output to a file.             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcPlateNu
      use param
      use local_arrays, only: temp
      use mpih
      implicit none
      integer :: j,i
      real :: grtlow,grtupp,fcder,fcdern
      real ::  my_anussupp,surf,fac2
      real :: del1,del,deln,del1n,udel1q,udelq,udel1qn,udelqn
  
!
!     COMPUTATION OF THE NUSSELT NUMBER AT THE 
!     LOWER AND UPPER WALLS
!

       if(myid.eq.0.OR.myid.eq.numtasks-1) then
         if(myid.eq.0) then
      anusslow = 0.d0
      del1 = zm(1)-zc(1)
      del  = zm(2)-zc(1)
      udel1q = 1.d0/del1**2
      udelq = 1.d0/del**2
      fcder = 1.d0/(1.d0/del1 - 1.d0/del)
      endif
         if(myid.eq.numtasks-1) then
      anussupp = 0.d0
      my_anussupp = 0.d0
      del1n = -zm(n3m)+zc(n3)
      deln  = -zm(n3m-1)+zc(n3)
      udel1qn = 1.d0/del1n**2
      udelqn = 1.d0/deln**2
      fcdern = 1.d0/(1.d0/del1n - 1.d0/deln)
        endif
      surf = rext*rext2
      fac2 = 1/(dx1*dx2)

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,i,grtlow,grtupp)
!$OMP& REDUCTION(+:anusslow,anussupp,my_anussupp)
      do j=1,n2m
         do i=1,n1m
         if(myid.eq.0) then
         grtlow = (temp(i,j,1)*udel1q-temp(i,j,2)*udelq &
                    -denbs(i,j)*(udel1q-udelq)) *fcder
           anusslow = anusslow + grtlow*fac2
      endif
         if(myid.eq.numtasks-1) then
           grtupp = (temp(i,j,n3m)*udel1qn-temp(i,j,n3m-1)*udelqn &
                    -denbn(i,j)*(udel1qn-udelqn) )*fcdern
           anussupp = anussupp + grtupp*fac2
           my_anussupp = my_anussupp + grtupp*fac2
         endif
        enddo
      end do
!$OMP  END PARALLEL DO

      if(myid.eq.0) then
      anusslow = -1.0*anusslow / surf
      call MPI_IRECV(anussupp,1,MDP,numtasks-1,1,MPI_COMM_WORLD,status &
       ,ierr)
      endif

      if(myid.eq.numtasks-1) then
      my_anussupp = my_anussupp / surf
      call MPI_SEND(my_anussupp,1,MDP,0,1,MPI_COMM_WORLD,ierr)
      endif
      endif

      if(myid.eq.0) then
       write(97,546) time, anusslow, anussupp
 546   format(4(1x,e14.6))
      endif

      return         
      end                                                               
