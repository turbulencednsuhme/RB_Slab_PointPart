!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpVXYZ_X.F90                             !
!    CONTAINS: subroutine SolveImpVXYZ_X                  !
!                                                         ! 
!    PURPOSE: Solve the implicit equation in the X        !
!     direction                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpVXYZ_X(betadx)
!EP   Solves tridiagonal system in i direction
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic
      real,intent(in) :: betadx
      real, dimension(n1):: amil,apil,acil,fil
      real :: ackl_b

      ackl_b = 1.0/(1.0+2.0*betadx)
      do ic=1,n1m
       apil(ic)=-betadx*ackl_b
       acil(ic)=1.0d0
       amil(ic)=-betadx*ackl_b
      end do

      do kc=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic,fil)
          do jc=1,n2m
             do ic=1,n1m
                fil(ic)=rhs(ic,jc,kc)*ackl_b
             enddo
                call tripvmy_line(amil,acil,apil,fil,1,n1m,n1)
             do ic=1,n1m
                rhs(ic,jc,kc) = fil(ic)  
             enddo
          end do
!$OMP  END PARALLEL DO
      end do 

      return
      end
