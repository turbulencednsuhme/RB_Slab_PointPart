!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpVXYZ_Y.F90                             !
!    CONTAINS: subroutine SolveImpVXYZ_Y                  !
!                                                         ! 
!    PURPOSE: Solve the implicit equation in the Y        !
!     direction                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine SolveImpVXYZ_Y(betadx)
!EP   Solves tridiagonal system in j direction
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic
      real,intent(in) :: betadx
      real, dimension(n2):: amjl,apjl,acjl,fjl
      real :: ackl_b

!     betadx=beta*al

      ackl_b = 1.0/(1.0+2.0*betadx)
      do jc=1,n2m
       apjl(jc)=-betadx*ackl_b
       acjl(jc)=1.0d0
       amjl(jc)=-betadx*ackl_b
      end do

      do kc=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic,fjl)
          do ic=1,n1m
             do jc=1,n2m
                fjl(jc)=rhs(ic,jc,kc)*ackl_b
             enddo
                call tripvmy_line(amjl,acjl,apjl,fjl,1,n2m,n2)
             do jc=1,n2m
                rhs(ic,jc,kc) = fjl(jc)  
             enddo
          end do
!$OMP  END PARALLEL DO
      end do 

      return
      end
