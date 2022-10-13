
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SetTempBCs.F90                                 !
!    CONTAINS: subroutine SetTempBCs                      !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     temperature boundary conditions at the plates       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetTempBCs
      use param
      implicit none
      integer :: j,i,stripesx
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i,j)
      do j=1,n2
      do i=1,n1
              denbn(i,j)=0.d0
              denbs(i,j)=1.d0
      enddo
      enddo
!$OMP  END PARALLEL DO

      ! number of stripes (meaning a set of 1 adiabatic and 1 conduction stripe)
      stripesx = 2
      
        do j=1,n2m
         do i=1,n1m
      
          ! Bottom boundary is fully adiabatic
          lwtb(i,j)=1.0
      
          ! Stripes in x direction
          ! uwtb(i,j)=0.5 + sign(0.5,sin(2*stripesx*pi*(tm(i)/rext)))
          ! Top plate fully conducting 
          uwtb(i,j)=1.0
      
          enddo
         enddo


      return
      end
