!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity and temperature             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CreateInitialConditions
      use param
      use local_arrays, only: vy,vz,temp,vx
      use mpi_param, only: kstart,kend
      use mpih
      implicit none
      integer :: j,k,i
      real :: xxx,yyy,eps
      integer :: kc

      eps=0.01d0
      do k=kstart-1,kend+1
      kc=k
      if(kc.lt.1) then
      kc=1
      elseif(kc.gt.n3) then
      kc=n3
      endif
        do j=1,n2
           do i=1,n1
           vx(i,j,k) = 0.0d0

           yyy=zm(kc) 
           xxx=rc(j)            
          vy(i,j,k)=(2.0d0*yyy-6.0d0*yyy**2+4.0d0*yyy**3)*sin(3*xxx)*eps
          yyy=zc(kc)          
          xxx=rm(j)
          vz(i,j,k)=-yyy**2*(1.0d0-yyy)**2*cos(3.1*xxx)*eps
         enddo
        enddo
      enddo


       do k=kstart,kend
        do j=1,n2m    
           do i=1,n1m
!            temp(i,j,k)= denbs(i,j) - (denbs(i,j) - denbn(i,j))
!    %                   *zm(k)
             temp(i,j,k)=0.0d0
           enddo
          end do 
        end do
      return                                                            
      end                                                               
