!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsTemp.F90                          !
!    CONTAINS: subroutine ExplicitTermsTemp               !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the temperature.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ExplicitTermsTemp
      use param
      use local_arrays, only: vy,vz,hro,temp,vx
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic
      integer :: kpp,km,kp,jp,jmm,jpp,ip,imm,ipp
      real    :: h32,h33,udx2,udx1,h31

      udx1=dx1*0.5
      udx2=dx2*0.5
      do kc=kstart,kend
      km=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jp,jmm,jpp,ic,ip,ipp,imm)
!$OMP& PRIVATE(h31,h32,h33)
      do jc=1,n2m
      jp=jc+1
      jmm=jmv(jc)
      jpp=jpv(jc)
      do ic=1,n1m
      ip=ic+1
      ipp=ipv(ic)
      imm=imv(ic)
!
!
!    rho vx term
!
!
!                d  rho q_t 
!             -----------
!                d   t      
!
      h31=(vx(ip,jc,kc)*(temp(ipp,jc,kc)+temp(ic,jc,kc))- &
           vx(ic,jc,kc)*(temp(ic,jc,kc)+temp(imm,jc,kc)) &
          )*udx1 
!
!
!    rho vy term
!
!
!                d  rho q_r 
!             -----------
!                d   r      
!
      h32=(vy(ic,jp,kc)*(temp(ic,jpp,kc)+temp(ic,jc,kc))- &
           vy(ic,jc,kc)*(temp(ic,jc,kc)+temp(ic,jmm,kc)) &
          )*udx2
!
!    rho vz term
!
!
!                 d  rho q_x 
!                -----------
!                 d   x      
!
      h33=(vz(ic,jc,kp)*(temp(ic,jc,kpp)+temp(ic,jc,kc))- &
           vz(ic,jc,kc)*(temp(ic,jc,kc)+temp(ic,jc,km)) &
          )*udx3m(kc)*0.5d0

      hro(ic,jc,kc)=-(h31+h32+h33)
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end
