!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVZ.F90                            !
!    CONTAINS: subroutine ExplicitTermsVZ                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the z (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ExplicitTermsVZ
      use param
      use local_arrays, only: vy,vz,qcap,temp,vx
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc
      integer :: km,kp,jp,jmm,jpp,ic,im,ip
      real    :: h32,h33,h31
      real    :: tempit,udx1,udx2
      integer :: kstartp

      if(kstart.eq.1) then
      kstartp=2
      else
      kstartp=kstart
      endif

      udx1=dx1*0.25
      udx2=dx2*0.25

      do kc=kstartp,kend
      km=kmv(kc)
      kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jmm,jpp,jp,im,ip)
!$OMP& PRIVATE(h31,h32,h33,tempit)
      do jc=1,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      jp=jc+1
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
!
!
!    vz vx term
!
!
!                d  q_x q_t 
!             -----------
!                d   t      
!
!
      h31=(((vx(ip,jc,kc)+vx(ip,jc,km)) &
           *(vz(ip,jc,kc)+vz(ic,jc,kc))) &
          -((vx(ic,jc,kc)+vx(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(im,jc,kc))))*udx1
!
!    vz vy term
!
!
!                d  q_x q_r 
!             -----------
!                d   r      
!
      h32=(((vy(ic,jp,kc)+vy(ic,jp,km)) &
           *(vz(ic,jpp,kc)+vz(ic,jc,kc))) &
          -((vy(ic,jc,kc)+vy(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(ic,jmm,kc))))*udx2
!
!    vz vz term
!
!
!                 d  q_x q_x 
!                -----------
!                 d   x      
!
      h33=((vz(ic,jc,kp)+vz(ic,jc,kc))*(vz(ic,jc,kp)+vz(ic,jc,kc)) &
          -(vz(ic,jc,kc)+vz(ic,jc,km))*(vz(ic,jc,kc)+vz(ic,jc,km)) &
          )*udx3c(kc)*0.25d0
 
!  add the buoyancy term
 
        tempit=(temp(ic,jc,kc)+temp(ic,jc,km))*0.5d0
 
        qcap(ic,jc,kc)=-(h31+h32+h33)+tempit
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end
