!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectVelocity.F90                            !
!    CONTAINS: subroutine CorrectVelocity                 !
!                                                         ! 
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CorrectVelocity
      use param
      use local_arrays, only: vy,vz,dph,vx
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,jm,kc,km,ic,im
      real    :: usukm,udx2,udx1,locdph

      udx1 = al*dt*dx1
      udx2 = al*dt*dx2
      do kc=kstart,kend
        km=kmv(kc)
        usukm = al*dt*udx3c(kc)
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,ic,im,locdph)
        do jc=1,n2m
          jm=jmv(jc)
          do ic=1,n1m
          im=imv(ic)
          locdph=dph(ic,jc,kc)
          vx(ic,jc,kc)=vx(ic,jc,kc)-(locdph-dph(im,jc,kc))*udx1
          vy(ic,jc,kc)=vy(ic,jc,kc)-(locdph-dph(ic,jm,kc))*udx2
          vz(ic,jc,kc)=vz(ic,jc,kc)-(locdph-dph(ic,jc,km))*usukm
        enddo 
       enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end
