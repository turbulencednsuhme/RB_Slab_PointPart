!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcLocalDivergence.F90                        !
!    CONTAINS: subroutine CalcLocalDivergence             !
!                                                         ! 
!    PURPOSE: Compute the divergence of the intermediate  !
!     velocity at every point for the pressure            !
!     correction step                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcLocalDivergence
      use param
      use local_arrays, only: vx,vy,vz,dph
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,jp,kc,kp,ic,ip
      real    :: usdtal,dqcap   

      usdtal = 1.d0/(dt*al)

      do kc=kstart,kend
        kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic,jp,ip,dqcap)
        do jc=1,n2m
          jp=jpv(jc)
            do ic=1,n1m
              ip=ipv(ic)
              dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
                    +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
                    +(vz(ic,jc,kp)-vz(ic,jc,kc))*udx3m(kc)
              dph(ic,jc,kc)=dqcap*usdtal
            enddo
         enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end
