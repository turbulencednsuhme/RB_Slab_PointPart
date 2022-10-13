!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectPressure.F90                            !
!    CONTAINS: subroutine CorrectPressure                 !
!                                                         ! 
!    PURPOSE: Apply the pressure correction to the        !
!     pressure                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CorrectPressure
      use param
      use local_arrays, only: pr,dph
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kp,km,jm,jp,jc,kc,ic,ip,im
      real    :: be,amm,acc,app

      be=al*beta
      do kc=kstart,kend
        kp=kpv(kc)
        km=kmv(kc)
        amm=amphk(kc)
        acc=acphk(kc)
        app=apphk(kc)
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,im,ip)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
              pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)-be*( &
              (dph(ip,jc,kc)-2.0*dph(ic,jc,kc)+dph(im,jc,kc))*dx1q+ &
              (dph(ic,jp,kc)-2.0*dph(ic,jc,kc)+dph(ic,jm,kc))*dx2q+ &
              (dph(ic,jc,kp)*app+dph(ic,jc,kc)*acc+dph(ic,jc,km)*amm))
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo
      return
      end
