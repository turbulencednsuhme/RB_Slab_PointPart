      subroutine PPForceVel
      USE param
      USE mpih
      USE mpi_param
      use pointparticle
      USE local_arrays

      IMPLICIT NONE
      integer :: ic,jc,kc,im,jm,km

      call MpiAddLowerGhost(for_xc_part,n1,n2)
      call MpiAddLowerGhost(for_yc_part,n1,n2)
      call MpiAddLowerGhost(for_zc_part,n1,n2)
      call MpiAddLowerGhost(for_tc_part,n1,n2)
      call MpiAddUpperGhost(for_xc_part,n1,n2)
      call MpiAddUpperGhost(for_yc_part,n1,n2)
      call MpiAddUpperGhost(for_zc_part,n1,n2)
      call MpiAddUpperGhost(for_tc_part,n1,n2)

      call MpiBarrier

      do kc=kstart,kend
       km=kc-1
       do jc=1,n2m
       jm=jmv(jc)
        do ic=1,n1m
         im=imv(ic)
         vx(ic,jc,kc)=vx(ic,jc,kc)+(for_xc_part(ic,jc,kc)+for_xc_part(im,jc,kc)) &
          *0.5d0
         vy(ic,jc,kc)=vy(ic,jc,kc)+(for_yc_part(ic,jc,kc)+for_yc_part(ic,jm,kc)) &
          *0.5d0
         vz(ic,jc,kc)=vz(ic,jc,kc)+(for_zc_part(ic,jc,kc)+for_zc_part(ic,jc,km)) &
          *0.5d0
         temp(ic,jc,kc)=temp(ic,jc,kc)+for_tc_part(ic,jc,kc)
        end do
       end do
      end do


      return
      end
