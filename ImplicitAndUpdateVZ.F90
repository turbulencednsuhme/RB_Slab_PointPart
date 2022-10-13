!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVZ.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVZ             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (vertical) dimension          !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in z has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateVZ
      use param
      use local_arrays, only: vz,qcap, pr,ru3,rhs
      use mpi_param, only: kstart,kend
      use ibm_param
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      integer :: n,ie,je,ke
      real    :: udx3
      real    :: dvz2,dvz3,dcvz,dpx33,dvz1
      real    :: app,acc,amm
      real    :: alre,udx1q,udx2q,usaldto,vze
      integer :: kstartp
      if(kstart.eq.1) then
      kstartp=2
      else
      kstartp=kstart
      endif
        
      alre=al/ren
      udx1q=dx1q
      udx2q=dx2q

      do kc=kstartp,kend
        km=kmv(kc)
        kp=kc+1
        udx3 = al*udx3c(kc)
        amm=am3ck(kc)
        acc=ac3ck(kc)
        app=ap3ck(kc)
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,ip,im)
!$OMP& PRIVATE(dvz1,dvz2,dvz3,dcvz,dpx33)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
            do ic=1,n1m
              im=imv(ic)
              ip=ipv(ic)
 
!  viscid terms
 
            dvz1=(vz(im,jc,kc)-2.0*vz(ic,jc,kc)+vz(ip,jc,kc))*udx1q
            dvz2=(vz(ic,jm,kc)-2.0*vz(ic,jc,kc)+vz(ic,jp,kc))*udx2q
            dvz3=vz(ic,jc,kp)*app+vz(ic,jc,kc)*acc+vz(ic,jc,km)*amm
            dcvz=dvz2+dvz3+dvz1
 
!  pressure gradient
 
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3

            rhs(ic,jc,kc)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc) &
                          +alre*dcvz-dpx33)*dt 

            ru3(ic,jc,kc)=qcap(ic,jc,kc)
         enddo
       enddo
!$OMP  END PARALLEL DO
      enddo


         usaldto = 1./aldto
         do n=1,npunfx(3)
           ic=indgeo(3,n,1)
           jc=indgeo(3,n,2)
           kc=indgeo(3,n,3)
           ie=indgeoe(3,n,1)
           je=indgeoe(3,n,2)
           ke=indgeoe(3,n,3)
           vze=((al*dt+aldto)*vz(ie,je,ke)-al*dt*vzbo(n))*usaldto
!
!      External boundaries
!
           rhs(ic,jc,kc) = -vz(ic,jc,kc) + vze*distb(3,n)  &
                  + (1-distb(3,n))*vzpr(n)
           vzbo(n)= vz(ie,je,ke)
         end do

        do n=1,npunifx(3)
         ic=indgeoee(3,n,1)
         jc=indgeoee(3,n,2)
         kc=indgeoee(3,n,3)
         rhs(ic,jc,kc) = -vz(ic,jc,kc)
        end do


      call SolveImpVXYZ_X(beta*al*dx1q)
      call SolveImpVXYZ_Y(beta*al*dx2q)
      call SolveImpVZ_Z

      return
      end
