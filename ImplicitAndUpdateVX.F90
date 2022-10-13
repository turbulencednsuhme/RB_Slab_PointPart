!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVX.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVX             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the X (vertical) direction and call !
!     the implicit solver. After this routine, the        !
!     vertical velocity has been updated to the new       !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateVX
      use param
      use local_arrays, only: pr,rhs,ru1,vx,dq
      use mpi_param, only: kstart,kend
      use ibm_param
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real    :: udx1,amm,acc,app
      real    :: dcvx,dpx11,usaldto
      real    :: d22vx,d33vx,d11vx
      real    :: alre,udx1q,udx2q,vxe
      integer :: n,ie,je,ke

      alre=al/ren

      udx1=dx1*al
      udx1q=dx1q
      udx2q=dx2q
!
!  compute the rhs of the factored equation
!  everything at i,j+1/2,k+1/2
!
!    points inside the flowfield
!
        do kc=kstart,kend
          km=kmv(kc)
          kp=kpv(kc)
          amm=am3sk(kc)
          acc=ac3sk(kc)
          app=ap3sk(kc)
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,ip,im)
!$OMP& PRIVATE(d11vx,d22vx,d33vx,dcvx,dpx11)
          do jc=1,n2m
           jm=jmv(jc)
           jp=jpv(jc)
            do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)

!
!   Viscid Terms
!
            d11vx=(vx(ip,jc,kc)-2.0*vx(ic,jc,kc)+vx(im,jc,kc))*udx1q
            d22vx=(vx(ic,jp,kc)-2.0*vx(ic,jc,kc)+vx(ic,jm,kc))*udx2q
            d33vx=vx(ic,jc,kp)*app+vx(ic,jc,kc)*acc+vx(ic,jc,km)*amm

            dcvx=d11vx+d22vx+d33vx


!   component of grad(pr) along direction

            dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1

            rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc) &
                          +alre*dcvx-dpx11)*dt

            ru1(ic,jc,kc)=dq(ic,jc,kc)
         enddo
       enddo
!$OMP  END PARALLEL DO
      enddo

       usaldto = 1./aldto
!
!      External boundaries
!
        do n=1,npunfx(1)
         ic=indgeo(1,n,1)
         jc=indgeo(1,n,2)
         kc=indgeo(1,n,3)
         ie=indgeoe(1,n,1)
         je=indgeoe(1,n,2)
         ke=indgeoe(1,n,3)
         vxe=((al*dt+aldto)*vx(ie,je,ke)-al*dt*vxbo(n))*usaldto
         rhs(ic,jc,kc) = -vx(ic,jc,kc) + vxe*distb(1,n)  &
                  + (1-distb(1,n))*vxpr(n)
         vxbo(n)= vx(ie,je,ke)
        end do

        do n=1,npunifx(1)
         ic=indgeoee(1,n,1)
         jc=indgeoee(1,n,2)
         kc=indgeoee(1,n,3)
         rhs(ic,jc,kc) = -vx(ic,jc,kc)
        end do

      call SolveImpVXYZ_X(beta*al*dx1q)

      call SolveImpVXYZ_Y(beta*al*dx2q)
      
      call SolveImpVXY_Z(vx)
      
      vx(n1,:,:) = vx(1,:,:)
     
      return
      end
