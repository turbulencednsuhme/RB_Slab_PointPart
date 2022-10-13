!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVY.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVY             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the y (horizontal) dimension        !
!     and call the implicit solver                        !
!     After this routine, the velocity field in y has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateVY
      use param
      use mpih
      use local_arrays, only: vy,pr,rhs,dph,ru2
      use mpi_param, only: kstart,kend
      use ibm_param
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      integer :: n,ie,je,ke
      real    :: udx2,amm,app,acc
      real    :: dcvy,dpx22,usaldto,vye
      real    :: d22vy,d33vy,d11vy
      real    :: alre,udx1q,udx2q
      real    :: mck2, cksum2
      integer :: i,j,k

      
      alre=al/ren
      udx2=dx2*al
      udx1q=dx1q
      udx2q=dx2q

        do kc=kstart,kend
          km=kmv(kc)
          kp=kpv(kc)
          amm=am3sk(kc)
          acc=ac3sk(kc)
          app=ap3sk(kc)
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,ip,im)
!$OMP& PRIVATE(d11vy,d22vy,d33vy,dcvy,dpx22)
          do jc=1,n2m
           jm=jmv(jc)
           jp=jpv(jc)
            do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
!   viscid terms
            d11vy=(vy(ip,jc,kc)-2.0*vy(ic,jc,kc)+vy(im,jc,kc))*udx1q
            d22vy=(vy(ic,jp,kc)-2.0*vy(ic,jc,kc)+vy(ic,jm,kc))*udx2q
            d33vy=vy(ic,jc,kp)*app+vy(ic,jc,kc)*acc+vy(ic,jc,km)*amm

            dcvy=d22vy+d33vy+d11vy
 
            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2

            rhs(ic,jc,kc)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc) &
                          +alre*dcvy-dpx22)*dt

            ru2(ic,jc,kc)=dph(ic,jc,kc)
         enddo
       enddo
!$OMP  END PARALLEL DO
      enddo

         usaldto = 1./aldto
         do n=1,npunfx(2)
             ic=indgeo(2,n,1)
             jc=indgeo(2,n,2)
             kc=indgeo(2,n,3)
             ie=indgeoe(2,n,1)
             je=indgeoe(2,n,2)
             ke=indgeoe(2,n,3)
             vye=((al*dt+aldto)*vy(ie,je,ke)-al*dt*vybo(n))*usaldto
!
!     External boundaries
!
             rhs(ic,jc,kc) = -vy(ic,jc,kc) + vye*distb(2,n)  &
                  + (1-distb(2,n))*vypr(n)
             vybo(n)= vy(ie,je,ke)

         end do

      call SolveImpVXYZ_X(beta*al*dx1q)
                       
      call SolveImpVXYZ_Y(beta*al*dx2q)
      
      call SolveImpVXY_Z(vy)

      vy(:,n2,:) = vy(:,1,:)
     
      return
      end
