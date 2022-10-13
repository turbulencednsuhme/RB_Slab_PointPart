!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateTemp.F90                      !
!    CONTAINS: subroutine ImplicitAndUpdateTemp           !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the temperature and call the implicit solver.       !
!     After this routine, the temperature has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateTemp
      use param
      use local_arrays, only: temp,hro,ruro,rhs
      use mpi_param, only: kstart,kend
      use ibm_param
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      integer :: n,ie,je,ke
      real    :: dq32,dq33,dcq3,dq31
      real    :: app,acc,amm
      real    :: alpec,udx1q,udx2q
      real    :: del1,del2, fcder
      real    :: usaldto,dne
 
      alpec=al/pec
      udx1q=dx1q
      udx2q=dx2q

      do kc=kstart,kend
      if( (kc.ge.2) .and. (kc.le.n3m-1) ) then
        km=kmv(kc)
        kp=kpv(kc)
        app=ap3ssk(kc)
        acc=ac3ssk(kc)
        amm=am3ssk(kc)
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,ip,im)
!$OMP& PRIVATE(dq31,dq32,dq33,dcq3)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
           im=imv(ic)
           ip=ipv(ic)

!  Conductive terms

           dq31=(temp(ip,jc,kc)-2.0*temp(ic,jc,kc)+temp(im,jc,kc))*udx1q
           dq32=(temp(ic,jp,kc)-2.0*temp(ic,jc,kc)+temp(ic,jm,kc))*udx2q
           dq33= temp(ic,jc,kp)*app+temp(ic,jc,kc)*acc+temp(ic,jc,km)*amm
           dcq3=dq32+dq33+dq31

           rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc) &
                    +alpec*dcq3)*dt

           ruro(ic,jc,kc)=hro(ic,jc,kc)
         enddo
        enddo
!$OMP  END PARALLEL DO
      endif

      if(kc.eq.1) then
      del1 = zm(1)-zc(1)
      del2 = zm(2)-zm(1)
      fcder = 2.d0/(del1*del2*(del1+del2))
        kp = kc + 1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,ip,im)
!$OMP& PRIVATE(dq31,dq32,dq33,dcq3)
        do jc=1,n2m
              jm=jmv(jc)
              jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
            dq31=(temp(ip,jc,kc)-2.0*temp(ic,jc,kc)+temp(im,jc,kc))*udx1q
            dq32=(temp(ic,jp,kc)-2.0*temp(ic,jc,kc)+temp(ic,jm,kc))*udx2q
            dq33=((temp(ic,jc,kp)-temp(ic,jc,kc))*del1 &
                 -(temp(ic,jc,kc)-denbs(ic,jc))*del2*lwtb(ic,jc)) &
                    *fcder
 
            dcq3=dq32+dq33+dq31

            rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc) &
                    +alpec*dcq3)*dt

            ruro(ic,jc,kc)=hro(ic,jc,kc)
            enddo
        enddo
!$OMP  END PARALLEL DO
      endif
!
!       UPPER COLD WALL
!     
      if(kc.eq.n3m) then
      del1 = zc(n3)-zm(n3m)
      del2 = zm(n3m)-zm(n3m-1)
      fcder = 2.d0/(del1*del2*(del1+del2))
        km = kc - 1 
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,ip,im)
!$OMP& PRIVATE(dq31,dq32,dq33,dcq3)
            do jc=1,n2m
              jm=jmv(jc)
              jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
            dq31=(temp(ip,jc,kc)-2.0*temp(ic,jc,kc)+temp(im,jc,kc))*udx1q
            dq32=(temp(ic,jp,kc)-2.0*temp(ic,jc,kc)+temp(ic,jm,kc))*udx2q
            dq33=((temp(ic,jc,km)-temp(ic,jc,kc))*del1 &
                 -(temp(ic,jc,kc)-denbn(ic,jc))*del2*uwtb(ic,jc)) &
                    *fcder
 
            dcq3=dq32+dq33+dq31

            rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc) &
                    +alpec*dcq3)*dt

            ruro(ic,jc,kc)=hro(ic,jc,kc)
            enddo
        enddo
!$OMP  END PARALLEL DO
      endif
      enddo

         usaldto = 1./aldto
         do n=1,npunfx(4)
           ic=indgeo(4,n,1)
           jc=indgeo(4,n,2)
           kc=indgeo(4,n,3)
           ie=indgeoe(4,n,1)
           je=indgeoe(4,n,2)
           ke=indgeoe(4,n,3)
           dne=((al*dt+aldto)*temp(ie,je,ke)-al*dt*dnbo(n))*usaldto
!
!      External boundaries
!
           rhs(ic,jc,kc) = -temp(ic,jc,kc) + dne*distb(4,n)  &
             +(1.-distb(4,n))*dnpr(n)
           dnbo(n)= temp(ie,je,ke)
         end do


      call SolveImpVXYZ_X(al*dt*0.5*dx1q/pec)

      call SolveImpVXYZ_Y(al*dt*0.5*dx2q/pec)

      call SolveImpTemp_Z

      return
      end
