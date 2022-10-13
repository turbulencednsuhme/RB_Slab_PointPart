!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcDissipationNu.F90                          !
!    CONTAINS: subroutine CalcDissipationNu               !
!                                                         ! 
!    PURPOSE: Calculate the Nusselt number through the    !
!     global balance equations relating dissipation and   !
!     heat transport.                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcDissipationNu
      use mpih
      use param
      use local_arrays,only: vx,vy,vz,temp
      use mpi_param, only: kstart,kend
      use stat_arrays

      implicit none
      integer :: i,j,k
      integer :: imm,ipp,jmm,jpp,kmm,kpp,kp
      real :: dissip2,dissip1
      real :: udx3_m
      real :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real :: nu_th,nu_mu,volt
      real :: udx1,udx2,dissipte,my_dissipte,dissipth

      
      nu_th = 0.0d0
      nu_mu = 0.0d0

      udx1=dx1
      udx2=dx2
      
       do k=kstart,kend
        kp=k+1
        kpp=kpv(k)
        kmm=kmv(k)

       udx3_m=1.0/(zm(kpp)-zm(kmm))

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,i,jmm,jpp,imm,ipp)
!$OMP& PRIVATE(h11,h12,h13)
!$OMP& PRIVATE(h21,h22,h23)
!$OMP& PRIVATE(h31,h32,h33)
!$OMP& PRIVATE(dissipte,dissipth,dissip1,dissip2)
!$OMP& REDUCTION(+:nu_th, nu_mu)
       do j=1,n2m       
        jmm=jmv(j)
        jpp=jpv(j)

        dissip1 = 0.
        dissip2 = 0.
        
       do i=1,n1m
       imm= imv(i)
       ipp= ipv(i)

       h11=(vx(ipp,j,k)-vx(i,j,k))*udx1
       h12=(vx(i,jpp,k)-vx(i,j,k))*udx2
       h13=(vx(i,j,kpp)-vx(i,j,k))*udx3m(k)

       h21=(vy(ipp,j,k)-vy(i,j,k))*udx1
       h22=(vy(i,jpp,k)-vy(i,j,k))*udx2
       h23=(vy(i,j,kpp)-vy(i,j,k))*udx3m(k)

       h31=(vz(ipp,j,k)-vz(i,j,k))*udx1
       h32=(vz(i,jpp,k)-vz(i,j,k))*udx2
       h33=(vz(i,j,kp)-vz(i,j,k))*udx3c(k)

       dissipte = 2.0*(h11**2+h22**2+h33**2)+ &
               (h21+h12)**2+(h31+h13)**2+(h32+h23)**2


       nu_mu = nu_mu+dissipte*g3rm(k)/ray

       dissip1 = dissip1 + dissipte

       h31=(temp(ipp,j,k)-temp(imm,j,k))*udx1*0.5
       h32=(temp(i,jpp,k)-temp(i,jmm,k))*udx2*0.5
       h33=(temp(i,j,kpp)-temp(i,j,kmm))*udx3_m

       dissipth  = h31*h31 + h32*h32 + h33*h33 
       dissip2 = dissip2 + dissipth

       nu_th = nu_th+dissipth*g3rm(k)

       end do

       disste(j,k) =  disste(j,k) + dissip1 / ray / real(n1m)
       dissth(j,k) =  dissth(j,k) + dissip2 / real(n1m)

       end do
!$OMP  END PARALLEL DO
       end do



       call MpiSumRealScalar(nu_th)
       call MpiSumRealScalar(nu_mu)
      
      
       volt = 1.d0/(real(n3m)*real(n1m)*real(n2m))
      if(myid.eq.0) then
      nu_mu = nu_mu*volt
      nu_th = nu_th*volt 
      write(92,*) time,nu_mu,nu_th
      endif

      return   
      end
