!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: GlobalQuantities.F90                           !
!    CONTAINS: subroutine GlobalQuantities                !
!                                                         ! 
!    PURPOSE: Calculate maximum velocity and temperature, !
!     volume averaged Nusselt number                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine GlobalQuantities
      use param
      use local_arrays, only: vy,vz,vx,temp
      use mpi_param, only: kstart,kend
      use mpih
      implicit none
      integer :: jc,kc,kp,ic
      real :: vol,vzcen,fac2,anusin

      vmax=-100.d0
      anusin=0.d0 
      tempm=0.d0 

      vol = 1.d0/(alx3*dx3*real(n1m)*real(n2m))
        do kc=kstart,kend
        kp = kc + 1
        fac2 = g3rm(kc)
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic)
!$OMP& REDUCTION(max: vmax1,vmax2,vmax3)
!$OMP& REDUCTION(+: anusin,tempm)
          do jc=1,n2m
            do ic=1,n1m
              vmax(1) = max(vmax(1),abs(vx(ic,jc,kc)))
              vmax(2) = max(vmax(2),abs(vy(ic,jc,kc)))
              vmax(3) = max(vmax(3),abs(vz(ic,jc,kc)))
              vzcen = (vz(ic,jc,kc)+vz(ic,jc,kp))*0.5d0
              anusin=anusin+temp(ic,jc,kc)*vzcen*fac2
             tempm=tempm+temp(ic,jc,kc)*fac2
       enddo
       enddo
!$OMP  END PARALLEL DO
       enddo



      call MpiSumRealScalar(tempm)
      call MpiSumRealScalar(anusin)
      call MpiMaxRealScalar(vmax(1))
      call MpiMaxRealScalar(vmax(2))
      call MpiMaxRealScalar(vmax(3))

      if(myid.eq.0) then
      anusin=1.d0 + dsqrt(pra*ray)*anusin*vol
      tempm=tempm*vol
      write(95,*) time, anusin, tempm
      endif

      return   
      end     
