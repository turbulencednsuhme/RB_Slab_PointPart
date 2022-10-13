!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpVZ_Z.F90                               !
!    CONTAINS: subroutine SolveImpVZ_Z                    !
!                                                         ! 
!    PURPOSE: Solve the implicit equation in the Z        !
!     direction and update velocities                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine SolveImpVZ_Z
      use param
      use local_arrays, only : vz,rhs
      use mpi_param
      implicit none
      real, dimension(n3) :: amkl,apkl,ackl, fkl
      real :: amkT(n3-1),apkT(n3-1)
      real :: appk(n3-2)
      real :: ackT(n3)
      integer :: jc,kc,info,ic
      integer :: ipkv(n3)
      real :: betadx,ackl_b
      real, allocatable, dimension(:,:,:) :: rhst

      allocate(rhst(1:n3,1:n1,jstart:jend))


!$OMP  PARALLEL
!$OMP  SINGLE
      call PackZ_UnpackR(rhs(:,:,kstart:kend),rhst(:,:,jstart:jend))
!$OMP  END SINGLE NOWAIT
!$OMP  SINGLE 

      betadx=beta*al

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(1)=1.d0
      do kc=2,n3m
        ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
        amkl(kc)=-am3ck(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ck(kc)*betadx*ackl_b
      enddo
      amkl(n3)=0.d0
      apkl(n3)=0.d0
      ackl(n3)=1.d0

      amkT=amkl(2:n3)
      apkT=apkl(1:(n3-1))
      ackT=ackl(1:n3)

      call dgttrf(n3,amkT,ackT,apkT,appk,ipkv,info)
!$OMP  END SINGLE
!$OMP  END PARALLEL

      do jc=jstart,jend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(ic,kc,fkl,ackl_b,info)
          do ic=1,n1m
            fkl(1)= 0.d0
          do kc=2,n3m
            ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
            fkl(kc) = rhst(kc,ic,jc)*ackl_b
          enddo
            fkl(n3)= 0.d0
          
          call dgttrs('N',n3,1,amkT,ackT,apkT,appk,ipkv,fkl,n3,info)
          do kc=1,n3
            rhst(kc,ic,jc)= fkl(kc)
          enddo
          enddo
!$OMP  END PARALLEL DO
      end do

      call PackR_UnpackZ(rhst(:,:,jstart:jend),rhs(:,:,kstart:kend)) 

      do kc=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic)
      do jc=1,n2m
      do ic=1,n1m
      vz(ic,jc,kc) = vz(ic,jc,kc) + rhs(ic,jc,kc)
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      if(allocated(rhst)) deallocate(rhst)

      if(kstart.eq.1) then
      vz(:,:,1)=0.0d0
      elseif(kend.eq.n3m) then
      vz(:,:,n3)=0.0d0
      endif

      return
      end
