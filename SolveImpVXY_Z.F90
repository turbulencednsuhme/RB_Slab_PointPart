!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpVXY_Z.F90                              !
!    CONTAINS: subroutine SolveImpVXY_Z                   !
!                                                         ! 
!    PURPOSE: Solve the implicit equation in the Z        !
!     direction and update velocities                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpVXY_Z(q)
      use param
      use local_arrays, only : rhs
      use mpi_param
      use mpih
      implicit none
      real, intent(inout) :: q(1:n1,1:n2,kstart-1:kend+1)
      real, dimension(n3) :: amkl,apkl,ackl,fkl
      integer :: jc,kc,info,ipkv(n3m),ic
      real :: betadx,ackl_b
      real :: amkT(n3m-1),ackT(n3m),apkT(n3m-1),appk(n3-3)
      real,allocatable :: rhst(:,:,:)

      allocate(rhst(1:n3,1:n1,jstart:jend))


!$OMP  PARALLEL
!$OMP  SINGLE
      call PackZ_UnpackR(rhs(:,:,kstart:kend),rhst(:,:,jstart:jend))
!$OMP  END SINGLE NOWAIT
!$OMP  SINGLE
       
      betadx=beta*al

      do kc=1,n3m
        ackl_b=1.0d0/(1.0d0-ac3sk(kc)*betadx)
        amkl(kc)=-am3sk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3sk(kc)*betadx*ackl_b
      enddo

      amkT=amkl(2:n3m)
      apkT=apkl(1:(n3m-1))
      ackT=ackl(1:n3m)

      call dgttrf(n3m,amkT,ackT,apkT,appk,ipkv,info)
!$OMP  END SINGLE
!$OMP  END PARALLEL

      do jc=jstart,jend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(ic,kc,ackl_b,fkl,info)
         do ic=1,n1m
          do kc=1,n3m
            ackl_b=1.0d0/(1.-ac3sk(kc)*betadx)
            fkl(kc)=rhst(kc,ic,jc)*ackl_b
          end do

          call dgttrs('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)

          do kc=1,n3m
            rhst(kc,ic,jc)=fkl(kc)
          end do
         enddo
!$OMP  END PARALLEL DO
      end do

      call PackR_UnpackZ(rhst(:,:,jstart:jend),rhs(:,:,kstart:kend))

      do kc=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(ic,jc)
      do jc=1,n2m
      do ic=1,n1m
      q(ic,jc,kc) = q(ic,jc,kc) + rhs(ic,jc,kc)
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      deallocate(rhst)

      return
      end
