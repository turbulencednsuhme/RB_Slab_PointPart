!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpTemp_Z.F90                             !
!    CONTAINS: subroutine SolveImpTemp_Z                  !
!                                                         ! 
!    PURPOSE: Solve the implicit equation in the Z        !
!     direction and update temperature                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpTemp_Z
      use param
      use local_arrays, only : temp,rhs
      use mpi_param
      implicit none
      real, dimension(n3) :: amkl,apkl,ackl, fkl
      integer :: jc,kc,info,ipkv(n3m),ic
      real :: betadx,ackl_b,ackl_a,a33m,a33p
      real :: amkT(n3m-1),ackT(n3m),apkT(n3m-1),appk(n3-3)
      real, allocatable, dimension(:,:,:) :: rhst

      allocate(rhst(1:n3,1:n1,jstart:jend))

!$OMP  PARALLEL
!$OMP  SINGLE
      call PackZ_UnpackR(rhs(:,:,kstart:kend),rhst(:,:,jstart:jend))
!$OMP  END SINGLE NOWAIT
!$OMP  SINGLE

      betadx=0.5d0*al*dt/pec
     

      do kc=2,n3m-1
        ackl_b=1.0d0/(1.-ac3ssk(kc)*betadx)
        amkl(kc)=-am3ssk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ssk(kc)*betadx*ackl_b
      enddo

      do jc=jstart,jend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(ic,kc,ackl_b,fkl,info)
         do ic=1,n1m

          a33m=dx3q/(g3rm(1)*g3rc(1))
          ackl_a=-(ap3ssk(1)+lwtb(ic,jc)*a33m*2.0d0)
          ackl_b=1.0d0/(1.-ackl_a*betadx)
          ackl(1)=1.0d0
          apkl(1)=-ap3ssk(1)*betadx*ackl_b
          fkl(1)=rhst(1,ic,jc)*ackl_b

          do kc=2,n3m-1
            ackl_b=1./(1.-ac3ssk(kc)*betadx)
            fkl(kc)=rhst(kc,ic,jc)*ackl_b
          end do

          a33p=dx3q/(g3rm(n3m)*g3rc(n3))
          ackl_a = -(am3ssk(n3m)+uwtb(ic,jc)*a33p*2.0d0)
          ackl_b = 1.0d0/(1.-betadx*ackl_a)
          amkl(n3m)=-am3ssk(n3m)*ackl_b*betadx
          ackl(n3m)=1.0d0
          fkl(n3m)=rhst(n3m,ic,jc)*ackl_b

          amkT=amkl(2:n3m)
          apkT=apkl(1:(n3m-1))
          ackT=ackl(1:n3m)

          call dgttrf(n3m,amkT,ackT,apkT,appk,ipkv,info)
          
          call dgttrs('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)
          
          do kc=1,n3m
            rhst(kc,ic,jc)= fkl(kc)
          end do
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
      temp(ic,jc,kc) = temp(ic,jc,kc) + rhs(ic,jc,kc)
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      if(allocated(rhst)) deallocate(rhst)

      return
      end
