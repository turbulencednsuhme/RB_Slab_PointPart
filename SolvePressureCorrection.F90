!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolvePressureCorrection.F90                    !
!    CONTAINS: subroutine SolvePressureCorrection,        !
!     CreateFFTTmpArrays, DestroyFFTTmpArrays             ! 
!                                                         ! 
!    PURPOSE: Compute the pressure correction by solving  !
!     a Poisson equation                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolvePressureCorrection
      use param
      use mpih
      use local_arrays, only: dph
      use mpi_param
      implicit none
      integer :: j,k,info,i
      real :: acphT_b
      real :: xr(n2m,n1m)
      complex :: xa(n2mh,n1m)
      real :: appph(n3m-2)
      real :: amphT(n3m-1), apphT(n3m-1)
      real, dimension(n3m) :: acphT,drhs,apph,amph
      integer :: phpiv(n3m)
      integer :: jmh
      real,allocatable,dimension(:,:,:) :: dpht

      allocate(dpht(1:n3,1:n1,jstartp:jendp))


      do k=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,i)
        do j=1,n2m
          do i=1,n1m
           xr(j,i)=dph(i,j,k)
          enddo
        enddo
!$OMP  END PARALLEL DO
        
        call dfftw_execute_dft_r2c(fwd_plan,xr,xa)

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,i)
        do j=1,n2mh
         do i=1,n1m
         dph(i,j,k)=real(xa(j,i)/(n1m*n2m))
         dph(i,j+n2mh,k)=aimag(xa(j,i)/(n1m*n2m))
        enddo
        enddo
!$OMP  END PARALLEL DO
      end do


      call PackZ_UnpackRP(dph(:,:,kstart:kend),dpht(:,:,jstartp:jendp))




!==================================================================
!     inversion of the matrix in the x2 and x3 directions (real part)
 
      do j=jstartp,jendp
        jmh=jmhv(j)
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i,k,amph,apph,drhs,info,amphT,apphT,acphT)
!$OMP& PRIVATE(acphT_b,appph,phpiv)
        do i=1,n1m
         do k = 1,n3m
          acphT_b=1.0/(acphk(k)-ak2(jmh)-ak1(i))
          drhs(k)=dpht(k,i,j)*acphT_b
          apph(k)=apphk(k)*acphT_b
          amph(k)=amphk(k)*acphT_b
          acphT(k)=1.0d0
         enddo
  
         amphT=amph(2:n3m)
         apphT=apph(1:(n3m-1))
  
         call dgttrf(n3m, amphT, acphT, apphT, appph, phpiv, info)
  
         call dgttrs('N',n3m,1,amphT,acphT,apphT,appph,phpiv,drhs, &
                       n3m, info)
  
          do k=1,n3m
            dpht(k,i,j) = drhs(k)
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo

      call PackR_UnpackZP(dpht(:,:,jstartp:jendp),dph(:,:,kstart:kend))


!
!================================================================
!   inverse fft applied to the phi x1 direction
!   from wave number space to physical space
!
      do k=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,i)
       do j=1,n2mh
        do i=1,n1m
          xa(j,i)=cmplx(dph(i,j,k),dph(i,j+n2mh,k))
        enddo
       end do
!$OMP  END PARALLEL DO

      call dfftw_execute_dft_c2r(bck_plan,xa,xr)

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,i)
       do j=1,n2m
         do i=1,n1m
           dph(i,j,k)=xr(j,i)
         enddo
       end do
!$OMP  END PARALLEL DO
      end do

      deallocate(dpht)

#ifdef _OPENMP
      call dfftw_cleanup_threads()
#endif

      return
      end
