!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitPressureSolver.F90                         !
!    CONTAINS: subroutine InitPressureSolver              !
!                                                         ! 
!    PURPOSE: Initialization routines. Compute the metric !
!     terms and modified wavenumbers for the pressure     !
!     correction                                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitPressureSolver
      use param
      implicit none
      integer FFTW_EXHAUSTIVE
      parameter(FFTW_EXHAUSTIVE=64)
      integer  :: kc,km,kp,i,j
      real :: ugmmm,a33icc,a33icp
      real, dimension(n2m,n1m) :: xr
      complex, dimension(n2mh,n1m) :: xa
#ifdef _OPENMP
      integer :: nt,fftw_info
#endif

      do kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3rc(kc)
        a33icp=kpc(kc)*dx3q/g3rc(kp)
        ugmmm=1.0d0/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo
    
      do i=1,n1mh
        ao(i)=(i-1)*2.d0*pi
      enddo
      do i=n1mp,n1m
        ao(i)=-(n1m-i+1)*2.d0*pi
      enddo
      do i=1,n1m
        ak1(i)=2.d0*(1.d0-dcos(ao(i)/n1m))*(float(n1m)/rext)**2
      enddo

      do j=1,n2mh
        ap(j)=(j-1)*2.d0*pi
      enddo
      do j=n2mp,n2m
        ap(j)=-(n2m-j+1)*2.d0*pi
      enddo
      do j=1,n2m
        ak2(j)=2.d0*(1.d0-dcos(ap(j)/n2m))*(float(n2m)/rext2)**2
      enddo

#ifdef _OPENMP
      nt = 0
!$OMP  PARALLEL
!$OMP$ REDUCTION(+:nt)
      nt = nt + 1
!$OMP  END PARALLEL
      call dfftw_init_threads(fftw_info)
      if(fftw_info.eq.0) write(*,*) "ERROR: FFTW THREAD INIT FAIL"
      call dfftw_plan_with_nthreads(nt)
#endif

      call dfftw_plan_dft_r2c_2d(fwd_plan,n2m,n1m,xr,xa,FFTW_EXHAUSTIVE)
      call dfftw_plan_dft_c2r_2d(bck_plan,n2m,n1m,xa,xr,FFTW_EXHAUSTIVE)

      return
      end
