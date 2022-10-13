!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: TimeMarcher.F90                                !
!    CONTAINS: subroutine TimeMarcher                     !
!                                                         ! 
!    PURPOSE: Main time integrating routine, which calls  !
!     other subroutines for calculating the Navier-Stokes !
!     equations and advancing velocity and temperature in !
!     time                                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TimeMarcher
      use param
      use local_arrays
      use mpih
      use local_aux
      use mpi_param, only: kstart,kend
      use pointparticle
      implicit none
      integer :: ns,i,j,k
      real :: cksum2, mck2

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Point particle routines

       if(ONparticle.eq.0)then
         call InitPointP
         ONparticle=1
       end if
       call UpdatePointParticlePosition
       call storeold
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      
      do ns=1,nsst                                                 
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)
        if(ntime.le.1) then
         aldto = alm(1)*dt
        else
         aldto = al*dt
        end if

        call DefineIBMObject

        call ExplicitTermsVX
        call ExplicitTermsVY
        call ExplicitTermsVZ
        call ExplicitTermsTemp

        call ImplicitAndUpdateVX
        call ImplicitAndUpdateVY
        call ImplicitAndUpdateVZ
        call ImplicitAndUpdateTemp

        call update_lower_ghost(n1,n2,vx)
        call update_upper_ghost(n1,n2,vz)
        call update_both_ghosts(n1,n2,temp,kstart,kend)

        call CalcLocalDivergence
        call SolvePressureCorrection

        call update_both_ghosts(n1,n2+1,dph,kstart,kend)

        call CorrectPressure
        call CorrectVelocity

        call update_both_ghosts(n1,n2,vx,kstart,kend)
        call update_both_ghosts(n1,n2,vy,kstart,kend)
        call update_both_ghosts(n1,n2,vz,kstart,kend)
        call update_lower_ghost(n1,n2,pr)

       enddo

      if(withppart) then

!CS   Calculate material derivative at staggered location
      call CalcMaterialDerivative
      call update_both_ghosts(n1,n2,matderx,kstart,kend)
      call update_both_ghosts(n1,n2,matdery,kstart,kend)
      call update_both_ghosts(n1,n2,matderz,kstart,kend)

!CS   Center velocity and material derivative
      call CentreVariables

      call CalcVorticity
      call update_upper_ghost(n1m,n2m,vorx)
      call update_upper_ghost(n1m,n2m,vory)
      call update_upper_ghost(n1m,n2m,vorz)

!CS   Now solve for velocity of particles
      if(ONparticle.eq.1) call CalcPointPVel

!CS   Now reset forcing
       for_x_part(:,:,:)=0.d0
       for_y_part(:,:,:)=0.d0
       for_z_part(:,:,:)=0.d0
      end if

!RO    Write to screen

        if(mod(time,tpin).lt.dt) then
        if(myid.eq.0) then
        write(6,*) ' ---------------------------------------- '
        write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
        endif
        endif

      return                                                            
      end                                                               
