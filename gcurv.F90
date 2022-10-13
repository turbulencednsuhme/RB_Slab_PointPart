      subroutine gcurv
      use mpih
      use mpi_param,only: kstart,kend
      use param
      use local_arrays, only: vy,vz,temp,pr,vx
      use hdf5
      use pointparticle
      use stat_arrays, only: timeint_cdsp

      implicit none
      integer :: ntstf, l, hdf_error
      real    :: cflm,dmax
      real    :: ti(2), tin(3)

      tin(1) = MPI_WTIME()
      
      call CreateGrid

      call h5open_f(hdf_error)
      call initstst

      if(myid.eq.0) ismaster = .true.
      if(myid.eq.0) then
      write(6,754)n1,n2,n3                                              
  754 format(/,5x,'grid resolution: ',' n1= ',i5,' n2= ',i5, &
       ' n3= ',i5/)                       
      write(6,755) 1.d0/dx1,1.d0/dx2,1.d0/dx3,dt,ntst                  
  755 format(/,2x,' dx1=',e10.3,' dx2=',e10.3,' dx3=',e10.3,' dt=' &
       ,e10.3,' ntst=',i7,/)
      endif

      
      time=0.d0
      vmax=0
!                                                                       
!   read or create initial fields                                       
!                                                                       

      call InitPressureSolver
      call SetTempBCs
!
!      create the initial conditions
!
      if(nread.eq.0) then
 
       if(ismaster) write(6,*)' nread=0 ---> Create Initial Conditions'
 
       ntime=0                                                         
       time=0.d0
       cflm=0.d0
         
       call CreateInitialConditions
       call update_both_ghosts(n1,n2,temp,kstart,kend)
       call update_both_ghosts(n1,n2,pr,kstart,kend)
       call update_both_ghosts(n1,n2,vx,kstart,kend)
       call update_both_ghosts(n1,n2,vy,kstart,kend)
       call update_both_ghosts(n1,n2,vz,kstart,kend)
        
      else

       if(ismaster) write(6,*)' nread=1 ---> Read Initial Conditions'

       call ReadInitCond
       call update_both_ghosts(n1,n2,temp,kstart,kend)
       call update_both_ghosts(n1,n2,pr,kstart,kend)
       call update_both_ghosts(n1,n2,vx,kstart,kend)
       call update_both_ghosts(n1,n2,vy,kstart,kend)
       call update_both_ghosts(n1,n2,vz,kstart,kend)

       tmax = tmax + time
       
      endif                                                             

      call CheckDivergence(dmax)
      if(ismaster) write(6,*)' initial divg dmax  ',dmax
 
      ntstf=ntst                                                   

      if(ismaster) then
       write(6,711) tprint,ntstf,tpin
711   format(3x,'check in cond : tprint =',f10.1,  &
             '  ntstf =',i8,2x,'tpin =',f10.1//)
      endif
!m================================ 
      if(idtv.eq.1) then

       if(ismaster) then 
        write(6,*)ntime,time,vmax(1),vmax(2),vmax(3), &
         dt,dmax,tempm,denmax,denmin
       endif

      else

       if(ismaster) then
         write(6,*)ntime,time,vmax(1),vmax(2),vmax(3), &
          cflm,dmax, tempm,denmax,denmin
       endif
     
       cflm=cflm*dt

      endif

      if(ismaster) then
       tin(2) = MPI_WTIME()
       write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
      endif
! ==================================
                                                                      
      do ntime=0,ntstf                                           
!
!     the calculation stops if the velocities are diverging for numerical
!     stability conditions (courant number restrictions)                
!
       ti(1) = MPI_WTIME()
!
       call CalcMaxCFL(cflm)
            
       if(idtv.eq.1) then

        if(ntime.ne.1) then
         dt=cflmax/cflm
         if(dt.gt.dtmax) dt=dtmax
        endif

        if(dt.lt. 0.00000001d0) go to 166

       else

        cflm=cflm*dt
        if(cflm.gt.cfllim) go to 165

       endif

       beta=dt/ren*0.5d0

       call TimeMarcher

       time=time+dt

       if(ntime.eq.1) go to 306 

       if(mod(time,tpin).lt.dt) go to 306  

       go to 305

  306 continue                                                          

       call GlobalQuantities
       if(vmax(1).gt.1000.d0.and.vmax(2).gt.1000.d0) go to 266
       call CalcMaxCFL(cflm)
       call CheckDivergence(dmax)
       call CalcPlateNu

       if(time.gt.tsta) then
        call CalcStats
        call CalcDissipationNu
       endif

       if(idtv.eq.1) then
       else
        cflm=cflm*dt
       endif

       if(dmax.gt.resid) go to 169

  305 continue                                                          

       if(mod(time,tframe).lt.dt) then
        call Dump3DFiles
       endif
       if(mod(time,toutpp).lt.dt) then
          call part_write_continua(.true.)
       end if
          
       if(time.gt.tmax) go to 333

       ti(2) = MPI_WTIME()

       if(mod(time,tpin).lt.dt.and.ismaster) then
        write(6,*) 'Iteration Time = ', ti(2) -ti(1), ' sec.'
       endif

       if( (ti(2) -tin(1)) .gt. walltimemax) then

        call WriteStats
        call mpi_write_continua

        if(ismaster) then
           write(6,*) 'Restart file updated at exit'
           write(6,*) 'exit time=',ti(2) -tin(1),'sec'
        endif

        stop ! Exit after time limit

       endif

      end do  ! End of time marching loop

 
333   continue
  
      tin(3) = MPI_WTIME()

      if(ismaster) then
        write(6,*) 'Total Iteration Time = ',tin(3) -tin(2),' sec.'
      endif

      call WriteStats
      call mpi_write_continua
      go to 167                                                         

  165 continue 

!m======================================
      if(ismaster) then
       write(6,164)                             
  164 format(10x,'cfl too large  ')
      endif
!m======================================                                  
      go to 167                                                         
  166 continue
!m======================================
      if(myid.eq.0) then
      write(6,168) dt 
  168 format(10x,'dt too small, DT= ',e14.7)
      endif
!m======================================   
      go to 167                                                         
  266 continue
!m======================================
      if(myid.eq.0) then
      write(6,268)                                                      
  268 format(10x,'velocities diverged')
      endif
!m======================================                                 
      go to 167                                                         
  169 continue
  
!m======================================
      if(myid.eq.0) then
      write(6,178) dmax                                 
  178 format(10x,'too large local residue for mass conservation : ' &
             ,e12.5,' at ')
     
      endif
      call LocateLargeDivergence
!m======================================

  167 continue                                                          

      call dfftw_destroy_plan(fwd_plan)
      call dfftw_destroy_plan(bck_plan)

      call h5close_f(hdf_error)
      
      return                                                            
      end                                                               
