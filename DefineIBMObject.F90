      subroutine DefineIBMObject
      use param
      use ibm_param
      use mpi_param, only: kstart,kend

      implicit none
      integer :: i,j,k,l,n,kstartp
      integer :: km,kp,jm,jp,im,ip

      real, allocatable, dimension(:,:) :: strct
      integer :: sy, sz,xtr

      real    :: xe, xem, xep
      real    :: ye, yem, yep
      real    :: ze, zem, zep
      real    :: delta1x, delta2x
      real    :: cylsp, cylr, yp, cylh, cyloff
      real    :: vy_cyls, tphas,periodmov, movfrac
      real    :: safety_margin, yc_cyls
      real    :: tstart


!    Geometrical Constants for Flute

      npunt=0
      npunr=0
      npunz=0
!
!     IDENTIFICATION OF THE GRID POINTS IN THE BODY
!     (BOUNDARY + INNER PART)
!
!
      if(.not.allocated(forclo)) allocate(forclo(1:n1,1:n2,kstart:kend))

      do l = 1,4 !{ start do over the 3 velocity components
      n=0
      forclo = 0.0d0

!     l = 1   Q_1 vel. component
!     l = 2   Q_2 vel. component
!     l = 3   Q_3 vel. component


      if(l.eq.1) then
        if(n.gt.mpun) &
        write(*,*) 'Dim max di indgeot e'' stata superata n=',n
        npunfx(1)= n
!        write(6,332)npunfx(1)
 332  format(5x,'For Vx N ='i7)
      end if
      if(l.eq.2) then
        if(n.gt.mpun) &
        write(*,*) 'Dim max di indgeor e'' stata superata n=',n
        npunfx(2)= n
!        write(6,331)npunfx(2)
 331  format(5x,'For Vy N ='i7)
      end if
      if(l.eq.3) then
        if(n.gt.mpun) &
        write(*,*) 'Dim max di indgeoz e'' stata superata n=',n
        npunfx(3)= n
!        write(6,330)npunfx(3)
 330  format(5x,'For Vz N ='i7)
      end if
      if(l.eq.4) then
        if(n.gt.mpun) &
        write(*,*) 'Dim max di indgeoz e'' stata superata n=',n
        npunfx(4)= n
!        write(6,333)npunfx(4)
 333  format(5x,'For Dens N ='i7)
      end if
      end do   !} end do over the 3 velocity components and scalar

!      write(*,*) 'IBM Object at zcyl=',zc_cyls

      return
      end

