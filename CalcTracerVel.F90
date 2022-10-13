!----------------------------------------------------------------------
!     explicit  solver for bubble velocity (drag gravity,lift,addded mass)
!     solved by each node from the position info received from master
!-----------------------------------------------------------------------

      subroutine CalcPointPVel
      USE param
      USE mpih
      USE mpi_param
      use pointparticle

      IMPLICIT NONE

      integer inp
      integer bubind(6)
      real pos(3),qInt(7),dmatInt(3)
      real uxn,uyn,uzn,vxn,vyn,vzn,omexn,omeyn,omezn
      real uxn1,uyn1,uzn1,vxn1,vyn1,vzn1,omexn1,omeyn1,omezn1
      real kalbxn,kalbyn,kalbzn,kalbxn1,kalbyn1,kalbzn1
      real tau_b,fdrag,renp_temp
      real Forxn,Foryn,Forzn,Forxn1,Foryn1,Forzn1
      real Adn1                 ! Diagonal matrix entries
      real Aoxn1,Aoyn1,Aozn1    ! Off-diagonal matrix entries
      real RHSx,RHSy,RHSz

      integer i,j,k,ist,jst,kst
      integer dum,inb,li,lj,lk,istwoway
      real fBox(8,3),matderpart(3),posp(3)
      real Volbub,dbub,Volcell
      real dummy


      do inp=1,Npointpart

      ! Get positions
      pos(1)=xp(inp,1)
      pos(2)=xp(inp,2)
      pos(3)=xp(inp,3)

      ! Get indices from indpos
      call indpos(pos,inp,bubind)

      ! Only specific processors will calculate particles
      if(bubind(6).ge.kstart-1.and.bubind(6).le.kend-1)then

      ! Interpolate velocity and vorticity onto bubble position
      call interpol(pos,bubind,inp,qInt)

      ! Store values for simplicity
      vxn1=qInt(1)    !u_{f,1,interp}
      vyn1=qInt(2)    !u_{f,2,interp}
      vzn1=qInt(3)    !u_{f,3,interp}

      q1p(inp)=vxn1
      q2p(inp)=vyn1
      q3p(inp)=vzn1

!     ----------------------------------------------------------------------

      else

      q1p(inp)=0.d0 ; q2p(inp)=0.d0 ; q3p(inp)=0.d0 
      
      end if

      end do


      call MpiBarrier

      call MpiAllSumReal1D(q1p,Npointpart)
      call MpiAllSumReal1D(q2p,Npointpart)
      call MpiAllSumReal1D(q3p,Npointpart)

      return
      end
