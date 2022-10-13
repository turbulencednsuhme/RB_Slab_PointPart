!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: OpenLogs.F90                                   !
!    CONTAINS: subroutine OpenLogs, CloseLogs             !
!                                                         ! 
!    PURPOSE: (1) Initialization routine. Open all log    !
!     files and reset if necessary. (2) Finalization      !
!     routine. Close all log files                        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine OpenLogs
      use param
      implicit none

!EP    nusse.out  in vmaxv.F
      open(95,file='nu_vol.out',status='unknown',access='sequential', &
       position='append')

!EP    nusse2.out  in densmc.F
       open(97,file="nu_plate.out",status='unknown',access='sequential', &
       position='append')

!EP   rms_vel.out in stst.F
       open(94,file='rms_vel.out',status='unknown',position='append', &
       access='sequential')

!EP   nusse3.out in balance.F
      open(92,file='nu_diss.out',status='unknown',access='sequential', &
       position='append')

      return
      end   
      

      subroutine CloseLogs
      implicit none
      close(95)
      close(97)
      close(94)
      close(92)
      return 
      end
