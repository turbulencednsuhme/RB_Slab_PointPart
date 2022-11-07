!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ReadPpartInput.F90                             !
!    CONTAINS: subroutine ReadPpartInput                  !
!                                                         ! 
!    PURPOSE: Read parameters from ppart.in file          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ReadPPartInput
      use param
      use mpih
      use pointparticle
      implicit none
      character(len=4) :: dummy

      open(unit=15,file='ppart.in',status='old')
        read(15,301) dummy
        read(15,*) Npointpart,timeONp,iresetp,toutpp
        read(15,301) dummy
        read(15,*) usfroude,stokp,rhohat,thstokp
        read(15,301) dummy
        read(15,*) twowayflag, dbd0, cppcpf
        read(15,301) dummy
        read(15,*) cpi(1),cpi(2),cpi(3)
        read(15,301) dummy
        read(15,*) cpf(1),cpf(2),cpf(3)
301     format(a4)                
      close(15)

      if(Npointpart.gt.0) withppart = .true.
      if(twowayflag.gt.0) istwoway = .true.

      return
      end
