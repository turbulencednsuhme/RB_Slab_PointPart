!===========================================================
! Declaration of global variables
!***********************************************************      
      module param
        implicit none

        integer   :: n2, n3,n1
        integer   :: nsst, nwrit, nread, ntst
        real      :: tprint,tpin,tmax
        real      :: alx3,str3
        integer   :: istr3
        real      :: rext,rext2
        real      :: ray,pra,dt,resid,cflmax
        integer   :: inslws,inslwn
        integer   :: starea,tsta
        real      :: dtmax,cfllim ,walltimemax
        integer   :: nson,idtv
        real   :: tframe
        logical :: ismaster = .false.

        real :: time
        real :: dx2,dx3,dx1
        real :: dx2q,dx3q,dx1q

        real, allocatable, dimension(:) :: tc,tm
        real, allocatable, dimension(:) :: rc,rm
        real, allocatable, dimension(:) :: zc,zm
        real, allocatable, dimension(:) :: g3rc,g3rm
        real, allocatable, dimension(:) :: udx3c, udx3m
        real, allocatable, dimension(:) :: etaz, etazm

        integer, allocatable, dimension(:) :: imv,ipv
        integer, allocatable, dimension(:) :: jmv,jpv
        integer, allocatable, dimension(:) :: kmv,kpv
        integer, allocatable, dimension(:) :: jmhv
        integer, allocatable, dimension(:) :: kup,kum
        integer, allocatable, dimension(:) :: kmc,kpc

        real, allocatable, dimension(:) :: ap3j,ac3j,am3j
        real, allocatable, dimension(:) :: ap3ck,ac3ck,am3ck
        real, allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
        real, allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk   
        real, dimension(13) :: ifx1
        real, allocatable, dimension(:) :: ak1, ao
        real, allocatable, dimension(:) :: ak2, ap
        real, allocatable, dimension(:) :: amphk, acphk, apphk
        real, allocatable, dimension(:) :: trigx1
        integer*8 :: fwd_plan,bck_plan
        
!===========================================================
!******* Other variables ***********************************
        integer  :: n2m, n3m, n1m
        integer  :: n1mh, n1mp, n2mh, n2mp
        integer  :: iaxsy
        real :: rint
        real :: ren, pec
        real :: pi
        real :: al,ga,ro
        real :: beta,aldto
        real :: qqmax,qqtot
        real :: re
        real :: anusslow,anussupp
        real :: denmax,denmin,tempm
        integer :: ntime
        integer, parameter:: ndv=3
        real, dimension(1:3) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        real, allocatable, dimension(:,:) :: denbs,denbn
        real, allocatable, dimension(:,:) :: uwtb, lwtb
      end module param
      
!************* End of param module******************************
!===============================================================
!******* 2D arrays, dynamically allocated by each process*******
      module local_arrays
      use param
        implicit none
        real,allocatable,dimension(:,:,:) :: vx,vy,vz
        real,allocatable,dimension(:,:,:) :: pr,temp,hro,rhs
        real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3,ruro
        real,allocatable,dimension(:,:,:) :: qcap
        real,allocatable,dimension(:,:,:) :: dph,dq
      end module local_arrays

!===============================================================
      module stat_arrays
       implicit none
       real,allocatable, dimension(:,:) :: vx_me,vx_rms 
       real,allocatable, dimension(:,:) :: vy_me,vz_me,vy_rms,vz_rms 
       real,allocatable, dimension(:,:) :: temp_me,temp_rms 
       real, allocatable,dimension(:,:) :: disste,dissth,tempvz_me
       integer :: timeint_cdsp
      end module stat_arrays
!=====================================================       
      module mpih
        implicit none
        include 'mpif.h'
        integer :: myid, numtasks, ierr
        integer, parameter :: master=0
        integer :: MDP = MPI_DOUBLE_PRECISION
        integer :: STATUS(MPI_STATUS_SIZE,4)
        integer :: req(1:4)
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
      end module mpih
      
      module mpi_param
        implicit none
        integer :: istart,iend,jstart,jend, kstart,kend
        integer :: jstartp,jendp
        integer :: dj,dk,mydata,mydatam
        integer :: djp
        integer, allocatable, dimension(:) :: offsetj,offsetk
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: countj,countk
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer(8), allocatable, dimension(:) :: offsetf 
      end module mpi_param
!====================================================
!====================================================
      module ibm_param
        use param
        implicit none
        integer :: npunt,npunr,npunz,mpun,npund
        integer :: nbfx,mpugeo
        integer :: niix, nifx, njix, njfx, nkix, nkfx
        parameter (mpun=2084352)
        parameter (mpugeo=45000)
        real,allocatable,dimension(:,:,:) :: forclo
        real,dimension(mpugeo,9) :: xyzbfx
        real,dimension(mpugeo,3) :: barfx, nxyzfx
        real,dimension(mpugeo) :: areafx
        real, dimension(4,mpun) :: vbfx, vbifx
        integer,dimension(4,mpun,3) :: indgeo, indgeoe, indgeoee
        integer,dimension(4) :: npunifx,npunfx
        real,dimension(4,mpun) :: distb
        integer, dimension(4,mpun):: ntribfx
        real,dimension(mpun) :: vxbo,vybo,vzbo,dnbo
        real, dimension(mpun) :: vxpr, vypr, vzpr, dnpr

      end module ibm_param

      module local_aux
       use param
       implicit none
       real, allocatable, dimension(:,:,:) :: vxc, vyc, vzc
       real, allocatable, dimension(:,:,:) :: vxo, vyo, vzo
       real, allocatable, dimension(:,:,:) :: matderxc, matderx
       real, allocatable, dimension(:,:,:) :: matderyc, matdery
       real, allocatable, dimension(:,:,:) :: matderzc, matderz
       real,allocatable,dimension(:,:,:) :: vorx, vory, vorz

      end module local_aux

!====================================================
      module pointparticle
       use param
       implicit none
       logical :: withppart = .false.
       logical :: istwoway = .false.
       integer :: Npointpart,Onparticle
       integer :: iresetp
       integer :: twowayflag
       real  :: timeONp,toutpp

       real, allocatable, dimension(:,:,:) :: for_xc_part,for_yc_part,for_zc_part
       real, allocatable, dimension(:,:,:) :: for_tc_part
       real, dimension(:), allocatable :: qVal1,qVal2,qVal3      ! for timestepping
       real, dimension(:), allocatable :: qValo1,qValo2,qValo3   ! for timestepping
       real, dimension(:), allocatable :: kalb1,kalb2,kalb3
       real, dimension(:), allocatable :: kalbo1,kalbo2,kalbo3
       real, dimension(:), allocatable :: vort1,vort2,vort3
       real, dimension(:), allocatable :: vorto1,vorto2,vorto3!
       real, dimension(:,:), allocatable :: xp,xpo
       real, dimension(:), allocatable :: vxp,vyp,vzp,dtempp    ! part vel new
       real, dimension(:), allocatable :: vxpo,vypo,vzpo,dtemppo ! part vel old
       real, dimension(:), allocatable :: renp       ! part Re number
       real, dimension(:,:), allocatable :: aap      ! part acc
       real, dimension(:,:), allocatable :: facc_for,drag_for,lift_for
       real, dimension(:,:), allocatable :: buoy_for

       real, dimension(:), allocatable :: dbd,stokes
       real, dimension(:), allocatable :: gammap,temptime
       real cpi(3), cpf(3)

       real :: usfroude,dbd0,rhohat,stokp,thstokp,cppcpf
       

 
      end module pointparticle
