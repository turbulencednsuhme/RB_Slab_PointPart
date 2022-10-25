!----------------------------------------------------
!     Routine for reading in particle position 
!     and velocity data.
!----------------------------------------------------

      subroutine part_read_continua

      USE hdf5
      USE param
      USE pointparticle

      IMPLICIT NONE
      character*30 :: filename,dsetnam

      filename = trim('outputdir/parthist.h5')
 
      if(ismaster) then

       dsetnam=trim('xp')
       call HdfSerialReadReal2D(dsetnam,filename,xp,Npointpart,4)
 
       dsetnam=trim('vxp')
       call HdfSerialReadReal1D(dsetnam,filename,vxp,Npointpart)
       dsetnam=trim('vyp')
       call HdfSerialReadReal1D(dsetnam,filename,vyp,Npointpart)
       dsetnam=trim('vzp')
       call HdfSerialReadReal1D(dsetnam,filename,vzp,Npointpart)

      end if

      call MpiBcastReal2DArray(xp,Npointpart,4)
      call MpiBcastReal1DArray(vxp,Npointpart)
      call MpiBcastReal1DArray(vyp,Npointpart)
      call MpiBcastReal1DArray(vzp,Npointpart)

      call MpiBarrier

      return
      end subroutine part_read_continua

!----------------------------------------------------
!     Routine for writing particle position 
!     and velocity data.
!----------------------------------------------------

      subroutine part_write_continua(ismidrun)

      USE hdf5
      USE param
      USE pointparticle

      IMPLICIT NONE
      character*50 :: filename,dsetnam
      character*5  :: ipfi
      integer hdf_error, itime
      real tprfi
      logical :: ismidrun

      if(ismaster) then 

      if(ismidrun) then 
        tprfi = 100.
        itime=nint(time*tprfi)
        write(ipfi,82)itime
   82   format(i5.5)
       filename = trim('outputdir/parthist_'//ipfi//'.h5')
      else 
       filename = trim('outputdir/parthist.h5')
      end if

      call HdfCreateBlankFile(filename)

      dsetnam=trim('xp')
      call HdfSerialWriteReal2D(dsetnam,filename,xp,Npointpart,4)
 
      dsetnam=trim('vxp')
      call HdfSerialWriteReal1D(dsetnam,filename,vxp,Npointpart)
      dsetnam=trim('vyp')
      call HdfSerialWriteReal1D(dsetnam,filename,vyp,Npointpart)
      dsetnam=trim('vzp')
      call HdfSerialWriteReal1D(dsetnam,filename,vzp,Npointpart)

      end if
      call MpiBarrier

      return
      end subroutine part_write_continua



!================================================================================

      subroutine HdfCreateBlankFile(filename)
      use hdf5
      implicit none
      character*30,intent(in) :: filename
      integer(HID_T) :: file_id
      integer :: hdf_error

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfCreateBlankFile

!====================================================================
      subroutine HdfSerialWriteReal2D(dsetname,filename,var,sx,sy)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sx,sy
      real, dimension(sx,sy), intent(in) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset, filespace
      integer :: hdf_error
      integer(HSIZE_T) :: dims(2)
      logical :: fileexists


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sx
      dims(2)=sy

      call h5screate_simple_f(2, dims, filespace, hdf_error)

      call h5lexists_f(file_id,dsetname,fileexists,hdf_error)

      if(fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE,  &
                      filespace, dset, hdf_error)


       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE,  &
         var(1:sx,1:sy), dims, hdf_error)


      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialWriteReal2D

!====================================================================
      subroutine HdfSerialReadReal2D(dsetname,filename,var,sx,sy)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sx,sy
      real, dimension(sx,sy), intent(out) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(2)


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sx
      dims(2)=sy

      call h5dopen_f(file_id, dsetname, dset, hdf_error)

      call h5dread_f(dset, H5T_NATIVE_DOUBLE, var, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialReadReal2D
!====================================================================

      subroutine HdfStart
      use hdf5
      implicit none
      integer :: hdf_error

      call h5open_f(hdf_error)

      end subroutine HdfStart

!====================================================================
      subroutine HdfClose
      use hdf5
      implicit none
      integer :: hdf_error

      call h5close_f(hdf_error)

      end subroutine HdfClose

!====================================================================
      subroutine HdfSerialWriteReal1D(dsetname,filename,var,sz)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sz
      real, dimension(sz), intent(in) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset, filespace
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)
      logical :: fileexists


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sz

      call h5screate_simple_f(1, dims, &
     &                        filespace, hdf_error)

      call h5lexists_f(file_id,dsetname,fileexists,hdf_error)

      if(fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, &
     &                filespace, dset, hdf_error)


       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
     &   var(1:sz), dims, hdf_error)


      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialWriteReal1D
!====================================================================

      subroutine HdfSerialReadReal1D(dsetname,filename,var,sz)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sz
      real, dimension(sz), intent(out) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sz

      call h5dopen_f(file_id, dsetname, dset, hdf_error)

      call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
     &   var, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialReadReal1D

!====================================================================
