subroutine  Readcsv(filename,nrow,ncol,dat)
implicit none
character(LEN=20),intent(in)      :: filename
integer,          intent(in)      :: nrow, ncol
character(LEN=10),dimension(ncol) :: header
integer,          parameter       :: stdout=6, funit = 9
real, dimension(nrow,ncol), intent(out):: dat
integer    :: err, ix
logical    :: I_opened
real       :: cff(ncol)

   dat(:,:) = 0.0
   ! Inquire whether the unit has been opened or not
   INQUIRE (funit, OPENED=I_opened) 

   if (I_opened) then
      print *, 'The file unit already open!'
      close(funit)
   endif

   OPEN(unit=funit,file=filename,status='OLD',iostat=err)

   IF (err /= 0) THEN
     write(stdout,*) 'open ', TRIM(filename),' fails'
     stop
     close(funit)

   ELSE
     read(funit,*,iostat=err) header  !read the first row
    ! write(stdout,*)  header


     if (err .gt. 0) then 
       write(stdout,*) 'The error is: ', err, 'at Row: ',ix, ' for ',TRIM(filename)
       CLOSE(funit)
       stop
     elseif (err .lt. 0) then
       write(stdout,*) 'The error is: ', err, 'at Row: ',ix, ' for ',TRIM(filename)
       write(stdout,*) 'End of file reached!'
       CLOSE(funit)
       stop
     else
       continue
     endif

     do ix = 1,nrow  
       read(funit,*,iostat=err) cff(:)
       if (err .gt. 0) then 
         write(stdout,*) 'The error is: ', err, 'at Row: ',ix, ' for ',TRIM(filename)
         CLOSE(funit)
         stop
       elseif (err .lt. 0) then
         write(stdout,*) 'The error is: ', err, 'at Row: ',ix, ' for ',TRIM(filename)  
         write(stdout,*) 'End of file reached!'
         CLOSE(funit)
         stop
       else
         dat(ix,:) = cff(:)
       endif
     enddo  

   ENDIF
   CLOSE(funit)
! End reading data
End subroutine Readcsv
