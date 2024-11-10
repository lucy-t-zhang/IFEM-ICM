      subroutine readinit
      implicit real*8 (a-h,o-z)
      include 'common'
 1000 read(61,*) ndum,x1,y1
      ndumtest=ndum-1
      if (ndumtest .ge. 0) then
         if (initdir .eq. 1) then
            dis(1,ndum,1)=x1
            dis(2,ndum,1)=y1
            xinvel(1,ndum)=x1
            xinvel(2,ndum)=y1
         endif
         if (initdir .eq. 2) then
            vel(1,ndum,1)=x1
            vel(2,ndum,1)=y1
            xinvel(1,ndum)=x1
            xinvel(2,ndum)=y1
         endif
         goto 1000
      endif
      return
      end	