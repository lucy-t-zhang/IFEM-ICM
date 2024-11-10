      subroutine belement(r,det,nes)
      implicit real*8 (a-h,o-z) 
      include 'common'
      if (niss .eq. 3) then
         hfs(1)=0.5d0*(r+r**2)
         hfs(2)=0.5d0*(-r+r**2)
         hfs(3)=1.0d0-r**2
         hfsp(1)=0.5d0*(1.0d0+2.0d0*r)
         hfsp(2)=0.5d0*(-1.0d0+2.0d0*r)
         hfsp(3)=-2.0d0*r
      endif     
      if (niss .eq. 2) then
         hfs(1)=0.5d0*(1.0d0+r)
         hfs(2)=0.5d0*(1.0d0-r)
         hfsp(1)=0.5d0
         hfsp(2)=-0.5d0
      endif
      det=0.0d0
      do i=1,niss
         nstem=nfsnd(nes,i)
         det=det+hfsp(i)*coor(nstem,1)
      enddo
      return
      end