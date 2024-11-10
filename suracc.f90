      subroutine suracc(nes,det)
      implicit real*8 (a-h,o-z)
      include 'common'
      xv1=0.0d0
      xv2=0.0d0
      xdsur=0.0d0
      xds=0.0d0
      xvms=0.0d0
      ns1=nfsnd(nes,1)
      ns2=nfsnd(nes,2)
      hrs(1)=dabs(coor(ns1,1)-coor(ns2,1))
      do i=1,niss
         nstem=nfsnd(nes,i)
         xv1=xv1+hfs(i)*vel(1,nstem,1)
         xv2=xv2+hfs(i)*vel(2,nstem,1)
         xdsur=xdsur+hfsp(i)*dis(2,nstem,1)/det
         xds=xds+hfs(i)*dis(2,nstem,1)
         xvms=xvms+hfs(i)*velm(2,nstem,1)
      enddo
      if (nsmoos .eq. 1) then
         if (xv1 .gt. 0.0d0) then
            stem=1.0d0
         else
            stem=-1.0d0
         endif
         if (dabs(xv1) .gt. 0.0d0) then
            upws=stem/dsqrt(15.0d0)*hrs(1)
         else
            upws=0.0d0
         endif
      endif
      if (nsmoos .eq. 4) then
         xuk(1,1)=xco*dabs(xv1)*hrs(1)
      endif
      return
      end