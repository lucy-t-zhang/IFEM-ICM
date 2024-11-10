      subroutine bload(x,si,nflag,wt,wt1,wt2,wg,ni,ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension hsf(3),psf(3),x(2,9)
      r=si
      s=si
      if (nflag .eq. 1) then
         hsf(3)=0.5d0*(1.0d0+s)-0.5d0*(1.0d0-s**2)
         hsf(2)=1.0d0-s**2
         hsf(1)=0.5d0*(1.0d0-s)-0.5d0*(1.0d0-s**2)
         psf(3)=0.5d0+s
         psf(2)=-2.0d0*s
         psf(1)=-0.5d0+s
         nbn(ne,1)=4
         nbn(ne,2)=8
         nbn(ne,3)=1
         xxbp=boupo(ne,1)*(1.0d0-s)/2.0d0+
     $        boupo(ne,2)*(1.0d0+s)/2.0d0
      endif
      if (nflag .eq. 2) then
         hsf(1)=0.5d0*(1.0d0+r)-0.5d0*(1.0d0-r**2)
         hsf(2)=1.0d0-r**2
         hsf(3)=0.5d0*(1.0d0-r)-0.5d0*(1.0d0-r**2)
         psf(1)=0.5d0+r
         psf(2)=-2.0d0*r
         psf(3)=-0.5d0+r
         nbn(ne,1)=1
         nbn(ne,2)=5
         nbn(ne,3)=2
         xxbp=boupo(ne,1)*(1.0d0+r)/2.0d0+
     $        boupo(ne,2)*(1.0d0-r)/2.0d0
      endif
      if (nflag .eq. 3) then
         hsf(1)=-0.5d0*(1.0d0+s)+0.5d0*(1.0d0-s**2)
         hsf(2)=-1.0d0+s**2
         hsf(3)=-0.5d0*(1.0d0-s)+0.5d0*(1.0d0-s**2)
         psf(1)=-0.5d0-s
         psf(2)=2.0d0*s
         psf(3)=0.5d0-s
         nbn(ne,1)=2
         nbn(ne,2)=6
         nbn(ne,3)=3
         xxbp=boupo(ne,1)*(1.0d0+s)/2.0d0+
     $        boupo(ne,2)*(1.0d0-s)/2.0d0
      endif
      if (nflag .eq. 4) then
         hsf(3)=-0.5d0*(1.0d0+r)+0.5d0*(1.0d0-r**2)
         hsf(2)=-1.0d0+r**2
         hsf(1)=-0.5d0*(1.0d0-r)+0.5d0*(1.0d0-r**2)
         psf(3)=-0.5d0-r
         psf(2)=2.0d0*r
         psf(1)=0.5d0-r
         nbn(ne,1)=3
         nbn(ne,2)=7
         nbn(ne,3)=4
         xxbp=boupo(ne,1)*(1.0d0-r)/2.0d0+
     $        boupo(ne,2)*(1.0d0+r)/2.0d0
      endif
      dbmx=0.0d0
      dbmy=0.0d0
      do i=1,3
         dbmx=dbmx+psf(i)*x(1,nbn(ne,i))
         dbmy=dbmy+psf(i)*x(2,nbn(ne,i))
      enddo
      detb=dsqrt(dbmx**2+dbmy**2)
      wt=wg*hsf(ni)*xxbp*detb
      wt1=-wt*dbmy/detb
	wt2=wt*dbmx/detb
	return
      end