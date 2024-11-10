      subroutine material
      implicit real*8 (a-h,o-z)
      include 'common'
      young=syoun
      pois=spois
      if (npss .eq. 1) then
         cc=young*(1.0d0-pois)/(1.0d0+pois)/(1.0d0-2.0d0*pois)
         cmat(1,1)=cc
         cmat(1,2)=cc*pois/(1.0d0-pois)
         cmat(1,3)=0.0d0
         cmat(2,1)=cmat(1,2)
         cmat(2,2)=cmat(1,1)
         cmat(2,3)=0.0d0
         cmat(3,1)=0.0d0
         cmat(3,2)=0.0d0
         cmat(3,3)=cc*(1.0d0-2.0d0*pois)/2.0d0/(1.0d0-pois)
      endif
      if (npss .eq. 0) then
         cc=young/(1.0d0-pois)/(1.0d0+pois)
         cmat(1,1)=cc
         cmat(1,2)=cc*pois
         cmat(1,3)=0.0d0
         cmat(2,1)=cmat(1,2)
         cmat(2,2)=cmat(1,1)
         cmat(2,3)=0.0d0
         cmat(3,1)=0.0d0
         cmat(3,2)=0.0d0
         cmat(3,3)=cc*(1.0d0-pois)/2.0d0
      endif
      return
      end