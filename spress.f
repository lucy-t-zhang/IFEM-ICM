      subroutine spress(rs,ne)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension rs(2)
      r=rs(1)
      s=rs(2)
      hp(1)=0.25d0*(1.0d0+r)*(1.0d0+s)
      hp(2)=0.25d0*(1.0d0-r)*(1.0d0+s)
      hp(3)=0.25d0*(1.0d0-r)*(1.0d0-s)
      hp(4)=0.25d0*(1.0d0+r)*(1.0d0-s)     
      cpre=0.0d0
      do i=1,nump
		ntt=neap(ne,i)
          cpre=cpre+epc(ntt,1)*hp(i)
      enddo
      return
      end