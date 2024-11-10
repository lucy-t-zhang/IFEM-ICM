      subroutine scalkuu(xkuu,ocuu,i,ni)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension ocuu(6,6)
      xkuu=0.0d0
      do k=1,3
        xkuu=xkuu+tos(k)*ddge(k,i,ni)
         do m=1,3
            xkuu=xkuu+ocuu(k,m)*dge(k,i,ni)*dge(m,i,ni)
		enddo
	enddo
      return
      end