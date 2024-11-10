      subroutine jacob(x,xj,xji,det)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension x(2,9),xj(2,2),xji(2,2)
      do i=1,2
         do j=1,2
            dum=0.0d0
            do k=1,nis(ig)
               dum=dum+p(i,k)*x(j,k)
            enddo
            xj(i,j)=dum
		enddo
	enddo
      det=xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
      if (det .lt. 1.0d-15) then
         write(*,100) 
         stop
      endif
  100 format(6x, 'error, zero or negative jacobian determinant')
      xji(1,1)=xj(2,2)/det
      xji(1,2)=-xj(1,2)/det
      xji(2,1)=-xj(2,1)/det
      xji(2,2)=xj(1,1)/det
      return
      end