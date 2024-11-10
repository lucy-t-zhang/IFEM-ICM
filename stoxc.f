      subroutine stoxc(xto,xot,xj,xji,toxj,toxji,toc)
      implicit real*8 (a-h,o-z) 
      dimension xto(2,2),xot(2,2),xj(2,2),xji(2,2),toc(3,3),
     $	toxj(2,2),toxji(2,2)
      do i=1,2
         do j=1,2
            xto(i,j)=0.0d0
            xot(i,j)=0.0d0
            do m=1,2
               xto(i,j)=xto(i,j)+xj(m,i)*toxji(j,m)
               xot(i,j)=xot(i,j)+toxj(m,i)*xji(j,m)
		  enddo
		enddo
	enddo
      do i=1,2
         do j=1,2
            toc(i,j)=0.0d0
            do m=1,2
               toc(i,j)=toc(i,j)+xto(m,i)*xto(m,j)
		  enddo
		enddo
	enddo
      return
      end