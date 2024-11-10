      subroutine sstrain(toc,xto,lx,ly,ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension xto(2,2),toc(3,3)
      do i=1,2
         ge(i,ne,lx,ly)=0.5d0*(toc(i,i)-1.0d0)
      enddo
      ge(3,ne,lx,ly)=toc(1,2)
      do i=1,2
         do k=1,nis(ig)
            do j=1,2
			dge(i,j,k)=xto(j,i)*bd(i,k)
			ddge(i,j,k)=bd(i,k)*bd(i,k)
		  enddo
	   enddo
      enddo
      do k=1,nis(ig)
		do j=1,2
              dge(3,j,k)=xto(j,1)*bd(2,k)+xto(j,2)*bd(1,k)
			ddge(3,j,k)=bd(1,k)*bd(2,k)+bd(2,k)*bd(1,k)
	    enddo
	enddo
      return
      end