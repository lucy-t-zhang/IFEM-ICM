      subroutine bdpd(xji)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension xji(2,2)
      do i=1,nis(ig)
         do k=1,2
            dumcd=0.0d0
            do j=1,2
               dumcd=dumcd+xji(k,j)*p(j,i)
            enddo
            bd(k,i)=dumcd
         enddo
	enddo
      return
      end