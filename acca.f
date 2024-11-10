      subroutine acca(ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      do i=1,2
		xac(i)=0.0d0
          do k=1,nis(ig)
			ntem=nea(ne,k)
              xac(i)=xac(i)+h(k)*acm(i,ntem,1)
          enddo
	enddo
      do i=1,2
         xvel(i)=0.0d0
         do k=1,nis(ig)
            ntem=nea(ne,k)
            xvel(i)=xvel(i)+h(k)*vel(i,ntem,1)
         enddo
	enddo
      defv(1)=xvel(1)
      defv(2)=xvel(2)
      do i=1,2
         do j=1,2
            xdve(i,j)=0.0d0
		  xdde(i,j)=0.0d0
            do k=1,nis(ig)
               ntem=nea(ne,k)
               xdve(i,j)=xdve(i,j)+bd(j,k)*vel(i,ntem,1)
               xdde(i,j)=xdde(i,j)+bd(j,k)*dis(i,ntem,1)
            enddo
		enddo
	enddo
	xac(1)=xac(1)-fbacco(1)
      xac(2)=xac(2)-fbacco(2)
      return
      end