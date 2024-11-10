      subroutine acc(ne)
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
      if (nale .eq. 1) then
		do i=1,2
            xvelm(i)=0.0d0
			do k=1,nis(ig)
               ntem=nea(ne,k)
               xvelm(i)=xvelm(i)+h(k)*velm(i,ntem,1)
			enddo
		enddo
      endif
      if (nale .eq. 1) then
         defv(1)=xvel(1)-xvelm(1)
         defv(2)=xvel(2)-xvelm(2)
      else
         defv(1)=xvel(1)
         defv(2)=xvel(2)
      endif
      do i=1,2
         do j=1,2
            xdve(i,j)=0.0d0
            do k=1,nis(ig)
               ntem=nea(ne,k)
               xdve(i,j)=xdve(i,j)+bd(j,k)*vel(i,ntem,1)
            enddo
		enddo
	enddo
	xac(1)=xac(1)-fbacco(1)
      xac(2)=xac(2)-fbacco(2)
      if ((nsmoo .ne. 5) .and. (ncon .ne. 1)) then
		xac(1)=xac(1)+defv(1)*xdve(1,1)+&
     		defv(2)*xdve(1,2)
          xac(2)=xac(2)+defv(1)*xdve(2,1)+&
            defv(2)*xdve(2,2)
      endif
      if (nfrees .eq. 1) then
         xac(2)=xac(2)+gravit
      endif
      return
      end