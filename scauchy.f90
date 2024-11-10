      subroutine scauchy(det,todet,xto,lx,ly,ne)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension xto(2,2)
      ss(1,1)=tos(1)
      ss(2,2)=tos(2)
      ss(1,2)=tos(3)
      ss(2,1)=ss(1,2)
      ss(3,3)=tos(4)
      do i=1,4
         cstr(i,ne,lx,ly)=0.0d0
      enddo
      do i=1,2
		do j=1,2
            tt(i,j)=xto(i,j)
		enddo
	enddo
      do m=1,3
		do n=1,3
            cstr(1,ne,lx,ly)=cstr(1,ne,lx,ly)+&
                todet/det*tt(1,m)*ss(m,n)*tt(1,n)
            cstr(2,ne,lx,ly)=cstr(2,ne,lx,ly)+&
                todet/det*tt(2,m)*ss(m,n)*tt(2,n)
            cstr(3,ne,lx,ly)=cstr(3,ne,lx,ly)+&
                todet/det*tt(1,m)*ss(m,n)*tt(2,n)
            cstr(4,ne,lx,ly)=cstr(4,ne,lx,ly)+&
                todet/det*tt(3,m)*ss(m,n)*tt(3,n)
		enddo
      enddo
	return
      end