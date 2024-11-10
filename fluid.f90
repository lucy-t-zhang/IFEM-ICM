      subroutine fluid(w,ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      tem1=w*fdensi
      tem3=w/fbulk
      tem5=fdensi*w/falph
      call acc(ne)
      bpres=0.0d0
      bpresv=0.0d0
	do i=1,nump
		ni=neap(ne,i)
          bpres=bpres+epc(ni,1)*hp(i)
          bpresv=bpresv+epcv(ni,1)*hp(i)
	enddo
	xelep=xdve(1,1)+xdve(2,2)+bpresv/fbulk
	xtp=w*bpres-2.0d0*tem5*bpresv/fbulk/3.0d0
	do i=1,nump
		xkpp(i,ne)=xkpp(i,ne)-tem3*hp(i)*hp(i)/beta/dt
      enddo
      if ((nsmoo .ne. 5) .and. (ncon .ne. 1)) then
		do ni=1,nis(ig)
              xfsk(ni)=xfsk(ni)+&
     			tem1*h(ni)*(defv(1)*bd(1,ni)+defv(2)*bd(2,ni))
              xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+&
     			tem1*h(ni)*(defv(1)*bd(1,ni)+defv(2)*bd(2,ni))
		enddo  
      endif
      do ni=1,nis(ig)
          xfsk(ni)=xfsk(ni)+tem5*(bd(1,ni)*bd(1,ni)+&
     		bd(2,ni)*bd(2,ni))
          xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+&
     		tem5*(bd(1,ni)*bd(1,ni)+bd(2,ni)*bd(2,ni))
      enddo
	do ni=1,nis(ig)
		xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+&
     		h(ni)*h(ni)*tem1/beta/dt
          xfsk(ni)=xfsk(ni)+h(ni)*h(ni)*tem1/beta/dt   
	enddo
      do i=1,nis(ig)
		xfu(i)=xfu(i)-(tem1*h(i)*xac(1)-xtp*bd(1,i)+&
     		tem5*(bd(1,i)*xdve(1,1)+bd(2,i)*xdve(1,2)))
		xfu(i+nis(ig))=xfu(i+nis(ig))-(tem1*h(i)*xac(2)-&
     		xtp*bd(2,i)+tem5*(bd(1,i)*xdve(2,1)+&
     		bd(2,i)*xdve(2,2)))
      enddo
      do i=1,nump
		xfp(i,ne)=xfp(i,ne)-w*hp(i)*xelep
      enddo
      return
      end