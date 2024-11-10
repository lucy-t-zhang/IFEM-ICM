      subroutine fluida(w,ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      tem1=w*fdensi
      tem3=w/fbulk
      tem5=fdensi*w/falph
      call acca(ne)
      bpres=0.0d0
      bpresv=0.0d0
	do i=1,nump
		ni=neap(ne,i)-(nndp-nndaim*npart)
          bpres=bpres+epco(ni,1)*hp(i)
          bpresv=bpresv+epcov(ni,1)*hp(i)
	enddo
	xtp=w*bpres-2.0d0*tem5*bpresv/fbulk/3.0d0
      do ni=1,nis(ig)
          xfska(ni)=xfska(ni)+tem5*(bd(1,ni)*bd(1,ni)+
     $		bd(2,ni)*bd(2,ni))
          xfska(ni+nis(ig))=xfska(ni+nis(ig))+
     $		tem5*(bd(1,ni)*bd(1,ni)+bd(2,ni)*bd(2,ni))
      enddo
	do ni=1,nis(ig)
		xfska(ni+nis(ig))=xfska(ni+nis(ig))+
     $		h(ni)*h(ni)*tem1/beta/dt
          xfska(ni)=xfska(ni)+h(ni)*h(ni)*tem1/beta/dt   
	enddo
      do i=1,nis(ig)
		xfua(i)=xfua(i)-(tem1*h(i)*xac(1)-xtp*bd(1,i)+
     $		tem5*(bd(1,i)*xdve(1,1)+bd(2,i)*xdve(1,2)))
		xfua(i+nis(ig))=xfua(i+nis(ig))-(tem1*h(i)*xac(2)-
     $		xtp*bd(2,i)+tem5*(bd(1,i)*xdve(2,1)+
     $		bd(2,i)*xdve(2,2)))
      enddo
      return
      end