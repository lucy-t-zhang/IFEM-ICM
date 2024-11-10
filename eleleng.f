      subroutine eleleng(ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension xt(2),xs(2),xy(2,9),xji(2,2)
      dimension rs(2),xj(2,2),u(2),sn(2)
      do nos=1,nis(ig)
         ntem=nea(ne,nos)
         do noj=1,2
            xy(noj,nos)=coor(ntem,noj)+dis(noj,ntem,0)
		enddo
	enddo
      if (nis(ig) .eq. 4) then
         rs(1)=1.0d0
         rs(2)=0.0d0
         call element(rs)
         xt(1)=0.0d0
         xt(2)=0.0d0
         do k=1,nis(ig)
            do i=1,2
               xt(i)=xt(i)+h(k)*xy(i,k)
		  enddo
	   enddo
         rs(1)=-1.0d0
         rs(2)=0.0d0
         call element(rs)
         xs(1)=0.0d0
         xs(2)=0.0d0
         do k=1,nis(ig)
            do i=1,2
               xs(i)=xs(i)+h(k)*xy(i,k)
            enddo
	   enddo
         hrs(1)=dsqrt((xt(2)-xs(2))**2+
     $        (xt(1)-xs(1))**2)
         rs(1)=0.0d0
         rs(2)=1.0d0
         call element(rs)
         xt(1)=0.0d0
         xt(2)=0.0d0
         do k=1,nis(ig)
            do i=1,2
               xt(i)=xt(i)+h(k)*xy(i,k)
            enddo
	   enddo
         rs(1)=0.0d0
         rs(2)=-1.0d0
         call element(rs)
         xs(1)=0.0d0
         xs(2)=0.0d0
         do k=1,nis(ig)
            do i=1,2
               xs(i)=xs(i)+h(k)*xy(i,k)
            enddo
	   enddo
         hrs(2)=dsqrt((xt(2)-xs(2))**2+
     $        (xt(1)-xs(1))**2)
         rs(1)=0.0d0
         rs(2)=0.0d0
         u(1)=0.0d0
         u(2)=0.0d0
         call element(rs)
         call jacob(xy,xj,xji,det)
         call bdpd(xji)
         do k=1,nis(ig)
            ntem=nea(ne,k)
            do i=1,2
               u(i)=u(i)+h(k)*(vel(i,ntem,0)-velm(i,ntem,0))
            enddo
	   enddo
         xx1=dsqrt(xj(1,1)**2+xj(1,2)**2)
         xx2=dsqrt(xj(2,1)**2+xj(2,2)**2)
         uw(1)=(xj(1,1)*u(1)+xj(1,2)*u(2))/xx1
         uw(2)=(xj(2,1)*u(1)+xj(2,2)*u(2))/xx2
      endif
      if (nis(ig) .eq. 9) then
         ntem1=nea(ne,5)
         xt(1)=coor(ntem1,1)+dis(1,ntem1,0)
         xt(2)=coor(ntem1,2)+dis(2,ntem1,0)
         ntem1=nea(ne,7)
         xs(1)=coor(ntem1,1)+dis(1,ntem1,0)
         xs(2)=coor(ntem1,2)+dis(2,ntem1,0)
         hrs(2)=dsqrt((xt(2)-xs(2))**2+
     $        (xt(1)-xs(1))**2)
         ntem1=nea(ne,6)
         xt(1)=coor(ntem1,1)+dis(1,ntem1,0)
         xt(2)=coor(ntem1,2)+dis(2,ntem1,0)
         ntem1=nea(ne,8)
         xs(1)=coor(ntem1,1)+dis(1,ntem1,0)
         xs(2)=coor(ntem1,2)+dis(2,ntem1,0)
         hrs(1)=dsqrt((xt(2)-xs(2))**2+
     $        (xt(1)-xs(1))**2)
         rs(1)=0.0d0
         rs(2)=0.0d0
         call element(rs)
         call jacob(xy,xj,xji,det)
         call bdpd(xji)
         ntem=nea(ne,9)
         xx1=dsqrt(xj(1,1)**2+xj(1,2)**2)
         xx2=dsqrt(xj(2,1)**2+xj(2,2)**2)
         u(1)=(vel(1,ntem,0)-velm(1,ntem,0))
         u(2)=(vel(2,ntem,0)-velm(2,ntem,0))
         uw(1)=(xj(1,1)*u(1)+xj(1,2)*u(2))/xx1
         uw(2)=(xj(2,1)*u(1)+xj(2,2)*u(2))/xx2
      endif
      if (nsmoo .eq. 3) then
         do m=1,2
            uwa(m,1)=hrs(m)/6.0d0
            uwa(m,2)=hrs(m)**2/36.0d0
         enddo
      endif
      if (nsmoo .eq. 1) then
         do i=1,2
            if (uw(i) .gt. 0.0d0) then
               sn(i)=1.0d0
            else 
               sn(i)=-1.0d0
            endif
         enddo
         cmod=u(1)**2+u(2)**2
         bkx=0.0d0
         do i=1,2
            bkx=bkx+sn(i)*hrs(i)*uw(i)/dsqrt(15.0d0)
         enddo
         if (cmod .gt. 1.0d-15) then
            upw=bkx/cmod
         else
            upw=0.0d0
         endif
      endif
      if (nsmoo .eq. 4) then
         xv=dsqrt(u(1)**2+u(2)**2)
         xh=dsqrt(hrs(1)**2+hrs(2)**2)
         if (niso .eq. 1) then
            xuk(1,1)=xco*xv*xh
         else
            do m=1,2
               do n=1,2
                  if (xv .gt. 1.0d-15) then
                     xuk(m,n)=xco*xh*u(m)*u(n)/xv
                  else
                     xuk(m,n)=0.0d0
                  endif
               enddo
	      enddo
         endif
      endif
      return
      end