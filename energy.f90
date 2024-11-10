      subroutine energy(xene,xvol)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension rs(2),x(2,9),wtem(2),xtem(2)
      dimension xj(2,2),xji(2,2)
      xene=0.0d0
      xvol=0.0d0
      do ig=1,ignum
         if (ig .eq. 1) then
            dd=fdensi
         else
            dd=sdensi
         endif
         do ne=1,numele(ig)
            if (nale .eq. 1) then
               do nos=1,nis(ig)
                  ntem=nea(ne,nos)
                  do noj=1,2
                     x(noj,nos)=coor(ntem,noj)+dis(noj,ntem,1)
                  enddo
			 enddo
            else
               do nos=1,nis(ig)
                  ntem=nea(ne,nos)
                  do noj=1,2
                     x(noj,nos)=coor(ntem,noj)
                  enddo
               enddo
            endif     
            do lx=1,nint
               rs(1)=xg(lx,nint)         
               do ly=1,nint
                  rs(2)=xg(ly,nint)
                  call element(rs)     
                  call jacob(x,xj,xji,det)
                  wp=wgt(lx,nint)*wgt(ly,nint)
                  w=wp*det
                  xvel(1)=0.0d0
                  xvel(2)=0.0d0
                  do k=1,nis(ig)
                     ntem=nea(ne,k)
                     do i=1,2
                        xvel(i)=xvel(i)+h(k)*vel(i,ntem,1)
                     enddo
				enddo
                  xene=xene+(xvel(1)**2+xvel(2)**2)*dd*0.5d0*w
               enddo
	      enddo
         enddo
      enddo
      if (nfrees .eq. 1) then
         wtem(1)=1.0d0
         wtem(2)=1.0d0
         xtem(1)=-1.0d0/dsqrt(3.0d0)
         xtem(2)=1.0d0/dsqrt(3.0d0)
         do nes=1,neles
            do lx=1,2
               r=xtem(lx)
               w=wtem(lx)
               ttem=gravit*w*fdensi
               call belement(r,det,nes)
               ydis=0.0d0
               do k=1,niss
                  ntem=nfsnd(nes,k)
                  ydis=ydis+hfs(k)*dis(2,ntem,1)
               enddo
               xene=xene+0.5d0*ttem*det*ydis**2
               xvol=xvol+w*det*ydis
            enddo
		enddo
      endif
      return
      end