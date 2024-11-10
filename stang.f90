      subroutine stang(x,ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension x(2,9),toc(3,3),xto(2,2),xot(2,2),&
     	xj(2,2),xji(2,2),rs(2),toxj(2,2),toxji(2,2)
      dimension xmj(3),xmi(3),dxmj(3,6),ddxmj(3,6,6),&
     	obc(6,6),ocuu(6,6),ocup(6)
      do lx=1,nint
         rs(1)=xg(lx,nint)
         do ly=1,nint
            rs(2)=xg(ly,nint)
            call element(rs)
            call jacob(y,xj,xji,det)
            call jacob(x,toxj,toxji,todet)
            call bdpd(toxji)
            call stoxc(xto,xot,xj,xji,toxj,toxji,toc)
            call smaterj(wto,toc,xmi,xmj,dxmj,ddxmj)
            call spress(rs,ne)
            call sbpress(dxmj,ddxmj,xmj)
            call sboc(obc,ocpp,ocuu,ocup,xmj,dxmj,ddxmj)
            call sstrain(toc,xto,lx,ly,ne)
            call spiola(ocpp,xmj,dxmj)
            wp=wgt(lx,nint)*wgt(ly,nint)*thic 
            w=wp*todet
            call sstif(ocpp,ocuu,ocup,ne,w,toxj)
		enddo
	enddo
      return
      end