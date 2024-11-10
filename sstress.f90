      subroutine sstress
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension x(2,9),toc(3,3),xto(2,2),xot(2,2),&
       xj(2,2),xji(2,2),rs(2),toxj(2,2),toxji(2,2),temss(6)
      dimension xmj(3),xmi(3),dxmj(3,6),ddxmj(3,6,6),&
       obc(6,6),ocuu(6,6),ocup(6),temdis(2,9),temsa(6)
      ig=2
      do ne=1,numele(ig)
         do nos=1,nis(ig)
            ntem=nea(ne,nos)
            do noj=1,2
               y(noj,nos)=coor(ntem,noj)+dis(noj,ntem,1)
               x(noj,nos)=coor(ntem,noj)
            enddo
	   enddo
         if (nstr .eq. 2) then
            do nos=1,nis(ig)
               ntem=nea(ne,nos)
               do noj=1,2
                  temdis(noj,nos)=dis(noj,ntem,1)
               enddo
		  enddo
         endif
         if (nstr .eq. 3) then
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
                  call scauchy(det,todet,xto,lx,ly,ne)
			enddo
		  enddo
         endif
         if (nstr .eq. 2) then
            do i=1,4
               rs(1)=xxrs(i,1)
               rs(2)=xxrs(i,2)
               call element(rs)
               call jacob(x,xj,xji,det)
               call bdpd(xji)
               temsa(1)=0.0d0
               do k=1,nis(ig)
                  temsa(1)=temsa(1)+bd(1,k)*temdis(1,k)
			 enddo
               temsa(2)=0.0d0
               do k=1,nis(ig)
                  temsa(2)=temsa(2)+bd(2,k)*temdis(2,k)
			 enddo
               temsa(3)=0.0d0
               do k=1,nis(ig)
                  temsa(3)=temsa(3)+bd(1,k)*temdis(2,k)+&
                      bd(2,k)*temdis(1,k)
               enddo
               sbulk=syoun/3.0d0/(1-2.0d0*spois)
               if (npss .eq. 1) then
                  cstr(i,ne,1,1)=-sbulk*(temsa(1)+temsa(2))
               endif
               if (npss .eq. 0) then
                  call material
                  do k=1,3
                     temss(k)=0.0d0
                     do m=1,3
                        temss(k)=temss(k)+cmat(k,m)*temsa(m)
                     enddo
				enddo
                  temsa(4)=-1.0d0/syoun*spois*(temss(1)+temss(2))
                  cstr(i,ne,1,1)=sbulk*(temsa(1)+temsa(2)+temsa(4))
               endif
		  enddo
         endif             
      enddo
      return
      end