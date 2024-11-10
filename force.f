	subroutine force
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension x(2,9),xj(2,2),xji(2,2),rs(2)
      dimension xfutem(8),coorc(2,2)
      dimension nloc(4),xfsktem(4,4)
      dimension n1d(2),n1c(2)
      logical ja,jas
      do i=1,nndtem1
         drf(i)=0.0d0
         sk(i)=0.0d0
      enddo
	do i=1,nndim*npart*2
	   drfo(i)=0.0d0
	enddo
      do ig=1,ignum
		nej=0
          do i=1,ig-1
               nej=nej+numele(i)
          enddo
		if (ig .eq. 1) then
			numeltemp=numele(1)
		else
			if (nimmer .eq. 1) then
				numeltemp=numele(2)-numeleim*npart
			else
				numeltemp=numele(2)
			endif
		endif
          do nei=1,numeltemp
			ne=nej+nei
              if ((nstr .eq. 3) .and. (ig .eq. 2)) then
				do nos=1,nis(ig)
					ntem=nea(ne,nos)
					do noj=1,2
						y(noj,nos)=coor(ntem,noj)+dis(noj,ntem,1)
						x(noj,nos)=coor(ntem,noj)
					enddo
				enddo
			else
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
              endif
              do i=1,2*nis(ig)
				xfsk(i)=0.0d0
              enddo
              jas=(ig .eq. 1) .and. 
     $              ((nfrees .eq. 0) .or. (nale .eq. 0))
              if (jas) then
				call eleinit(ne)
                  do lx=1,nint
                     rs(1)=xg(lx,nint)
                     do ly=1,nint
                        rs(2)=xg(ly,nint)
                        call element(rs)
                        call press(rs)
                        call jacob(x,xj,xji,det)
                        call bdpd(xji)
                        wp=wgt(lx,nint)*wgt(ly,nint)
                        w=wp*det
                        call fluid(w,ne)
					enddo
				enddo
                  if ((nsmoo .eq. 5) .or. (ncon .eq. 1)) then
                     if (nis(ig) .eq. 4) then
                        do i=1,4
                           nloc(i)=nea(ne,i)		
		              enddo
		              call con_raw(nloc,xfsktem,xfutem)
                        do m=1,4
                           xfu(m)=xfu(m)-xfutem(m)
                           xfu(m+nis(1))=xfu(m+nis(1))-xfutem(m+4)
                           xfsk(m)=xfsk(m)+xfsktem(m,m)
                           xfsk(m+nis(1))=xfsk(m+nis(1))+xfsktem(m,m)
                        enddo
                     endif
                     if (nis(ig) .eq. 9) then
                        do k=1,4
                           do n=1,4
                              nloc(n)=nea(ne,nloct(k,n))
						 enddo
                           call con_raw(nloc,xfsktem,xfutem)
                           do m=1,4
                              ni=nloct(k,m)
                              xfu(ni)=xfu(ni)-xfutem(m)
                              xfu(ni+nis(1))=xfu(ni+nis(1))-xfutem(m+4)
                              xfsk(ni)=xfsk(ni)+xfsktem(m,m)
                              xfsk(ni+nis(1))=xfsk(ni+nis(1))+
     $							xfsktem(m,m)
						 enddo
	                   enddo
                     endif
                  endif
                  call asmpr(ne)
               endif
               jas=(ig .eq. 1) .and. (nfrees .eq. 1) .and.
     $              (nale .eq. 1) 
               if (jas) then
                  if ((nsmoo .ne. 0) .and. (nsmoo .ne. 5)) then
                     call eleleng(ne)
                  endif
                  call eleinit(ne)
                  do lx=1,nint
                     rs(1)=xg(lx,nint)
                     do ly=1,nint
                        rs(2)=xg(ly,nint)
                        call element(rs)
                        call press(rs)
                        call jacob(x,xj,xji,det)
                        call bdpd(xji)
                        wp=wgt(lx,nint)*wgt(ly,nint)
                        w=wp*det
                        call surfluid(w,ne)
                     enddo
				enddo
                  if ((nsmoo .eq. 5) .or. (ncon .eq. 1)) then
                     if (nis(1) .eq. 4) then
                        do i=1,4
                           nloc(i)=nea(ne,i)
					  enddo
                        call con_raw(nloc,xfsktem,xfutem)
                        do m=1,4
                           xfu(m)=xfu(m)-xfutem(m)
                           xfu(m+nis(1))=xfu(m+nis(1))-xfutem(m+4)
                           xfsk(m)=xfsk(m)+xfsktem(m,m)
                           xfsk(m+nis(1))=xfsk(m+nis(1))+xfsktem(m,m)
                        enddo
                     endif
                     if (nis(1) .eq. 9) then
                        do k=1,4
                           do n=1,4
                              nloc(n)=nea(ne,nloct(k,n))
                           enddo
                           call con_raw(nloc,xfsktem,xfutem)
                           do m=1,4
                              ni=nloct(k,m)
                              xfu(ni)=xfu(ni)-xfutem(m)
                              xfu(ni+nis(1))=xfu(ni+nis(1))-xfutem(m+4)
                              xfsk(ni)=xfsk(ni)+xfsktem(m,m)
                              xfsk(ni+nis(1))=xfsk(ni+nis(1))+
     $							xfsktem(m,m)
						 enddo
	                  enddo
                      endif
                  endif
                  call asmpr(ne)
              endif
              if (ig .eq. 2) then
				if (nis(2) .eq. 2) then
					call strut(w,ne)
                  else
					if (nstr .ne. 3) then
						do lx=1,nint
							rs(1)=xg(lx,nint)
							do ly=1,nint
								rs(2)=xg(ly,nint)
								call element(rs)
								call jacob(x,xj,xji,det)
								call bdpd(xji)
								wp=wgt(lx,nint)*wgt(ly,nint)
								w=wp*det                     
								call strut(w,ne)
							enddo
						enddo
					endif
					if (nstr .eq. 3) then
						call stang(x,ne)
					endif
                  endif
              endif
			if (ig .eq. 1) then
				do ni=1,nis(ig)
					do m=1,2
						nu1=2*(nea(ne,ni)-1)+m
						ne1=(m-1)*nis(ig)+ni
						km=id(nu1)
						if (km .ne. 0) then
							drf(km)=drf(km)+xfu(ne1)
							sk(km)=sk(km)+xfsk(ne1)
						endif
					enddo
				enddo
			endif
			if (ig .eq. 2) then
				do ni=1,nis(ig)
					do m=1,2
						nu1=2*(nea(ne,ni)-1)+m
						ne1=(m-1)*nis(ig)+ni
						km=id(nu1)
						if (km .ne. 0) then
							drf(km)=drf(km)+xfu(ne1)
							sk(km)=sk(km)+
     $							alpha*dt/beta*xfsk(ne1)+
     $							xfsm(ne1)/beta/dt
						endif
					enddo
				enddo
			endif
		enddo
	enddo
	if (nimmer .eq. 1) then
		ig=2
		nej=numele(1)+numele(2)-numeleim*npart
          do nei=1,numeleim*npart
			ne=nej+nei
              if (nstr .eq. 3) then
				do nos=1,nis(ig)
					ntem=nea(ne,nos)
					do noj=1,2
						y(noj,nos)=coor(ntem,noj)+dis(noj,ntem,1)
						x(noj,nos)=coor(ntem,noj)
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
              do i=1,2*nis(ig)
				xfsk(i)=0.0d0
              enddo
			do i=1,nis(ig)
				xfu(i)=0.0d0
                  xfu(i+nis(ig))=0.0d0
                  xfsm(i)=0.0d0
				xfsm(i+nis(ig))=0.0d0
			enddo
			if (nstr .eq. 2) then
				do lx=1,nint
					rs(1)=xg(lx,nint)
					do ly=1,nint
						rs(2)=xg(ly,nint)
						call element(rs)
						call jacob(x,xj,xji,det)
						call bdpd(xji)
						wp=wgt(lx,nint)*wgt(ly,nint)
						w=wp*det                     
						call strut(w,ne)
					enddo
				enddo
			endif
			if (nstr .eq. 3) then
				call stanga(x,ne)
			endif
		enddo
		if (nstr .eq. 4) then
			coorc(1,1)=0.015d0
			coorc(1,2)=0.1d0
			coorc(2,1)=0.025d0
			coorc(2,2)=0.1d0
			do ni=1,npart*nndim
				ncm1=2*(ni-1)+1
				ncm2=2*(ni-1)+2
				nim=nnd-npart*nndim+ni
				r1=dsqrt((coor(nim,1)+dis(1,nim,1)-coorc(1,1))**2+
     $           (coor(nim,2)+dis(2,nim,1)-coorc(1,2))**2)
				tension1=xstif(ni)*(r1-0.005d0)
				r2=dsqrt((coor(nim,1)+dis(1,nim,1)-coorc(2,1))**2+
     $           (coor(nim,2)+dis(2,nim,1)-coorc(2,2))**2)
				tension2=xstif(ni)*(r2-0.005d0) 
				if (r1 .ne. 0.0d0) then
					unitvector(1,1)=(coor(nim,1)+dis(1,nim,1)-
     $					coorc(1,1))/r1
 					unitvector(1,2)=(coor(nim,2)+dis(2,nim,1)-
     $					coorc(1,2))/r1
				endif              
				if (r2 .ne. 0.0d0) then
					unitvector(2,1)=(coor(nim,1)+dis(1,nim,1)-
     $					coorc(2,1))/r2
 					unitvector(2,2)=(coor(nim,2)+dis(2,nim,1)-
     $					coorc(2,2))/r2
				endif              
				drfo(ncm1)=drfo(ncm1)-
     $				xmass(ni)*(-pargravo(1)+acm(1,nim,1))-
     $				tension1*unitvector(1,1)-tension2*unitvector(2,1)
				drfo(ncm2)=drfo(ncm2)-
     $				xmass(ni)*(-pargravo(2)+acm(2,nim,1))-
     $				tension1*unitvector(1,2)-tension2*unitvector(2,2)
			enddo
		endif
		ig=2
		nej=numele(1)+numele(2)-numeleim*npart
          do nei=1,numeleim*npart
			ne=nej+nei
			do nos=1,nis(ig)
				ntem=nea(ne,nos)
				do noj=1,2
					x(noj,nos)=coor(ntem,noj)+dis(noj,ntem,1)
				enddo
			enddo
              do i=1,2*nis(ig)
				xfska(i)=0.0d0
              enddo
			do i=1,nis(ig)
				xfua(i)=0.0d0
				xfua(i+nis(ig))=0.0d0
			enddo
              do lx=1,nint
				rs(1)=xg(lx,nint)
				do ly=1,nint
					rs(2)=xg(ly,nint)
                      call element(rs)
                      call press(rs)
                      call jacob(x,xj,xji,det)
                      call bdpd(xji)
                      wp=wgt(lx,nint)*wgt(ly,nint)
                      w=wp*det
                      call fluida(w,ne)
				enddo
			enddo
			do ni=1,nis(ig)
				do m=1,2
					nu1=2*(nea(ne,ni)-1-(nnd-nndim*npart))+m
					ne1=(m-1)*nis(ig)+ni
					drfo(nu1)=drfo(nu1)-xfua(ne1)
				enddo
			enddo
		enddo
	endif
      if ((nfrees .eq. 1) .and. (nale .eq. 1)) then
		if (nsur .eq. 1) then
			call surdif
          else
              do nes=1,neles
				do lx=1,nint
                     r=xg(lx,nint)
                     call belement(r,det,nes)
                     wp=wgt(lx,nint)
                     w=wp*det
                     call surele(w,nes,det)
	            enddo
                  if ((nsmoos .eq. 5) .or. (ncons .eq. 1)) then
                     if (niss .eq. 2) then
                        n1d(1)=nndtem+neas(nes,1)
                        n1d(2)=nndtem+neas(nes,2)
                        n1c(1)=nnn(neas(nes,1))
                        n1c(2)=nnn(neas(nes,2))
                        det=(coor(n1c(1),1)-
     $                       coor(n1c(2),1))/2.0d0
                        call sur_1d(det,n1d,n1c)
                     endif
                     if (niss .eq. 3) then
                        n1d(1)=nndtem+neas(nes,1)
                        n1d(2)=nndtem+neas(nes,3)
                        n1c(1)=nnn(neas(nes,1))
                        n1c(2)=nnn(neas(nes,3))
                        det=(coor(n1c(1),1)-
     $                       coor(n1c(2),1))/2.0d0
                        call sur_1d(det,n1d,n1c)
                        n1d(1)=nndtem+neas(nes,3)
                        n1d(2)=nndtem+neas(nes,2)
                        n1c(1)=nnn(neas(nes,3))
                        n1c(2)=nnn(neas(nes,2))
                        det=(coor(n1c(1),1)-
     $                       coor(n1c(2),1))/2.0d0
                        call sur_1d(det,n1d,n1c)
					endif
				endif          
			enddo
		endif
      endif
      do ne=1,numeb
		nflag=nface(ne)
          if ((nale .eq. 1).or. (nfsi .gt. 1)) then
			do noj=1,2
				do nos=1,nisa
                     ntem=nea(nbe(ne),nos)
                     x(noj,nos)=coor(ntem,noj)+dis(noj,ntem,1)
				enddo
               enddo
		else
			do noj=1,2
				do nos=1,nisa
                     ntem=nea(nbe(ne),nos)
                     x(noj,nos)=coor(ntem,noj)
				enddo
			enddo
          endif
          do ni=1,3
			do iy=1,nint
				si=xg(iy,nint)
                  wg=wgt(iy,nint)
                  call bload(x,si,nflag,wt,wtt1,wtt2,wg,ni,ne)
                  nee=nea(nbe(ne),nbn(ne,ni))
                  call assbf(wt,wtt1,wtt2,nflag,nee)
			enddo
		enddo
	enddo
	if (nimmer .eq. 1) then
		nvelfor=0
		call delta
	endif
      if (numfn .gt. 0) then
		call nodalf
      endif
      np=nnd
      call skew
      call inter
	return
	end