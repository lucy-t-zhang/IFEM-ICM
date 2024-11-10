      subroutine delta      
      implicit real*8 (a-h,o-z)
      include 'common'            
	dimension xcent(2),xcenta(2)
	dimension xpp(3,3),xpj(3,4),xps(3),xpsa(3)
	dimension ntemp(2),ntempa(2)
c
c	nvelfor 1 interpolation of velocity/pressure
c	nvelfor 0 distribution of force
c
	if (nvelfor .eq. 1) then
		do i=1,nndim*npart
			ni=nnd-nndim*npart+i
			do j=1,2
				du(j,ni)=0.0d0
			enddo
		enddo
		do i=1,nndaim*npart
			epcoi(i)=0.0d0
		enddo	 
  		do i=1,nndim*npart
			ni=nnd-nndim*npart+i
			xps(1)=1.0d0
			xcent(1)=coor(ni,1)+dis(1,ni,1)
			xcent(2)=coor(ni,2)+dis(2,ni,1)
			xps(2)=xcent(1)
			xps(3)=xcent(2)
			ntemp(1)=dint(xcent(1)/xinc)+1
			ntemp(2)=dint(xcent(2)/yinc)+1
			ntot=0
			if (ndelta .eq. 1) then
				do m=-2,4
					do n=-2,4 
					do k=1,nmapp(ntemp(1)+m,ntemp(2)+n)
						nt1=npoints(ntemp(1)+m,ntemp(2)+n,k)
						xdis1=dabs(coor(nt1,1)-xcent(1))
      					xdis2=dabs(coor(nt1,2)-xcent(2))
						xrad=dsqrt(xdis1**2+xdis2**2)
						xrad=xrad/xrf
						if (xrad .le. 2.0d0) then
							ntot=ntot+1
                              weip(ntot)=0.25d0*(1+dcos(pi*xrad/2.0d0))
							npp(ntot)=nt1
							xefg(1,ntot)=1.0d0
							xefg(2,ntot)=coor(nt1,1)
							xefg(3,ntot)=coor(nt1,2)
						endif								
					enddo	
				enddo
				enddo
			endif
			if (ndelta .eq. 2) then
				do m=-2,4
					do n=-2,4 
						do k=1,nmapp(ntemp(1)+m,ntemp(2)+n)
							nt1=npoints(ntemp(1)+m,ntemp(2)+n,k)
							xdis1=coor(nt1,1)-xcent(1)
      						xdis2=coor(nt1,2)-xcent(2)
							xrad1=xdis1/xinc
							xrad2=xdis2/yinc
							if ((dabs(xrad1) .le. 2.0d0) .and. 
     $						(dabs(xrad2) .le. 2.0d0)) then
								ntot=ntot+1
						if((xrad1.ge.two1).and.(xrad1.lt.one1)) then
					awx   = (1.0d0/6.0d0)*(2.0d0 + xrad1)**3
					elseif ((xrad1.ge.one1).and.(xrad1.lt.zero)) then
					awx   =   2.0d0/3.0d0 - xrad1**2 - 0.50d0*xrad1**3
					elseif ((xrad1.ge.zero).and.(xrad1.lt.one2)) then
					awx   =   2.0d0/3.0d0 - xrad1**2 + 0.50d0*xrad1**3
					elseif ((xrad1.ge.one2).and.(xrad1.le.two2)) then
					awx   =   (1.0d0/6.0d0)*(2.0d0 - xrad1)**3
					else
					awx   = 0.00d0
					endif     
					if((xrad2.ge.two1).and.(xrad2.lt.one1)) then
					awy   = (1.0d0/6.0d0)*(2.0d0 + xrad2)**3
					elseif ((xrad2.ge.one1).and.(xrad2.lt.zero)) then
					awy   =   2.0d0/3.0d0 - xrad2**2 - 0.50d0*xrad2**3
					elseif ((xrad2.ge.zero).and.(xrad2.lt.one2)) then
					awy   =   2.0d0/3.0d0 - xrad2**2 + 0.50d0*xrad2**3
					elseif ((xrad2.ge.one2).and.(xrad2.le.two2)) then
					awy   =   (1.0d0/6.0d0)*(2.0d0 - xrad2)**3.
					else
					awy   = 0.00d0
					endif   
					weip(ntot)=awx*awy/xinc/yinc 
							npp(ntot)=nt1
							xefg(1,ntot)=1.0d0
							xefg(2,ntot)=coor(nt1,1)
							xefg(3,ntot)=coor(nt1,2)
					endif								
				enddo	
				enddo
				enddo
			endif
			if (ntot .ge. 1) then
				do l=1,3
					xpj(l,1)=0.0d0
					xpj(l,2)=0.0d0
				enddo
				do j=1,3
					do k=1,3
						xpp(j,k)=0.0d0
					enddo
				enddo				
				do m=1,ntot
					nt=npp(m)
					do j=1,3
					xpj(j,1)=xpj(j,1)+weip(m)*dui(1,nt)*xefg(j,m)
					xpj(j,2)=xpj(j,2)+weip(m)*dui(2,nt)*xefg(j,m)
						do k=1,3
							xpp(j,k)=xpp(j,k)+
     $						xefg(k,m)*xefg(j,m)*weip(m)
						enddo
					enddo
				enddo
				call gaussj(xpp,3,3,xpj,2,2)
				do n=1,3
					du(1,ni)=du(1,ni)+xps(n)*xpj(n,1)
					du(2,ni)=du(2,ni)+xps(n)*xpj(n,2)
				enddo
				dui(1,ni)=du(1,ni)
				dui(2,ni)=du(2,ni)
				vel(1,ni,1)=dui(1,ni)+vel(1,ni,1)
				vel(2,ni,1)=dui(2,ni)+vel(2,ni,1)
				acm(1,ni,1)=(-(1.0d0-beta)*acm(1,ni,0)+
     $				1.0d0/dt*(vel(1,ni,1)-vel(1,ni,0)))/beta
				acm(2,ni,1)=(-(1.0d0-beta)*acm(2,ni,0)+
     $				1.0d0/dt*(vel(2,ni,1)-vel(2,ni,0)))/beta
			endif
		enddo
		do i=1,nndaim*npart
			ni=nndp-nndaim*npart+i
			xpsa(1)=1.0d0
			xcenta(1)=coor(nmap(ni),1)+dis(1,nmap(ni),1)
			xcenta(2)=coor(nmap(ni),2)+dis(2,nmap(ni),1)
			xpsa(2)=xcenta(1)
			xpsa(3)=xcenta(2)
			ntempa(1)=dint(xcenta(1)/xinc)+1
			ntempa(2)=dint(xcenta(2)/yinc)+1
			ntot=0
			if (ndelta .eq. 1) then
				do m=-2,4
				do n=-2,4 
					do k=1,nmappa(ntempa(1)+m,ntempa(2)+n)
						nt1=npointsa(ntempa(1)+m,ntempa(2)+n,k)
						xdis1=dabs(coorp(nt1,1)-xcenta(1))
      					xdis2=dabs(coorp(nt1,2)-xcenta(2))
						xrad=dsqrt(xdis1**2+xdis2**2)
						xrad=xrad/xrf
						if (xrad .le. 2.0d0) then
							ntot=ntot+1
							weipa(ntot)=0.25d0*(1+dcos(pi*xrad/2.0d0))
							nppa(ntot)=nt1	
							xefga(1,ntot)=1.0d0
							xefga(2,ntot)=coorp(nt1,1)
							xefga(3,ntot)=coorp(nt1,2)
						endif								
					enddo	
				enddo
				enddo
			endif
			if (ndelta .eq. 2) then
				do m=-2,4
					do n=-2,4 
						do k=1,nmappa(ntempa(1)+m,ntempa(2)+n)
						nt1=npointsa(ntempa(1)+m,ntempa(2)+n,k)
						xdis1=coorp(nt1,1)-xcenta(1)
      					xdis2=coorp(nt1,2)-xcenta(2)
						xrad1=xdis1/xinc
						xrad2=xdis2/yinc
						if ((dabs(xrad1) .le. 2.0d0) .and. 
     $						(dabs(xrad2) .le. 2.0d0)) then
							ntot=ntot+1
				if((xrad1.ge.two1).and.(xrad1.lt.one1)) then
				awx   = (1.0d0/6.0d0)*(2.0d0 + xrad1)**3
				elseif ((xrad1.ge.one1).and.(xrad1.lt.zero)) then
				awx   =   2.0d0/3.0d0 - xrad1**2 - 0.50d0*xrad1**3
				elseif ((xrad1.ge.zero).and.(xrad1.lt.one2)) then
				awx   =   2.0d0/3.0d0 - xrad1**2 + 0.50d0*xrad1**3
				elseif ((xrad1.ge.one2).and.(xrad1.le.two2)) then
				awx   =   (1.0d0/6.0d0)*(2.0d0 - xrad1)**3
				else
				awx   = 0.00d0
				endif     
				if((xrad2.ge.two1).and.(xrad2.lt.one1)) then
				awy   = (1.0d0/6.0d0)*(2.0d0 + xrad2)**3
				elseif ((xrad2.ge.one1).and.(xrad2.lt.zero)) then
				awy   =   2.0d0/3.0d0 - xrad2**2 - 0.50d0*xrad2**3
				elseif ((xrad2.ge.zero).and.(xrad2.lt.one2)) then
				awy   =   2.0d0/3.0d0 - xrad2**2 + 0.50d0*xrad2**3
				elseif ((xrad2.ge.one2).and.(xrad2.le.two2)) then
				awy   =   (1.0d0/6.0d0)*(2.0d0 - xrad2)**3.
				else
				awy   = 0.00d0
				endif   
				weipa(ntot)=awx*awy/xinc/yinc 
							nppa(ntot)=nt1	
							xefga(1,ntot)=1.0d0
							xefga(2,ntot)=coorp(nt1,1)
							xefga(3,ntot)=coorp(nt1,2)
						endif								
					enddo	
				enddo
				enddo
			endif
			if (ntot .ge. 1) then
				do l=1,3
					xpj(l,1)=0.0d0
					xpj(l,2)=0.0d0
				enddo
				do j=1,3
					do k=1,3
						xpp(j,k)=0.0d0
					enddo
				enddo				
				do m=1,ntot
					nt=nppa(m)
					do j=1,3
					xpj(j,1)=xpj(j,1)+weipa(m)*epci(nt)*xefga(j,m)
						do k=1,3
							xpp(j,k)=xpp(j,k)+
     $						xefga(k,m)*xefga(j,m)*weipa(m)
						enddo
					enddo
				enddo
				call gaussj(xpp,3,3,xpj,1,1)
				do n=1,3
					epcoi(i)=epcoi(i)+xpsa(n)*xpj(n,1)
				enddo
			epco(i,1)=epco(i,1)+epcoi(i)
			epcov(i,1)=epcov(i,1)+epcoi(i)/beta/dt
			endif
		enddo
	endif	
	if (nvelfor .eq. 0) then
		do i=1,nndim*npart
			ni=nnd-nndim*npart+i
			xps(1)=1.0d0
			xcent(1)=coor(ni,1)+dis(1,ni,1)
			xcent(2)=coor(ni,2)+dis(2,ni,1)
			xps(2)=xcent(1)
			xps(3)=xcent(2)
			ntemp(1)=dint(xcent(1)/xinc)+1
			ntemp(2)=dint(xcent(2)/yinc)+1
			ntot=0
			if (ndelta .eq. 1) then
				do m=-2,4
				do n=-2,4
					do k=1,nmapp(ntemp(1)+m,ntemp(2)+n)
						nt1=npoints(ntemp(1)+m,ntemp(2)+n,k)
						xdis1=dabs(coor(nt1,1)-xcent(1))
      					xdis2=dabs(coor(nt1,2)-xcent(2))
						xrad=dsqrt(xdis1**2+xdis2**2)
						xrad=xrad/xrf
						if (xrad .le. 2.0d0) then
							ntot=ntot+1
                              weip(ntot)=0.25d0*(1+dcos(pi*xrad/2.0d0))
							npp(ntot)=nt1	
							xefg(1,ntot)=1.0d0
							xefg(2,ntot)=coor(nt1,1)
							xefg(3,ntot)=coor(nt1,2)
						endif								
					enddo	
				enddo
				enddo
			endif
			if (ndelta .eq. 2) then
				do m=-2,4
					do n=-2,4 
						do k=1,nmapp(ntemp(1)+m,ntemp(2)+n)
						nt1=npoints(ntemp(1)+m,ntemp(2)+n,k)
						nt1=npoints(ntemp(1)+m,ntemp(2)+n,k)
						xdis1=coor(nt1,1)-xcent(1)
      					xdis2=coor(nt1,2)-xcent(2)
						xrad1=xdis1/xinc
						xrad2=xdis2/yinc
						if ((dabs(xrad1) .le. 2.0d0) .and. 
     $						(dabs(xrad2) .le. 2.0d0)) then
							ntot=ntot+1
				if((xrad1.ge.two1).and.(xrad1.lt.one1)) then
				awx   = (1.0d0/6.0d0)*(2.0d0 + xrad1)**3
				elseif ((xrad1.ge.one1).and.(xrad1.lt.zero)) then
				awx   =   2.0d0/3.0d0 - xrad1**2 - 0.50d0*xrad1**3
				elseif ((xrad1.ge.zero).and.(xrad1.lt.one2)) then
				awx   =   2.0d0/3.0d0 - xrad1**2 + 0.50d0*xrad1**3
				elseif ((xrad1.ge.one2).and.(xrad1.le.two2)) then
				awx   =   (1.0d0/6.0d0)*(2.0d0 - xrad1)**3
				else
				awx   = 0.0d0
				endif     
				if((xrad2.ge.two1).and.(xrad2.lt.one1)) then
				awy   = (1.0d0/6.0d0)*(2.0d0 + xrad2)**3
				elseif ((xrad2.ge.one1).and.(xrad2.lt.zero)) then
				awy   =   2.0d0/3.0d0 - xrad2**2 - 0.50d0*xrad2**3
				elseif ((xrad2.ge.zero).and.(xrad2.lt.one2)) then
				awy   =   2.0d0/3.0d0 - xrad2**2 + 0.50d0*xrad2**3
				elseif ((xrad2.ge.one2).and.(xrad2.le.two2)) then
				awy   =   (1.0d0/6.0d0)*(2.0d0 - xrad2)**3
				else
				awy   = 0.0d0
				endif   
				weip(ntot)=awx*awy/xinc/yinc 
							npp(ntot)=nt1
							xefg(1,ntot)=1.0d0
							xefg(2,ntot)=coor(nt1,1)
							xefg(3,ntot)=coor(nt1,2)
						endif								
					enddo	
				enddo
				enddo
			endif
			if (ntot .ge. 1) then
				do j=1,3
					do k=1,3
						xpp(j,k)=0.0d0
					enddo
				enddo	
				do m=1,ntot
					nt=npp(m)
					do j=1,3
						do k=1,3
							xpp(j,k)=xpp(j,k)+
     $						xefg(k,m)*xefg(j,m)*weip(m)
						enddo
					enddo
				enddo							
				do m=1,ntot
					nt=npp(m)
		            do j=1,3
						xpm(j,m)=0.0d0
					enddo
					do j=1,3
						xpm(j,m)=xpm(j,m)+weip(m)*xefg(j,m)
					enddo
				enddo
				call gaussj(xpp,3,3,xpm,ntot,ntot)
				nv1=2*(i-1)+1
				nu1=2*(i-1)+2
				do m=1,ntot
					xtemp=0.0d0
					do n=1,3
						xtemp=xtemp+xps(n)*xpm(n,m)
					enddo				
					nf1=id(2*(npp(m)-1)+1)
					nf2=id(2*(npp(m)-1)+2)
					drf(nf1)=drf(nf1)+xtemp*drfo(nv1)
					drf(nf2)=drf(nf2)+xtemp*drfo(nu1)
				enddo
			endif
		enddo	
	endif
      return
      end 