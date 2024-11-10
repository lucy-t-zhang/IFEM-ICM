      subroutine rebound
      implicit real*8 (a-h,o-z)
      include 'common'
      logical ja1
	if ((nfsi .eq. 3) .or. (nfsi .eq. 1)) then
		do i=1,nnd
			do j=1,2
				dui(j,i)=0.0d0
			enddo
		enddo 
		if (nfrees .eq. 1) then
			 do i=1,nnd
				do j=1,2
					duc(j,i)=0.0d0
				enddo
			enddo
		endif
		do in=1,nndtem
			if (idu(in) .gt. 2*nnd) then
				nj=in-2*nnd
				duc(2,nnn(nj))=drf(in)
				velm(2,nnn(nj),1)=velm(2,nnn(nj),1)+drf(in)
			else
				if (mod(idu(in),2) .eq. 0) then
					ntem=idu(in)/2
					ni=2
				endif
				if (mod(idu(in),2) .eq. 1) then
					ntem=idu(in)/2+1
					ni=1
				endif
				dui(ni,ntem)=drf(in)
			endif
		enddo
		do i=1,numct
			dui(islavdir(i),nodesl(i))=amct(i)*&
     			dui(imasdir(i),nodema(i))
		enddo
		call antiskew
		if (nimmer .eq. 1) then
			nndt=nnd-nndim*npart
		else
			nndt=nnd
		endif
         	erru=0.0d0
		do ni=1,2
			do ntem=1,nndt
				erru=erru+dui(ni,ntem)**2
				vel(ni,ntem,1)=vel(ni,ntem,1)+dui(ni,ntem)
				acm(ni,ntem,1)=(-(1.0d0-beta)*acm(ni,ntem,0)+&
     				1.0d0/dt*(vel(ni,ntem,1)-vel(ni,ntem,0)))/beta
			enddo
		enddo
		if (nale .eq. 1) then
			do i=1,nndfsi1
				dum=0.0d0
				do k=1,n1link(i)
					ntem=nod1link(k,i)
					dum=dum+coe1link(k,i)*vel(1,ntem,1)
				enddo
				velm(1,i,1)=dum
				dum=0.0d0
				do k=1,n2link(i)
					ntem=nod2link(k,i)
					if (nfrees .eq. 1) then
						dum=dum+coe2link(k,i)*velm(2,ntem,1)
					else
						dum=dum+coe2link(k,i)*vel(2,ntem,1)
					endif
				enddo
				velm(2,i,1)=dum
			enddo
		endif
		if (nsmoo .eq. 2) then
			call smooth
		endif
		if (nale .eq. 1) then
			do i=1,2
				do ntem=1,nndfsi1
					acmm(i,ntem,1)=(-(1.0d0-beta)*acmm(i,ntem,0)+&
     				1.0d0/dt*(velm(i,ntem,1)-velm(i,ntem,0)))/beta
					dis(i,ntem,1)=dis(i,ntem,0)+&
     					dt*velm(i,ntem,0)+dt**2*&
     					((0.5d0-alpha)*acmm(i,ntem,0)+&
     					alpha*acmm(i,ntem,1))
				enddo
				do j=1,nndfsi2
					ntem=mapnea(j)
					dis(i,ntem,1)=dis(i,ntem,0)+&
     					dt*vel(i,ntem,0)+dt**2*&
     					((0.5d0-alpha)*acm(i,ntem,0)+&
     					alpha*acm(i,ntem,1))
				enddo
			enddo
		else
			do i=1,2
				do j=1,nndfsi2-nndim*npart
					ntem=mapnea(j)
					dis(i,ntem,1)=dis(i,ntem,0)+&
     					dt*vel(i,ntem,0)+dt**2*&
     					((0.5d0-alpha)*acm(i,ntem,0)+&
     					alpha*acm(i,ntem,1))
				enddo
			enddo
		endif
		if (nimmer .eq. 1) then
			nvelfor=1
			call delta
		    do ni=1,2
				do ntem=nndt+1,nnd
					dis(ni,ntem,1)=dis(ni,ntem,0)+&
     					dt*vel(ni,ntem,0)+dt**2*&
     					((0.5d0-alpha)*acm(ni,ntem,0)+&
     					alpha*acm(ni,ntem,1))
				enddo
			enddo

		endif
	endif
	if (nfsi .eq. 3) then
		rnormuf=0.0d0
		do ni=1,2
			do ntem=1,nndfsi1-icount3
				rnormuf=rnormuf+dui(ni,ntem)**2
			enddo
		enddo
		rnormus=0.0d0
		do ni=1,2
			do nte=1,nndfsi2
				ntem=mapnea(nte)
				rnormus=rnormus+dui(ni,ntem)**2
			enddo
		enddo
	endif
	if (nfsi .eq. 2) then
		do i=1,nnd
			do j=1,2
				dui(j,i)=0.0d0
			enddo
		enddo 
		do in=1,nndtem
			if (mod(idu(in),2) .eq. 0) then
				ntem=idu(in)/2
				ni=2
			endif
			if (mod(idu(in),2) .eq. 1) then
				ntem=idu(in)/2+1
				ni=1
			endif
			dui(ni,ntem)=drf(in)
		enddo
		do i=1,numct
			dui(islavdir(i),nodesl(i))=amct(i)*&
     			dui(imasdir(i),nodema(i))
		enddo
		call antiskew
		do ni=1,2
			do ntem=1,nnd
				erru=erru+dui(ni,ntem)**2
				vel(ni,ntem,1)=vel(ni,ntem,1)+dui(ni,ntem)
				acm(ni,ntem,1)=(-(1.0d0-beta)*acm(ni,ntem,0)+&
     				1.0d0/dt*(vel(ni,ntem,1)-vel(ni,ntem,0)))/beta
				dis(ni,ntem,1)=dis(ni,ntem,0)+&
     				dt*vel(ni,ntem,0)+dt**2*&
     				((0.5d0-alpha)*acm(ni,ntem,0)+&
     				alpha*acm(ni,ntem,1))
			enddo
		enddo
	endif
      return
      end