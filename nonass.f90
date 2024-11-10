      subroutine nonass(nonfl)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension x(2,9),xj(2,2),xji(2,2),rs(2)
      dimension xfutem(8)
      dimension nloc(4),xfsktem(4,4)
      dimension n1d(2),n1c(2)
      logical ja,jas
      nai=0
      na=0
      if (nrestart .eq. 1) then
         read(25,*) itin
         itin=itin+1
      endif
      if (nrestart .eq. 0) then
         itin=1
      endif
	iti=0
	call timefun
      if (ntprint .eq. 1) then
		do i=1,numtfun
			write(4,130) iti+1,i,tfun(i)
 130			format(1x, 'tfun',i1,'(',i4,')=',e23.7,';'/)
		enddo
		write(4,131) iti+1,0.0d0
 131		format(1x, 'time(',i4,')=',e23.7,';'/)
      endif
	if (nisa .eq. 9) then
		call printnine
	else
		call print
	endif
      do iti=itin,nts
		call timefun
		call load
		write(*,102)
 102		format(1x,'iteration dis. increment  pre.',&
     		&'increment   lambda increment'/,&
                &i3,3x,3(e23.16,2x))
		if (initdir .eq. 2) then
			do k=1,numgb
				if (ndirgb(k) .eq. 121111) then
					do m=1,numdir(k) 
					vel(1,nodegb(k,m),1)=tfun(1)*xinvel(1,nodegb(k,m))              
					enddo
				endif
				if (ndirgb(k) .eq. 112111) then
					do m=1,numdir(k) 
					vel(2,nodegb(k,m),1)=tfun(1)*xinvel(2,nodegb(k,m))              
					enddo
				endif
			enddo
		endif
		if (initdir .eq. 1) then
			do k=1,numgb
				if (ndirgb(k) .eq. 121111) then
					do m=1,numdir(k) 
					dis(1,nodegb(k,m),1)=tfun(1)*xinvel(1,nodegb(k,m))              
					enddo
				endif
				if (ndirgb(k) .eq. 112111) then
					do m=1,numdir(k) 
					dis(2,nodegb(k,m),1)=tfun(1)*xinvel(2,nodegb(k,m))              
					enddo
				endif
			enddo
		endif
		do i=1,nnd
			do j=1,2
				vel(j,i,0)=vel(j,i,1)
				acm(j,i,0)=acm(j,i,1)
			enddo
		enddo
		if (nale .eq. 1) then
			do i=1,nndfsi1
				do j=1,2
					velm(j,i,0)=velm(j,i,1)
					acmm(j,i,0)=acmm(j,i,1)
				enddo
			enddo
		endif
		if (nfsi .ge. 2) then
			do i=1,nnd
				do j=1,2
					dis(j,i,0)=dis(j,i,1)
				enddo
			enddo
		endif
		write(*,101) iti,(tfun(i),i=1,numtfun)
 101		format(2x,'time step'/,i5,2x,e23.16,2x,e23.16)
		do ni=1,nndp
			epc(ni,0)=epc(ni,1)
			epcv(ni,0)=epcv(ni,1)
		enddo
	    do i=1,nndaim*npart
		    epco(i,0)=epco(i,1)
		    epcov(i,0)=epcov(i,1)
	    enddo
		if (nisa .eq. 9) then
			call printnine
		else
			call print
		endif
		ii=0
 1013		errp=0.0d0
		erru=0.0d0
		call force
		do i=1,ndim
			vloc(i)=drf(i)
		enddo 
		do i=1,nnd
			do j=1,2
				vel(j,i,2)=vel(j,i,1)
				acm(j,i,2)=acm(j,i,1)
			enddo
		enddo
          if (nfsi .ge. 2) then
			do i=1,nnd
				do j=1,2
					dis(j,i,2)=dis(j,i,1)
				enddo
			enddo
			if (nale .eq. 1) then
				do i=1,nndfsi1
					do j=1,2
						velm(j,i,2)=velm(j,i,1)
						acmm(j,i,2)=acmm(j,i,1)
					enddo
				enddo
			endif
          endif
		do ni=1,nndp
			epc(ni,2)=epc(ni,1)
			epcv(ni,2)=epcv(ni,1)
		enddo
	    do i=1,nndaim*npart
		    epco(i,2)=epco(i,1)
		    epcov(i,2)=epcov(i,1)
	    enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Newton-Krylov (GMRES-Newton Matrix Free)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		eps = 1.0d-16
		epp = 1.0d-6
	    ept = 1.0d-6
		ein = dsqrt(eps)
		if (ii .eq. 1) then
			xnorm1=0.0d0
			do i=1,ndim1
				xnorm1=xnorm1+drf(i)**2
			enddo
			xnorm1=dsqrt(xnorm1/ndim1)
			xnorm2=0.0d0
			do i=1+ndim1,ndim
				xnorm2=xnorm2+drf(i)**2
			enddo
			xnorm2=dsqrt(xnorm2/ndim2)
			erru0=xnorm1
			errp0=xnorm2
		else
			erru0=1.0d0
			errp0=1.0d0
		endif
		do i=1,ndim1
			smg(i)=sk(i)
		enddo
		do i=ndim1+1,ndim
			smg(i)=sk(i)
		enddo
		iscaling=3
		if(iscaling.eq.0) then
			do i=1,ndim  
				wgres(i)=1.0d0
			enddo
		else if(iscaling.eq.1) then
			do i=1,ndim
				if (wgres(i).lt.0.0d0) smg(i)=-1.0d0
				wgres(i)=1.0d0/dsqrt(dabs(wgres(i)))
			enddo
		else if (iscaling.eq.3) then
			do i=1,ndim
				wgres(i)=1.0d0/smg(i)
			enddo
		else if (iscaling.eq.2) then
			do i=1,ndim
				if (wgres(i).lt.0.0d0) smg(i)=-1.0d0
				wgres(i)=1.0d0/dsqrt(dabs(wgres(i)))
			enddo
		endif
		do i=1,ndim
			dgres(i)=0.0d0
			do j=1,inner+1
				vgres(i,j)=0.0d0
			enddo
		enddo
		do i=1,ndim
			vgres(i,1)=drf(i)
			vgres(i,1)=wgres(i)*vgres(i,1)
		enddo
		igmres = 0 
 741		igmres = igmres+1
		rnormtt=0.0d0
		do i=1,ndim
			rnormtt=rnormtt+vgres(i,1)**2
		enddo
		rnormtt=dsqrt(rnormtt)
		print *, '  gmres', igmres,rnormtt
		if (rnormtt .le. epp .or. igmres .gt. iouter) goto 740
		do i=1,ndim
			vgres(i,1)=vgres(i,1)/rnormtt
		enddo
		do j=1,inner
			do i=1,ndim
			  zgres(i,j)=vgres(i,j)
			enddo
			do i=1,ndim
			  drf(i)=zgres(i,j)
			enddo
			do i=1,ndim
				drf(i)=ein*zgres(i,j)
			enddo
			do i=1,ndim2
			ik=idu(i+nndtem)-nndtem
			ik=idpu(ik)
			epci(ik)=drf(nndtem+i)
			epc(ik,1)=epc(ik,1)+epci(ik)
			epcv(ik,1)=epcv(ik,1)+epci(ik)/beta/dt
      		enddo
			call rebound
			call force
			do i=1,nnd
				do jj=1,2
					vel(jj,i,1)=vel(jj,i,2)
					acm(jj,i,1)=acm(jj,i,2)
				enddo
			enddo
			if (nfsi .ge. 2) then
				do i=1,nnd
					do jj=1,2
						dis(jj,i,1)=dis(jj,i,2)
					enddo
				enddo
				if (nale .eq. 1) then
					do i=1,nndfsi1
						do jj=1,2
							acmm(jj,i,1)=acmm(jj,i,2)
							velm(jj,i,1)=velm(jj,i,2)
						enddo
					enddo
				endif
			endif
			do ni=1,nndp
				epc(ni,1)=epc(ni,2)
				epcv(ni,1)=epcv(ni,2)
			enddo
	        do i=1,nndaim*npart
		    epco(i,1)=epco(i,2)
		    epcov(i,1)=epcov(i,2)
	        enddo
 			do i=1,ndim
				avloc(i)=-(drf(i)-vloc(i))/ein
			enddo
			do i=1,ndim
				avg(i)=avloc(i)
			enddo
     			do i=1,ndim
				avg(i)=wgres(i)*avg(i)
				vgres(i,j+1)=avg(i)
			enddo
			do i=1,j
				tmpo=0.0d0
				do k=1,ndim
					tmpo=tmpo+vgres(k,j+1)*vgres(k,i)
				enddo
				hgres(i,j)=tmpo
				do k=1,ndim
					vgres(k,j+1)=vgres(k,j+1)-tmpo*vgres(k,i) 
				enddo              
			enddo
			tmpo=0.0d0
			do i=1,ndim
				tmpo=tmpo+vgres(i,j+1)*vgres(i,j+1)
			enddo
			tmpo=dsqrt(tmpo)
			hgres(j+1,j)=tmpo
			do i=1,ndim
				vgres(i,j+1)=vgres(i,j+1)/tmpo
			enddo
		enddo		
		ygres(1)=rnormtt
   	    do i=2,inner+1
			ygres(i)=0.0d0
		enddo
	    do j=1,inner
			j1=j+1
			do i=2,j
				i1=i-1
				hsave=hgres(i1,j)
				hgres(i1,j)=+cgres(i1)*hsave+sgres(i1)*hgres(i,j)
				hgres(i,j)=-sgres(i1)*hsave+cgres(i1)*hgres(i,j)
			enddo
			gam=dsqrt(hgres(j,j)**2+hgres(j1,j)**2)
			cgres(j)=hgres(j,j)/gam
			sgres(j)=hgres(j1,j)/gam
			hgres(j,j)=cgres(j)*hgres(j,j)+sgres(j)*hgres(j1,j)
			ygres(j1)=-sgres(j)*ygres(j)
			ygres(j)=+cgres(j)*ygres(j)
	    enddo
		j=inner	
		ygres(j)=ygres(j)/hgres(j,j)
		do jj=2,j
			k = j - jj + 1
			k1 = k + 1
			ysave=ygres(k)
			do l=k1,j
				ysave=ysave-hgres(k,l)*ygres(l)
			enddo
			ygres(k)=ysave/hgres(k,k)
		enddo
		j=inner	
		do jj=1,j
			tmpo=ygres(jj)
			do i=1,ndim
				dgres(i)=dgres(i)+tmpo*zgres(i,jj)
			enddo
		enddo
		if (igmres .le. iouter) then
			do jj=1,j
				ij = j - jj + 2
				ygres(ij-1)=-sgres(ij-1)*ygres(ij)
				ygres(ij)=cgres(ij-1)*ygres(ij)
			enddo
			do jj=1,j+1
				tmpo=ygres(jj)
				if (jj.eq.1) tmpo=tmpo-1.0d0
				do i=1,ndim
					vgres(i,1)=vgres(i,1)+tmpo*vgres(i,jj)
				enddo
			enddo
		endif
	    goto 741
 740		erru=0.0d0
		do i=1,ndim
			drf(i)=dgres(i)
		enddo
		errp=0.0d0
		do i=1,ndim2
			ik=idu(i+nndtem)-nndtem
			ik=idpu(ik)
			epci(ik)=drf(nndtem+i)
			epc(ik,1)=epc(ik,1)+epci(ik)
			epcv(ik,1)=epcv(ik,1)+epci(ik)/beta/dt
             	errp=errp+drf(nndtem+i)**2
		enddo
		call rebound
		if (nfsi .eq. 3) then
			if (nimmer .eq. 1) then
				rnormps=0.0d0
				do i=1,nndaim*npart
					rnormps=rnormps+drf(nndtem1-i)**2
				enddo
				rnormpf=0.0d0
				do i=1,neqp-nndaim*npart-icount4
					rnormpf=rnormpf+drf(nndtem+i)**2
				enddo
				do i=1,icount4+nndaim*npart
					ntt=neqp-nndaim*npart-icount4+i
					rnormps=rnormps+drf(nndtem+ntt)**2
				enddo
			else
				rnormpf=0.0d0
				do i=1,neqp-icount4
					rnormpf=rnormpf+drf(nndtem+i)**2
				enddo
				rnormps=0.0d0
				do i=1,icount4
					rnormps=rnormps+drf(nndtem+neqp-i+1)**2
				enddo
			endif
		endif
		errp=dsqrt(errp)/pnorm
		erru=dsqrt(erru)/unorm
!		errp=dsqrt(errp)/errp0
!		erru=dsqrt(erru)/erru0
		if (nfsi .eq. 3) then
			rnormpf=dsqrt(rnormpf)/(nndfsia1-icount3)
			rnormps=dsqrt(rnormps)/nndfsia2
			rnormuf=dsqrt(rnormuf)/(nndfsi1-icount1)
			rnormus=dsqrt(rnormus)/nndfsi2
             ja=rnormps .ge. 1.0d-3
		else
			if (nnorm .eq. 1) then
! 				ja=erru .ge. tolu
				ja=errp .ge. tolp
			endif
		endif            
		if (ja) then
			if (ii .lt. miter) then
				ii=ii+1
				if (nfsi .eq. 3) then
					write(*,212) ii,rnormuf,rnormus,rnormpf,rnormps	
				else
					write(*,*) ii,erru,errp
				endif
				goto 1013
			endif
			if (ii .ge. miter) then
				do i=1,nnd
					write(22,*) i,dis(1,i,1),dis(2,i,1)
				enddo
				write(22,*)
				nonfl=0
				return
			endif
		endif
		ii=ii+1
 9998		if (nfsi .eq. 3) then    
			write(*,212) ii,rnormuf,rnormus,rnormpf,rnormps
		else
			write(*,*) ii,erru,errp
		endif
 212		format(1x,'uf us pf ps',i5,2x,4(e12.5,2x))
		xene=0.0d0
 100		format(1x,'Energy =',e23.16)
		if (nfrees .eq. 1) then
			call energy(xene,xvol)
		endif
		write(*,103) xene,xvol
 103		format(1x,'energy=',e23.16,1x,'vol=',e23.16)
		if (iflag .eq. 1) then
			do i=1,nndtem1
				drf(i)=0.0d0
			enddo
			call sstress
		endif
		if (nreact .eq. 1) then
			do ne=1,numeb
				nflag=nface(ne)
				if ((nale .eq. 1) .or. (nfsi .gt. 1)) then
					do noj=1,2
						do nos=1,nis(2)
							ntem=nea(nbe(ne),nos)
							x(noj,nos)=coor(ntem,noj)+dis(noj,ntem,1)
						enddo
					enddo
				else
					do noj=1,2
						do nos=1,nis(2)
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
			if (numfn .gt. 0) then
				call nodalf
			endif
			do i=1,nrtp
				nm=ndraf(i)+2*(nraf(i)-1)
				raf(nraf(i),ndraf(i))=-drf(nm)
			enddo
		endif
		if ((nstr .eq. 3) .and. (iflag .eq. 1)) then
			do ne=1,numel
				write(31,621)
 621				format(1x,'ele  gauss pts   stress strain-yy', &
     				&'stress strain-zz  stress strain-yz', &
     				&'stress strain-xx'/)
				do j=1,3
					do k=1,3
						write(31,602) ne,j,k,(cstr(m,ne,j,k),m=1,4)
						write(31,601) ne,j,k,(ge(m,ne,j,k),m=1,4)
 601						format(1x,i2,2x,2(i1,1x),4(e16.8,1x))
 602						format(1x,i2,2x,2(i1,1x),4(e16.8,1x))
						write(31,*)
					enddo
				enddo
				write(31,*)
			enddo
		endif
		if ((iti .eq. nts) .or. (iti .eq. 1)) then
			if (nisa .eq. 9) then
				call printnine
			else	
				call print
			endif 
		endif
	enddo
      if (nresave .eq. 1) then
		do i=1,nnd
			write(25,*) i,(dis(j,i,1),j=1,2),(vel(m,i,1),m=1,2),&
     			(velm(n,i,1),n=1,2),(acm(k,i,1),k=1,2)
		enddo
          do i=1,nndp
			write(25,*) i,epc(i,1),epcv(i,1)
		enddo
		write(25,*) nts
		close(25)
      endif
      close(22)
      return
      end