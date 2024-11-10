      subroutine inputdat
      implicit real*8 (a-h,o-z)
      include 'common'
      logical jj
	nimmer=0
      read(1,*) nsolcon,isys
      read(1,*) ntest,ncon,nrestart,nresave
      read(1,*) nsmoo,smoo,nsmfun
      if ((nsmoo .eq. 4) .or. (nsmoos .eq. 4)) then
         read(1,*) xco
         read(1,*) niso
      endif
      read(1,*) nout,ninit,initdir
      read(1,*) npr,nprp,ntprint
	if (npr .ne. 0) then
		read(1,*) (nprint(i),i=1,npr)
	endif
	if (nprp .ne. 0) then
		read(1,*) (nprintp(i),i=1,nprp)
	endif
      read(1,*) nint
	read(1,*) nfsi
	if (nfsi .eq. 1) then
		ignum=1
		nndfsi2=0
		nndfsia2=0
		numele(2)=0
		read(1,*) nndfsi1
	    read(1,*) nump,nndfsia1
	    read(1,*) nis(1)
		read(1,*) numele(1)
		do j=1,numele(1)
			read(53,*) ntt,(nea(j,k),k=1,nis(1))
		enddo
		do j=1,numele(1)
			read(54,*) ntt,(neap(j,k),k=1,4),ntx1,ntx2,&
     		ntx3,ntx4,ntx5
	    enddo
		do j=1,nndfsi1
			read(51,*) ntxx1,xtxx1,coor(j,1),coor(j,2),ntxx2
		enddo
		do j=1,nndfsia1
			read(52,*) ntxx1,xtxx1,coorp(j,1),coorp(j,2),ntxx2
		enddo
		nnd=nndfsi1
		numel=numele(1)
 		nndp=nndfsia1
		numelep=0
		if (nis(1) .eq. 9) then
			do j=1,numele(1)
				nn1=j
				numtt=numelep+(j-1)*4+1
				numep(numtt,1)=nea(nn1,1)
				numep(numtt,2)=nea(nn1,5)
				numep(numtt,3)=nea(nn1,9)
				numep(numtt,4)=nea(nn1,8)
				numtt=numelep+(j-1)*4+2
				numep(numtt,1)=nea(nn1,5)
				numep(numtt,2)=nea(nn1,2)
				numep(numtt,3)=nea(nn1,6)
				numep(numtt,4)=nea(nn1,9)
				numtt=numelep+(j-1)*4+3
				numep(numtt,1)=nea(nn1,8)
				numep(numtt,2)=nea(nn1,9)
				numep(numtt,3)=nea(nn1,7)
				numep(numtt,4)=nea(nn1,4)
				numtt=numelep+(j-1)*4+4
				numep(numtt,1)=nea(nn1,9)
				numep(numtt,2)=nea(nn1,6)
				numep(numtt,3)=nea(nn1,3)
				numep(numtt,4)=nea(nn1,7)
				do k=1,4
					nmap(neap(nn1,k))=nea(nn1,k)
				enddo
			enddo
          	numelep=numelep+4*numele(1)
		endif
	endif
	if (nfsi .eq. 2) then
		ignum=2
		nndfsi1=0
		nndfsi1a=0
		numele(1)=0
		nis(1)=0
		read(1,*) nndfsi2
	    read(1,*) nump
		read(1,*) nndfsia2
	    read(1,*) nis(2)
		read(1,*) numele(2)
		do j=1,numele(2)
			read(57,*) ntt,(nea(j,k),k=1,nis(2))
          enddo
		do j=1,numele(2)
			read(58,*) ntt,(neap(j,k),k=1,4),ntx1,ntx2,&
     			ntx3,ntx4,ntx5
          enddo
		do j=1,nndfsi2
			read(55,*) ntxx1,xtxx1,coor(j,1),coor(j,2),ntxx2
		enddo
		do j=1,nndfsia2
			read(56,*) ntxx1,xtxx1,coorp(j,1),coorp(j,2),ntxx2
		enddo
		nnd=nndfsi2
		numel=numele(2)
 		nndp=nndfsia2
		numelep=0
		if (nis(2) .eq. 9) then
			do j=1,numele(2)
				nn1=j
				numtt=numelep+(j-1)*4+1
				numep(numtt,1)=nea(nn1,1)
				numep(numtt,2)=nea(nn1,5)
				numep(numtt,3)=nea(nn1,9)
				numep(numtt,4)=nea(nn1,8)
				numtt=numelep+(j-1)*4+2
				numep(numtt,1)=nea(nn1,5)
				numep(numtt,2)=nea(nn1,2)
				numep(numtt,3)=nea(nn1,6)
				numep(numtt,4)=nea(nn1,9)
				numtt=numelep+(j-1)*4+3
				numep(numtt,1)=nea(nn1,8)
				numep(numtt,2)=nea(nn1,9)
				numep(numtt,3)=nea(nn1,7)
				numep(numtt,4)=nea(nn1,4)
				numtt=numelep+(j-1)*4+4
				numep(numtt,1)=nea(nn1,9)
				numep(numtt,2)=nea(nn1,6)
				numep(numtt,3)=nea(nn1,3)
				numep(numtt,4)=nea(nn1,7)
				do k=1,4
					nmap(neap(nn1,k))=nea(nn1,k)
				enddo
			enddo
          	numelep=numelep+4*numele(2)
		endif
	endif
	if (nfsi .eq. 3) then
		ignum=2
		read(1,*) nimmer
		if (nimmer .eq. 1) then
			read(1,*) ndelta
			read(1,*) ndiv1,xmax,xmin
	        read(1,*) ndiv2,ymax,ymin
			read(1,*) nparticle
			read(1,*) npart
			read(1,*) pargrav(1),pargrav(2)
			if (nparticle .ne. 1) then
				do ip=1,npart
					read(1,*) rratio(ip),rcenter(ip,1),rcenter(ip,2)
				enddo	
			else
				do ip=1,npart
					read(1,*) rcenter(ip,1),rcenter(ip,2)
				enddo
			endif	
			read(1,*) numeleim,nndim,nndaim				
		endif
		read(1,*) nndfsi1,nndfsi2
		read(1,*) nump,nndfsia1
		read(1,*) nndfsia2
	    read(1,*) (nis(i),i=1,ignum)
		read(1,*) (numele(i),i=1,ignum)
	    do ig=1,ignum
			if (ig .eq. 1) then
				do j=1,numele(ig)
					read(53,*) ntt,(nea(j,k),k=1,nis(ig))
				enddo
				do j=1,numele(ig)
					read(54,*) ntt,(neap(j,k),k=1,4),ntx1,ntx2,&
     				ntx3,ntx4,ntx5
				enddo
			endif
			if (ig .eq. 2) then
				do j=1,numele(ig)
					read(57,*) ntt,(neatemp(j,k),k=1,nis(ig))
				enddo
				do j=1,numele(ig)
					read(58,*) ntt,(neaptemp(j,k),k=1,4),nt1,nt2,&
     				ntx3,ntx4,ntx5
				enddo
			endif
		enddo
		do j=1,nndfsi1
			read(51,*) ntxx1,xtxx1,coor(j,1),coor(j,2),ntxx2
		enddo
		do j=1,nndfsi2
			read(55,*) ntxx1,xtxx1,coortemp(j,1),coortemp(j,2),ntxx2
		enddo
		do j=1,nndfsia1
			read(52,*) ntx1,xtx1,coorp(j,1),coorp(j,2),ntxx2
		enddo
		do j=1,nndfsia2
			read(56,*) ntx1,xtx1,coorptemp(j,1),coorptemp(j,2),ntxx2
		enddo
		nnd=nndfsi1
		icount1=0
		icount2=0
		do j=1,nndfsi2
			itest=0
			do m=1,nndfsi1
				xtemp=dsqrt((coortemp(j,1)-coor(m,1))**2+&
     			(coortemp(j,2)-coor(m,2))**2)
				if (xtemp .le. 1.0d-14) then
					mapnea(j)=m
					icount1=icount1+1
					itest=1
					goto 337 
				endif
 			enddo		
			if (itest .eq. 0) then
				icount2=icount2+1
				nnd=nnd+1
				mapnea(j)=nnd
				coor(nnd,1)=coortemp(j,1)
				coor(nnd,2)=coortemp(j,2)
			endif		
 337		enddo			
		do j=1,numele(2)
			nn1=numele(1)+j
			do k=1,nis(2)
				nxtemp=neatemp(j,k)
				nea(nn1,k)=mapnea(nxtemp)
			enddo
		enddo
		write(*,*) icount1,nndfsi1+nndfsi2-nnd
		write(*,*) icount2,nndfsi2-icount1
 		nndp=nndfsia1
		icount3=0
		icount4=0
		do j=1,nndfsia2
			itest=0
			do m=1,nndfsia1
				xtemp=dsqrt((coorptemp(j,1)-coorp(m,1))**2+&
     				(coorptemp(j,2)-coorp(m,2))**2)
				if (xtemp .le. 1.0d-14) then
					icount3=icount3+1
					mapneap(j)=m
					itest=1
					goto 336 
				endif
 			enddo		
			if (itest .eq. 0) then
				icount4=icount4+1
				nndp=nndp+1
				mapneap(j)=nndp
				coorp(nndp,1)=coorptemp(j,1)
				coorp(nndp,2)=coorptemp(j,2)
			endif			
 336		enddo			
		do j=1,numele(2)
			nn1=numele(1)+j
			do k=1,4
				nxtemp=neaptemp(j,k)
				neap(nn1,k)=mapneap(nxtemp)
			enddo
		enddo
		write(*,*) icount3,nndfsia1+nndfsia2-nndp
		write(*,*) icount4,nndfsia2-icount3
		numelt=numele(1)+numele(2)
		if (nimmer .eq. 1) then
			do j=1,numeleim
				nn1=numelt+j
				read(72,*) ntt,(neatemp(nn1,k),k=1,nis(2))
				do k=1,nis(2)
					nea(nn1,k)=neatemp(nn1,k)+nnd
				enddo
			enddo
			do ip=1,npart-1
				do i=1,numeleim
					do k=1,nis(2)
						ntt1=numelt+ip*numeleim+i
						ntt2=numelt+i
						nshi=nndim*ip+nnd	
						nea(ntt1,k)=neatemp(ntt2,k)+nshi
					enddo
				enddo
			enddo
			do j=1,numeleim
				nn1=numelt+j
				read(73,*) ntt,(neaptemp(nn1,k),k=1,4),ntx1,ntx2,&
     					ntx3,ntx4,ntx5
				do k=1,4
					neap(nn1,k)=neaptemp(nn1,k)+nndp
				enddo
			enddo
			do ip=1,npart-1
				do i=1,numeleim
					do k=1,4
						ntt1=numelt+ip*numeleim+i
						ntt2=numelt+i	
						nshi=nndaim*ip+nndp		
						neap(ntt1,k)=neaptemp(ntt2,k)+nshi
					enddo
				enddo
			enddo
			do k=1,nndim
				j=nndfsi2+k
				read(71,*) ntx1,xtx1,coortemp(j,1),coortemp(j,2),ntx2
			enddo
			do k=1,nndaim
				j=nndfsia2+k
				read(74,*) nx1,x1,coorptemp(j,1),coorptemp(j,2),nx2
			enddo				
			if (nparticle .ne. 1) then
				do ip=1,npart
					if (ip .eq. 1) then
					do j=1,nndim
						ntt1=j+nndfsi2
					coortempo(ntt1,1)=coortemp(ntt1,1)
					coortempo(ntt1,2)=coortemp(ntt1,2)
					enddo
					do j=1,nndaim
						nt1=j+nndfsia2
					coorptempo(nt1,1)=coorptemp(nt1,1)
					coorptempo(nt1,2)=coorptemp(nt1,2)
					enddo
					endif
					do j=1,nndim
						ntt1=j+(ip-1)*nndim+nndfsi2
						ntt2=j+nndfsi2
					coortemp(ntt1,1)=coortempo(ntt2,1)/rratio(ip)
					coortemp(ntt1,1)=coortemp(ntt1,1)+rcenter(ip,1)
					coortemp(ntt1,2)=coortempo(ntt2,2)/rratio(ip)
					coortemp(ntt1,2)=coortemp(ntt1,2)+rcenter(ip,2)
					enddo
					do j=1,nndaim
						nt1=j+(ip-1)*nndaim+nndfsia2
						nt2=j+nndfsia2
					coorptemp(nt1,1)=coorptempo(nt2,1)/rratio(ip)
					coorptemp(nt1,1)=coorptemp(nt1,1)+rcenter(ip,1)
					coorptemp(nt1,2)=coorptempo(nt2,2)/rratio(ip)
					coorptemp(nt1,2)=coorptemp(nt1,2)+rcenter(ip,2)
					enddo
				enddo
			else
				do ip=1,npart
					do j=1,nndim
						ntt1=j+(ip-1)*nndim+nndfsi2
						ntt2=j+nndfsi2
					coortemp(ntt1,1)=coortemp(ntt1,1)+rcenter(ip,1)
					coortemp(ntt1,2)=coortemp(ntt1,2)+rcenter(ip,2)
					enddo
				enddo				
			endif
			do ip=1,npart
			do i=1,nndim
				ni=i+nnd+(ip-1)*nndim
				nj=i+(ip-1)*nndim+nndfsi2
				coor(ni,1)=coortemp(nj,1)
				coor(ni,2)=coortemp(nj,2)
			enddo
			do i=1,nndaim
				ni=i+nndp+(ip-1)*nndim
				nj=i+(ip-1)*nndaim+nndfsia2
				coorp(ni,1)=coorptemp(nj,1)
				coorp(ni,2)=coorptemp(nj,2)
			enddo
			enddo
			do j=1,npart*nndim
				nj=nndfsi2+j
				mapnea(nj)=nnd+j
 			enddo		
			do j=1,npart*nndaim
				nj=nndfsia2+j
				mapneap(nj)=nndp+j
 			enddo
			nnd=nnd+nndim*npart
			nndp=nndp+nndaim*npart
			nndfsi2=nndfsi2+nndim*npart
			nndfsia2=nndfsia2+nndaim*npart
			numele(2)=numele(2)+numeleim*npart
      		numel=numele(1)+numele(2)							
	    endif
		numelep=0
		do ig=1,ignum
			nel11=0
			do i=1,ig-1
				nel11=nel11+numele(i)
			enddo
			if (nis(ig) .eq. 9) then
				do j=1,numele(ig)
						nn1=nel11+j
						numtt=numelep+(j-1)*4+1
						numep(numtt,1)=nea(nn1,1)
						numep(numtt,2)=nea(nn1,5)
						numep(numtt,3)=nea(nn1,9)
						numep(numtt,4)=nea(nn1,8)
						numtt=numelep+(j-1)*4+2
						numep(numtt,1)=nea(nn1,5)
						numep(numtt,2)=nea(nn1,2)
						numep(numtt,3)=nea(nn1,6)
						numep(numtt,4)=nea(nn1,9)
						numtt=numelep+(j-1)*4+3
						numep(numtt,1)=nea(nn1,8)
						numep(numtt,2)=nea(nn1,9)
						numep(numtt,3)=nea(nn1,7)
						numep(numtt,4)=nea(nn1,4)
						numtt=numelep+(j-1)*4+4
						numep(numtt,1)=nea(nn1,9)
						numep(numtt,2)=nea(nn1,6)
						numep(numtt,3)=nea(nn1,3)
						numep(numtt,4)=nea(nn1,7)
						do k=1,4
							nmap(neap(nn1,k))=nea(nn1,k)
							coorp(neap(nn1,k),1)=coor(nea(nn1,k),1)
							coorp(neap(nn1,k),2)=coor(nea(nn1,k),2)
						enddo
				enddo
          		numelep=numelep+4*numele(ig)
			endif
		enddo
	endif
      do i=1,2*nnd
         id(i)=1
      enddo
	do i=1,nndp
	   idp(i)=1
	enddo
      read(1,*) nfrees
      if (nfrees .eq. 1) then
         read(1,*) ncons,nsmoos,nsur
         read(1,*) neles,niss,gravit
         do i=1,neles
            read(1,*) (nfsnd(i,j),j=1,niss)
	   enddo
      endif         
      read(1,*) numskew,numgb,intnum,numct
      do j=1,numgb
         read(1,*) ndirgb(j)
         read(1,*) numdir(j)
		do k=1,numdir(j)
		read(1,*) nodegb(j,k)
		enddo
	enddo 
	if ((nfsi .eq. 3) .or. (nfsi .eq. 1)) then
		ntemb=numgb
		icount0=0
		do j=1,ntemb
			if (ndirgb(j) .eq. 333333) then
				do i=1,numdir(j)
					ntema=nodegb(j,i)
					itest=0
					do k=1,numdir(1)
						if (mapnea(ntema) .eq. nodegb(1,k)) then
							itest=1
							icount0=icount0+1
						endif
					enddo	
					if (itest .eq. 0) then
						numdir(1)=numdir(1)+1
						nodegb(1,numdir(1))=mapnea(ntema)
					endif
				enddo
				numgb=numgb-1
			endif
			if (ndirgb(j) .eq. 222222) then
				do i=1,numdir(j)
					idp(nodegb(j,i))=0
				enddo
				numgb=numgb-1
			endif
		enddo
!		do j=1,npr
!			ntemm=nprint(j)
!	write(*,*) npr,j,nprint(j),mapnea(ntemm)
!			nprint(j)=mapnea(ntemm)
!		enddo
	endif
	do i=1,intnum
         read(1,*) numint(i),ninsk(i)
      enddo
      do i=1,numskew
         read(1,*) xang(i)
      enddo
      do i=1,numct
         read(1,*) nodesl(i),nodema(i),islavdir(i),&
            imasdir(i),amct(i)
      enddo
      read(1,*) fbulk,falph,fdensi
      if (nfsi .ne. 1) then
		read(1,*) nstr
		if (nstr .eq. 1) then
			read(1,*) syoun,sdensi,spois
			read(1,*) sarea,xl
			read(1,*) nodecon,npss 
		endif
		if (nstr .eq. 2) then
			read(1,*) syoun,sdensi,spois
			read(1,*) npss 
		endif
		if (nstr .eq. 3) then
			read(1,*) rc1,rc2,rk,sdensi
			read(1,*) iflag,nales
		endif
		if (nstr .eq. 4) then
			do ni=1,npart*nndim
				read(1,*) xmass(ni),xstif(ni)
			enddo	
		endif
      endif
	read(1,*) numtfun
	read(1,*) (ntfun(i),i=1,numtfun)
	do i=1,numtfun
		if (ntfun(i) .eq. 2) then
			read(1,*) xome(i)
		endif
	enddo
	read(1,*) dt,nts
	read(1,*) nina
      read(1,*) alpha,beta
	read(1,*) numfn,numeb     
      do i=1,numeb
		read(1,*) nbe(i),nface(i),boup(i,1),boup(i,2),ntb(i),nbfsi(i)
	enddo
	if (nfsi .eq. 3) then
		do i=1,numeb
			if (nbfsi(i) .eq. 2) then
				nbe(i)=nbe(i)+numele(1)
			endif
		enddo
	endif
      do i=1,numfn
		read(1,*) nodefn(i),ndirfn(i),ftem
		fnod(nodefn(i),ndirfn(i))=ftem
	enddo
      read(1,*) fbacc(1),fbacc(2)
      do i=1,2
		read(1,*) nfuns(i)
		if (nfuns(i) .eq. 2) then
			read(1,*) xome(i)
          endif
	enddo
      read(1,*) nale,nnorm
      read(1,*) miter
      read(1,*) tolu,tolp,unorm,pnorm
      read(1,*) thic,nreact
      if (nreact .eq. 1) then
		read(1,*) nrtp
          do i=1,nrtp
			read(1,*) ntt,mtt
			ndraf(i)=ntt
			nraf(i)=mtt
          enddo
      endif
      if (nale .eq. 1) then
         do i=1,nndfsi1
            n1link(i)=1
            n2link(i)=1
            nod1link(1,i)=i
            nod2link(1,i)=i
            coe1link(1,i)=0.0d0
            coe2link(1,i)=0.0d0
          enddo
          i=0
 1000    read(59,*) ndum
         i=i+1
         if (ndum .ne. 0) then
            read(59,*) n1link(ndum),n2link(ndum)
            if (n1link(ndum) .ge. 1) then
               do j=1,n1link(ndum)
                  read(59,*) nod1link(j,ndum),coe1link(j,ndum)
               enddo
            endif
            if (n2link(ndum) .ge. 1) then
               do k=1,n2link(ndum)
                  read(59,*) nod2link(k,ndum),coe2link(k,ndum)
               enddo
            endif
            goto 1000
         endif
      endif
	neqt=0
      do i=1,numgb
         if (ndirgb(i) .eq. 111111) then
             do j=1,numdir(i)
               id(2*(nodegb(i,j)-1)+1)=0
			 neqt=neqt+1
	         nodeb(neqt)=2*(nodegb(i,j)-1)+1
	         id(2*(nodegb(i,j)-1)+2)=0
			 neqt=neqt+1
			 nodeb(neqt)=2*(nodegb(i,j)-1)+2
		   enddo
	   endif
         if (ndirgb(i) .eq. 110111) then
             do j=1,numdir(i)
               id(2*(nodegb(i,j)-1)+1)=0
			 neqt=neqt+1
      		 nodeb(neqt)=2*(nodegb(i,j)-1)+1
		   enddo
	   endif
         if (ndirgb(i) .eq. 101111) then
             do j=1,numdir(i)
	         id(2*(nodegb(i,j)-1)+2)=0
			 neqt=neqt+1
			 nodeb(neqt)=2*(nodegb(i,j)-1)+2
		   enddo
	   endif
         if (ndirgb(i) .eq. 121111) then
             do j=1,numdir(i)
               id(2*(nodegb(i,j)-1)+1)=0
			 neqt=neqt+1
      		 nodeb(neqt)=2*(nodegb(i,j)-1)+1
		   enddo
	   endif
         if (ndirgb(i) .eq. 112111) then
             do j=1,numdir(i)
	         id(2*(nodegb(i,j)-1)+2)=0
			 neqt=neqt+1
			 nodeb(neqt)=2*(nodegb(i,j)-1)+2
		   enddo
	   endif
      enddo
      return
      end