      subroutine printnine
      implicit real*8 (a-h,o-z)
      include 'common'
      write(3,105) iti+1,iti*dt
 105  format(1x, 'time(',i4,')=',e23.7,';'/)
      do ip=1,npr
		xtt1a=dis(1,nprint(ip),1)+&
     		coor(nprint(ip),1)
		xtt2a=dis(2,nprint(ip),1)+&
     		coor(nprint(ip),2)
		xtt1=dis(1,nprint(ip),1)
          xtt2=dis(2,nprint(ip),1) 
          write(3,124) ip,iti+1,xtt1a
 124      format(1x, 'coorx(',i4,',',i4,')=',e23.7,';'/)
          write(3,125) ip,iti+1,xtt2a
 125      format(1x, 'coory(',i4,',',i4,')=',e23.7,';'/)
          write(3,104) ip,iti+1,xtt1
 104      format(1x, 'disx(',i4,',',i4,')=',e23.7,';'/)
          write(3,135) ip,iti+1,xtt2
 135      format(1x, 'disy(',i4,',',i4,')=',e23.7,';'/)
		xtt1v=vel(1,nprint(ip),1)
		xtt2v=vel(2,nprint(ip),1)
		write(3,114) ip,iti+1,xtt1v
 114		format(1x, 'velx(',i5,',',i4,')=',e23.7,';'/)
		write(3,115) ip,iti+1,xtt2v
 115		format(1x, 'vely(',i5,',',i4,')=',e23.7,';'/)
	enddo
	do ip=1,nprp
		xttp=epc(nprintp(ip),1)
		write(3,116) ip,iti+1,xttp
 116		format(1x, 'xpres(',i5,',',i4,')=',e23.7,';'/)
		xttx=coorp(nprintp(ip),1)
		write(3,117) ip,iti+1,xttx
 117		format(1x, 'xpcoorx(',i5,',',i4,')=',e23.7,';'/)
		xtty=coorp(nprintp(ip),2)
		write(3,118) ip,iti+1,xtty
 118		format(1x, 'xpcoory(',i5,',',i4,')=',e23.7,';'/)
	enddo
      na=iti/nina-nai
      if (na .ne. 0) then
		write(22,*) na,iti+1,nts,nina
		write(23,*) na,iti+1,nts,nina
		do i=1,nndp
			write(42,*) i,dis(1,nmap(i),1),dis(2,nmap(i),1)
			write(43,*) i,vel(1,nmap(i),1),vel(2,nmap(i),1)
		enddo
		do i=1,nnd
			write(22,*) i,dis(1,i,1),dis(2,i,1)
			write(23,*) i,vel(1,i,1),vel(2,i,1)
	    enddo
		write(24,*) na,iti+1,nts,nina
		write(84,*) na,iti+1,nts,nina
		write(94,*) na,iti+1,nts,nina
		do ni=1,nndp-nndaim*npart
			write(24,*) ni,epc(ni,1)
			write(84,*) ni,epcv(ni,1)
			write(94,*) ni,epc(ni,1)
          enddo
		do ni=1,nndaim*npart
			nj=nndp-nndaim*npart+ni
			write(24,*) nndp-nndaim*npart+ni,epco(ni,1)
			write(84,*) nndp-nndaim*npart+ni,epcov(ni,1)
			write(94,*) nj,epc(nj,1)
		enddo
		if (nreact .eq. 1) then
			do k=1,nrtp
			write(33,607) ndraf(k),nraf(k),iti+1,raf(nraf(k),ndraf(k))
 607			format(1x, 'rforce',i1,i3,'(',i4,')=',e23.7,';'/)
			enddo
		endif
		write(22,*)
		nai=nai+1     
	endif
      return
      end