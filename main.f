c-------------------------------------*
c     Fluid-Structure Interaction     *
c     Navier-Stokes for fluids        *
c     Control Volume Upwinding        *
c     Rubber like material   ALE      *
c     Mixed Elements Inf-Sup          *
c     9/4c                            *
c     by X. Sheldon Wang              *
c-------------------------------------*
      program newfsi
      implicit real*8 (a-h,o-z)
      include 'common'
      integer time
	external time
      data xg/ 0.d0, 0.d0, 0.d0, 0.d0, -0.577350269189626d0,
     1     0.577350269189626d0, 0.d0, 0.d0, -0.774596669241483d0,0.d0, 
     2     0.774596669241483d0,0.d0, -0.861136311594053d0,
     3     -0.339981043584856d0, 0.339981043584856d0, 
     $     0.861136311594053d0/
      data wgt/ 2.d0, 0.d0,0.d0,0.d0,1.d0,1.d0,
     1     0.d0,0.d0,0.5555555555555556d0,0.8888888888888889d0,
     2     0.5555555555555556d0,0.d0, 0.347854845137454d0,
     $     0.652145154862546d0,
     3     0.652145154862546d0,0.347854845137454d0/
      open(1,file='coortable',status='old')
      open(71,file='imcoor.txt',status='old')
      open(72,file='imcont.txt',status='old')
      open(73,file='imcont4.txt',status='old')
      open(74,file='imcoor4.txt',status='old')
      open(61,file='initv.txt',status='old')
      open(51,file='finfcoor.txt',status='old')
      open(52,file='finfcoor4.txt',status='old')
      open(53,file='finfcont.txt',status='old')
      open(54,file='finfcont4.txt',status='old')
      open(55,file='finscoor.txt',status='old')
      open(56,file='finscoor4.txt',status='old')
      open(57,file='finscont.txt',status='old')
      open(58,file='finscont4.txt',status='old')
      open(59,file='outputmesh.txt',status='old')
      open(2,file='output',status='unknown')
      open(3,file='plotp.m',status='unknown')
      open(4,file='time.m',status='unknown')
      open(25,file='restart',status='unknown')
      open(19,file='node',status='unknown')
      open(20,file='coor',status='unknown')
      open(40,file='coor4',status='unknown')
      open(21,file='cont',status='unknown')
      open(41,file='cont4',status='unknown')
      open(22,file='disp',status='unknown')
      open(23,file='vel',status='unknown')
      open(42,file='disp4',status='unknown')
      open(43,file='vel4',status='unknown')
      open(24,file='pres',status='unknown')
      open(84,file='presv',status='unknown')
      open(94,file='presc',status='unknown')
      open(31,file='ssee',status='unknown')
      open(33,file='react',status='unknown')
      write(*,100)
	npresb=0
	naxx1=time()
  100 format(6x,'Welcome to Fluid-Structure Program')
      two1 = -2.0d0
      one1 = -1.0d0
      zero =  0.0d0
      one2 =  1.0d0
      two2 =  2.0d0
      tt(1,3)=0.0d0
      tt(2,3)=0.0d0
      tt(3,2)=0.0d0
      tt(3,1)=0.0d0
      tt(3,3)=1.0d0
      ss(1,3)=0.0d0
      ss(3,1)=0.0d0
      ss(2,3)=0.0d0
      ss(3,2)=0.0d0
	call con_cont
	call con_set
      pi=datan(1.0d0)*4.0d0
      call inputdat
	if (nimmer .eq. 1) then
		call ghost
	endif
	call map
      if (iflag .eq. 1) then
         xxrs(1,1)=1.0d0
         xxrs(1,2)=1.0d0
         xxrs(2,1)=-1.0d0
         xxrs(2,2)=1.0d0
         xxrs(3,1)=-1.0d0
         xxrs(3,2)=-1.0d0
         xxrs(4,1)=1.0d0
         xxrs(4,2)=-1.0d0
      endif
      if (nfrees .eq. 1) then
         xmat=fdensi*gravit
         call surcon
      endif
      if ((nfrees .eq. 1) .and. (nale .eq. 1)) then
         nndtem=neqv+nnds
      else
         nndtem=neqv
      endif
	if (nimmer .eq. 1) then
		nndtem=nndtem-nndim*npart*2
	endif
      nndtem1=nndtem+neqp
	ndim=nndtem1
	ndim2=neqp
	ndim1=nndtem
	do i=nndtem+1,nndtem1
		idu(i)=i
	enddo
      do i=1,nnd
         write(20,*) i,coor(i,1),coor(i,2)
      enddo
	if (nfsi .eq. 1) then
		nisa=nis(1)
	endif
	if (nfsi .eq. 2) then
		nisa=nis(2)
	endif
	if (nfsi .eq. 3) then
		nisa=nis(1)
	endif
	if (nisa .eq. 9) then
		do i=1,nndp
		   write(40,*) i,coor(nmap(i),1),coor(nmap(i),2)
		enddo	
	endif
	if (nisa .eq. 9) then
	    do i=1,numelep
		   write(21,*) i,(numep(i,j),j=1,4)
		enddo
	else
		do i=1,numel
			write(21,*) i,(nea(i,j),j=1,4)
		enddo
	endif
	if (nisa .eq. 9) then
		do i=1,numel
			write(41,*) i,(neap(i,j),j=1,4)
		enddo
	endif
      close(20)
      close(21)
	if (nrestart .ne. 1) then
		call init
	else
		if (nale .eq. 1) then
			do i=1,nnd
				read(25,*) ntti,(dis(j,i,1),j=1,2),
     $				(vel(m,i,1),m=1,2),
     $				(velm(n,i,1),n=1,2),(acm(k,i,1),k=1,2)
			enddo
		else           
			do i=1,nnd
				read(25,*) ntti,(dis(j,i,1),j=1,2),
     $				(vel(m,i,1),m=1,2),
     $				(acm(k,i,1),k=1,2)
			enddo
		endif
		do i=1,nndp
			read(25,*) ntti,epc(i,1),epcv(i,1)
		enddo
	endif
      do i=1,nnd
         write(16,*) i,n2link(i),nod2link(1,i),coe2link(1,i)
      enddo
      close(1) 
      nonfl=1
      if (nstr .eq. 3) then
		x13=1.0d0/3.0d0
		x23=2.0d0/3.0d0
		x43=4.0d0/3.0d0
		x53=5.0d0/3.0d0
		x73=7.0d0/3.0d0
		x83=8.0d0/3.0d0
		x49=4.0d0/9.0d0
		x109=10.0d0/9.0d0
      endif
      call nonass(nonfl)
      if (nonfl .eq. 0) then
		goto 1021
      endif
      goto 1022
 1021 write(*,225)
 225  format(6x,'failure in iteration')
      goto 9999
 1022 write(*,226)
 226  format(1x, 'calculation finished normally')
 9999 close(23)
      close(24)
      naxx2=time()
      write(*,*) naxx2-naxx1
      stop
      end