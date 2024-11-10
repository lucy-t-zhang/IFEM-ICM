      subroutine init
      implicit real*8 (a-h,o-z)
      include 'common'
      if (nale .eq. 1) then
		do i=1,nndfsi1
			do j=1,2
				velm(j,i,1)=0.0d0
				acmm(j,i,1)=0.0d0
			enddo
		enddo
	endif
	do i=1,nnd
		do j=1,2
			dis(j,i,1)=0.0d0
              acm(j,i,1)=0.0d0
			vel(j,i,1)=0.0d0
		enddo
	enddo
      do ni=1,nndp
		epc(ni,1)=0.0d0
		epcv(ni,1)=0.0d0
      enddo
	do i=1,nndaim*npart
		epco(i,1)=0.0d0
		epcov(i,1)=0.0d0
	enddo
      if (ninit .eq. 1) then
         call readinit
      endif
      return
      end	