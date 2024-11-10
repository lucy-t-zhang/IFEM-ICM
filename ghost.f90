      subroutine ghost      
      implicit real*8 (a-h,o-z)
      include 'common'      
	dimension ntempa(2),ntemp(2)
	xinc=(xmax-xmin)/ndiv1
	yinc=(ymax-ymin)/ndiv2 
	xrf=dsqrt(xinc**2+yinc**2)  
	do i=1,ndiv1
		do j=1,ndiv2
			nmapp(i,j)=0
	        nmappa(i,j)=0
		enddo
	enddo
	do i=1,nnd-nndim*npart
		ntemp(1)=1+dint(coor(i,1)/xinc)
		ntemp(2)=1+dint(coor(i,2)/yinc)
		nmapp(ntemp(1),ntemp(2))=nmapp(ntemp(1),ntemp(2))+1
		npoints(ntemp(1),ntemp(2),nmapp(ntemp(1),ntemp(2)))=i
	enddo
	do i=1,nndp-nndaim*npart
		ntempa(1)=1+dint(coorp(i,1)/xinc)
		ntempa(2)=1+dint(coorp(i,2)/yinc)
		nmappa(ntempa(1),ntempa(2))=nmappa(ntempa(1),ntempa(2))+1
		npointsa(ntempa(1),ntempa(2),nmappa(ntempa(1),ntempa(2)))=i
	enddo
      return
      end 