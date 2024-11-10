      subroutine assbf(wtt,wtt1,wtt2,nflag,nee)
      implicit real*8 (a-h,o-z)
      include 'common'
      nee1=id(2*(nee-1)+1)
      nee2=id(2*(nee-1)+2)
	if (nee1 .ne. 0) then
		if (nflag .eq. 1) then
			drf(nee1)=drf(nee1)-wtt1
		endif  
		if (nflag .eq. 3) then
			drf(nee1)=drf(nee1)+wtt1
		endif
		if (nflag .eq. 2) then
			drf(nee2)=drf(nee2)+wtt1
		endif
		if (nflag .eq. 4) then
			drf(nee2)=drf(nee2)+wtt1
		endif
	endif
	if (nee2 .ne. 0) then 
		if (nflag .eq. 1) then
			drf(nee1)=drf(nee1)-wtt2
		endif  
		if (nflag .eq. 3) then
			drf(nee1)=drf(nee1)+wtt2
		endif
		if (nflag .eq. 2) then
			drf(nee2)=drf(nee2)+wtt2
		endif
		if (nflag .eq. 4) then
			drf(nee2)=drf(nee2)+wtt2
		endif
	endif
      return
      end