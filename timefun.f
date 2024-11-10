      subroutine timefun
      implicit real*8 (a-h,o-z)
      include 'common'
	do i=1,numtfun
		nttfun=ntfun(i)
		if (nttfun .eq. 1) then
			tmid=10.0d0
			ntmid=dnint(tmid/dt)
			if (iti .lt. ntmid) then
				tfun(i)=dt*iti/tmid*1.0d0
			else
				tfun(i)=1.0d0
			endif
		endif
		if (nttfun .eq. 2) then
			tfun(i)=0.01d0*dsin(2*pi*xome(i)*iti*dt)
		endif
		if (nttfun .eq. 3) then
			tfun(i)=1.0d0*iti/nts
		endif
		if (nttfun .eq. 4) then
			tfun(i)=1.0d0
		endif
		if (nttfun .eq. 5) then
			if (iti*dt .lt. 35.0d0) then
				tfun(i)=0.93d0-(0.93d0+2.325d0)*iti*dt/35.0d0
			else
				tfun(i)=1.0d0
			endif
			tfun(i)=3.0d-1*tfun(i)
		endif
		if (nttfun .eq. 6) then
			tmid1=0.1d0
			ntmid1=dnint(tmid1/dt)
			tmid2=0.4d0
	        ntmid2=dnint(tmid2/dt)
			if (iti .lt. ntmid1) then
				tfun(i)=2.0d0*dt*iti/tmid1
			else
				if (iti .lt. ntmid2) then
					tfun(i)=2.0d0-(1.0d0+2.0d0)*iti*dt/tmid2
				else
					tfun(i)=-1.0d0
				endif
			endif
			tfun(i)=2.0d1*tfun(i)
		endif
		if (nttfun .eq. 7) then
			tmid1=0.2d0
			ntmid1=dnint(tmid1/dt)
			tmid2=0.4d0
	        ntmid2=dnint(tmid2/dt)
			if (iti .lt. ntmid1) then
				tfun(i)=0.0d0
			else
				if (iti .lt. ntmid2) then
					tfun(i)=(iti-ntmid1)*dt/(tmid2-tmid1)
				else
					tfun(i)=1.0d0
				endif
			endif
		endif
	enddo
      return
      end