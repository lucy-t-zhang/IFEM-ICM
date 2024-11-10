      subroutine load
      implicit real*8 (a-h,o-z)
      logical ja
      include 'common'
      ja=((dabs(fbacc(1)) .gt. 0.0d0) .or. 
     $     (dabs(fbacc(2)) .gt. 0.0d0))      
      if (ntprint .eq. 1) then
		do i=1,numtfun
			write(4,100) iti+1,i,tfun(i)
 100			format(1x, 'tfun',i1,'(',i4,')=',e23.7,';'/)
		enddo
		write(4,101) iti+1,iti*dt
 101		format(1x, 'time(',i4,')=',e23.7,';'/)
      endif
      do k=1,numeb
		boupo(k,1)=boup(k,1)*tfun(ntb(k))
		boupo(k,2)=boup(k,2)*tfun(ntb(k))
      enddo
      do i=1,numfn
		fnodo(nodefn(i),ndirfn(i))=fnod(nodefn(i),ndirfn(i))*
     $		tfun(1)
      enddo
	pargravo(1)=pargrav(1)*tfun(2)
	pargravo(2)=pargrav(2)*tfun(2)
	if (ja) then
		do i=1,2
               fbacco(i)=fbacc(i)*tfun(1)
		enddo
      endif
      return
      end