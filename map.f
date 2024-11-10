      subroutine map
      implicit real*8 (a-h,o-z)
      include 'common'
      logical ja,jb
      neqv=0
      do n=1,2*nnd
         if (id(n) .eq. 1) then
            neqv=neqv+1
            id(n)=neqv
            idu(neqv)=n
         endif
	enddo
	neqp=0
      do n=1,nndp
         if (idp(n) .eq. 1) then
            neqp=neqp+1
            idp(n)=neqp
            idpu(neqp)=n
         endif
	enddo
      return
      end