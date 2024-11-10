      subroutine eleinit(ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      do i=1,nump
		xkpp(i,ne)=0.0d0
      enddo
      if (nale .eq. 1) then
         if (nfrees .eq. 1) then
            ntem=nnds+2*nnd
         else
            ntem=2*nnd
         endif
      endif
      do i=1,nis(ig)
         xfu(i)=0.0d0
         xfu(i+nis(ig))=0.0d0
      enddo
      do i=1,nump
         xfp(i,ne)=0.0d0
      enddo
      return
      end