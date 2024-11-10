      subroutine surcon
      implicit real*8 (a-h,o-z)
      include 'common'
      if (niss .eq. 3) then
         nnds=2*neles+1
      endif
      if (niss .eq. 2) then
         nnds=neles+1
      endif
      do i=1,neles
         if (niss .eq. 3) then
            neas(i,1)=2*(i-1)+1
            neas(i,2)=2*i+1
            neas(i,3)=2*i
         endif
         if (niss .eq. 2) then
            neas(i,1)=i
            neas(i,2)=i+1
         endif
         do j=1,niss
            ntem=nfsnd(i,j)
            mtem=neas(i,j)
            nnn(mtem)=ntem
            nnni(ntem)=mtem
          enddo
      enddo
      return
      end