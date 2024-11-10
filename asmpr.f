      subroutine asmpr(ne)
      implicit real*8 (a-h,o-z)
      include 'common'
      do i=1,nump
		ni=nndtem+idp(neap(ne,i))
          drf(ni)=drf(ni)-xfp(i,ne)
          sk(ni)=sk(ni)+xkpp(i,ne)
      enddo
      return
      end