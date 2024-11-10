      subroutine nodalf
      implicit real*8 (a-h,o-z)
      include 'common'
	do i=1,numfn
		ni=ndirfn(i)+2*(nodefn(i)-1)
		drf(id(ni))=drf(id(ni))+fnodo(nodefn(i),ndirfn(i))
	enddo
      return
      end