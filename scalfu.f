      subroutine scalfu(fu,i,ni)
      implicit real*8 (a-h,o-z)
      include 'common'
      fu=0.0d0
      do m=1,3
         fu=fu+tos(m)*dge(m,i,ni)
      enddo
      return
      end