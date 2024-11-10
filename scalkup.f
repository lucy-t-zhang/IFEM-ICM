c     
c     caculation of kup
c     
      subroutine scalkup(xkup,ocup,i,k,ni)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension ocup(6)
      xkup=0.0d0
      do 10 m=1,3
         xkup=xkup+ocup(m)*dge(m,i,ni)*hp(k)
   10 continue
      return
      end


