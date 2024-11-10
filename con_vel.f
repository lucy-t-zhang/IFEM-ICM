      subroutine con_vel(nloc)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension nloc(4)
      do 24 i=1,2
         fvel(i)=0.0d0
         do 25 k=1,4
            ntem=nloc(k)
            fvel(i)=fvel(i)+
     $           h(k)*(vel(i,ntem,1)-velm(i,ntem,1))
 25      continue
 24   continue
      return
      end