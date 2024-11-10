      subroutine tridag(r,u,n)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension r(n),u(n)
      if (sk(1) .eq. 0.0d0) pause 'tridag: rewrite equations'
      bet=sk(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        sgam(j)=sc(j-1)/bet
        bet=sk(j)-sa(j)*sgam(j)
        if (bet .eq. 0.0d0) pause 'tridag failed'
        u(j)=(r(j)-sa(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-sgam(j+1)*u(j+1)
12    continue
      return
      end