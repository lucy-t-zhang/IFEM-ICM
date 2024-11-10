      subroutine sur_1d(det,n1d,n1c)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension n1d(2),n1c(2)
      nu1=n1d(1)
      nu2=n1d(2)
      nc1=n1c(1)
      nc2=n1c(2)
      vip1=0.5d0*(vel(1,nc1,1)+vel(1,nc2,1))
      dip1=0.5d0*(dis(2,nc1,1)+dis(2,nc2,1))
      ddip1=0.5d0*(dis(2,nc1,1)-dis(2,nc2,1))/det
      ddip1=dsqrt(1.0d0+ddip1**2)
      if (nsmoo .eq. 5) then
         if (vip1 .gt. 0.0d0) then
            dip1=dis(2,nc2,1)
            sk(nu2)=sk(nu2)+alpha*dt/beta*vip1*ddip1
         endif
         if (vip1 .lt. 0.0d0) then
            dip1=dis(2,nc1,1)
            sk(nu1)=sk(nu1)-alpha*dt/beta*vip1*ddip1
         endif
      else
         sk(nu2)=sk(nu2)+0.5d0*alpha*dt/beta*vip1*ddip1
         sk(nu1)=sk(nu1)-0.5d0*alpha*dt/beta*vip1*ddip1
      endif
      drf(nu2)=drf(nu2)-vip1*dip1*ddip1
      drf(nu1)=drf(nu1)+vip1*dip1*ddip1
      return
      end