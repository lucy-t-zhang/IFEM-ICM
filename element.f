      subroutine element(rs)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension rs(2)
      r=rs(1)
      s=rs(2)
      h(9)=0.0d0
      if (nis(ig) .eq. 9) then 
         h(9)=(1.0d0-r**2)*(1.0d0-s**2)
      endif
      if (nis(ig) .ne. 4) then
         h(5)=0.5d0*(1.0d0-r**2)*(1.0d0+s)-0.5d0*h(9)
         h(6)=0.5d0*(1.0d0-s**2)*(1.0d0-r)-0.5d0*h(9)
         h(7)=0.5d0*(1.0d0-r**2)*(1.0d0-s)-0.5d0*h(9)
         h(8)=0.5d0*(1.0d0-s**2)*(1.0d0+r)-0.5d0*h(9)
         h(1)=0.25d0*(1.0d0+r)*(1.0d0+s)-0.5d0*h(5)-
     $        0.5d0*h(8)-0.25d0*h(9)
         h(2)=0.25d0*(1.0d0-r)*(1.0d0+s)-0.5d0*h(5)-
     $        0.5d0*h(6)-0.25d0*h(9)
         h(3)=0.25d0*(1.0d0-r)*(1.0d0-s)-0.5d0*h(6)-
     $        0.5d0*h(7)-0.25d0*h(9)
         h(4)=0.25d0*(1.0d0+r)*(1.0d0-s)-0.5d0*h(7)-
     $        0.5d0*h(8)-0.25d0*h(9)
      endif
      if (nis(ig) .eq. 4) then 
         h(1)=0.25d0*(1.0d0+r)*(1.0d0+s)
         h(2)=0.25d0*(1.0d0-r)*(1.0d0+s)
         h(3)=0.25d0*(1.0d0-r)*(1.0d0-s)
         h(4)=0.25d0*(1.0d0+r)*(1.0d0-s)
      endif
      p(1,9)=0.0d0
      if (nis(ig) .eq. 9) then 
         p(1,9)=-2.0d0*r*(1.0d0-s**2)
      endif
      if (nis(ig) .ne. 4) then      
         p(1,5)=-r*(1.0d0+s)-0.5d0*p(1,9)
         p(1,6)=-0.5d0*(1.0d0-s**2)-0.5d0*p(1,9)
         p(1,7)=-r*(1.0d0-s)-0.5d0*p(1,9)
         p(1,8)=0.5d0*(1.0d0-s**2)-0.5d0*p(1,9)
         p(1,1)=0.25d0*(1.0d0+s)-0.5d0*p(1,5)-
     $        0.5d0*p(1,8)-0.25d0*p(1,9)
         p(1,2)=-0.25d0*(1.0d0+s)-0.5d0*p(1,5)-
     $        0.5d0*p(1,6)-0.25d0*p(1,9)
         p(1,3)=-0.25d0*(1.0d0-s)-0.5d0*p(1,6)-
     $        0.5d0*p(1,7)-0.25d0*p(1,9)
         p(1,4)=0.25d0*(1.0d0-s)-0.5d0*p(1,7)-
     $        0.5d0*p(1,8)-0.25d0*p(1,9)
      endif
      if (nis(ig) .eq. 4) then
         p(1,1)=0.25d0*(1.0d0+s)
         p(1,2)=-0.25d0*(1.0d0+s)
         p(1,3)=-0.25d0*(1.0d0-s)
         p(1,4)=0.25d0*(1.0d0-s)
      endif
      p(2,9)=0.0d0
      if (nis(ig) .eq. 9) then 
         p(2,9)=-2.0d0*s*(1.0d0-r**2)     
      endif
      if (nis(ig) .ne. 4) then
         p(2,5)=0.5d0*(1.0d0-r**2)-0.5d0*p(2,9)
         p(2,6)=-s*(1.0d0-r)-0.5d0*p(2,9)
         p(2,7)=-0.5d0*(1.0d0-r**2)-0.5d0*p(2,9)
         p(2,8)=-s*(1.0d0+r)-0.5d0*p(2,9)
         p(2,1)=0.25d0*(1.0d0+r)-0.5d0*p(2,5)-
     $        0.5d0*p(2,8)-0.25d0*p(2,9)
         p(2,2)=0.25d0*(1.0d0-r)-0.5d0*p(2,5)-
     $        0.5d0*p(2,6)-0.25d0*p(2,9)
         p(2,3)=-0.25d0*(1.0d0-r)-0.5d0*p(2,6)-
     $        0.5d0*p(2,7)-0.25d0*p(2,9)
         p(2,4)=-0.25d0*(1.0d0+r)-0.5d0*p(2,7)-
     $        0.5d0*p(2,8)-0.25d0*p(2,9)
      endif
      if (nis(ig) .eq. 4) then
         p(2,1)=0.25d0*(1.0d0+r)
         p(2,2)=0.25d0*(1.0d0-r)
         p(2,3)=-0.25d0*(1.0d0-r)
         p(2,4)=-0.25d0*(1.0d0+r)         
      endif
      return
      end