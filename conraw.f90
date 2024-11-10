      subroutine conraw(nloc,xfsktem,xfutem)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension x(2,9),xj(2,2),xji(2,2),rs(2)
      dimension pt(12,2),st(12,2),xm(12)
      dimension nloc(4),xfsktem(4,4)
      dimension xfutem(8),utem(12,2),htem(12,4)
      dimension xcoe(4,4),ycoe(4,4)
      nistem=nis(ig)
      do 28 m=1,4
         xfutem(m)=0.0d0
         xfutem(m+4)=0.0d0
         do 24 n=1,4
            xfsktem(m,n)=0.0d0
 24      continue
 28   continue
      if (nis(ig) .eq. 4) then
         if (nale .eq. 1) then
            do 331 nos=1,4
               do 332 j=1,2
                  ntem=nloc(nos)
                  x(j,nos)=coor(ntem,j)+dis(j,ntem,1)
 332           continue
 331        continue
         else
            do 431 nos=1,4
               do 432 j=1,2
                  ntem=nloc(nos)
                  x(j,nos)=coor(ntem,j)
 432           continue
 431        continue
         endif
      endif
      if (nis(ig) .eq. 9) then
         if (nale .eq. 1) then
            do 1 nos=1,4
               do 2 j=1,2
                  ntem=nloc(nos)
                  x(j,nos)=coor(ntem,j)+dis(j,ntem,1)
 2             continue
 1          continue
            nis(1)=4
         else
            do 301 nos=1,4
               do 302 j=1,2
                  ntem=nloc(nos)
                  x(j,nos)=coor(ntem,j)
 302           continue
 301        continue
            nis(1)=4
         endif
      endif
      do 5 i=1,2
         pt(1,i)=0.5d0*x(i,1)+0.5d0*x(i,2)
         pt(2,i)=0.5d0*x(i,2)+0.5d0*x(i,3)
         pt(3,i)=0.5d0*x(i,3)+0.5d0*x(i,4)
         pt(4,i)=0.5d0*x(i,4)+0.5d0*x(i,1)
         pt(5,i)=0.25d0*x(i,1)+0.25d0*x(i,2)+&
             0.25d0*x(i,3)+0.25d0*x(i,4)
         pt(6,i)=x(i,1)
         pt(7,i)=x(i,2)
         pt(8,i)=x(i,3)
         pt(9,i)=x(i,4)
 5    continue
      do 6 i=1,4
         st(i,1)=pt(5,2)-pt(i,2)
         st(i,2)=-(pt(5,1)-pt(i,1))
 6    continue
      st(5,1)=pt(1,2)-pt(6,2)
      st(5,2)=-(pt(1,1)-pt(6,1))
      st(6,1)=pt(7,2)-pt(1,2)
      st(6,2)=-(pt(7,1)-pt(1,1))     
      st(7,1)=pt(2,2)-pt(7,2)
      st(7,2)=-(pt(2,1)-pt(7,1))
      st(8,1)=pt(8,2)-pt(2,2)
      st(8,2)=-(pt(8,1)-pt(2,1))     
      st(9,1)=pt(3,2)-pt(8,2)
      st(9,2)=-(pt(3,1)-pt(8,1))
      st(10,1)=pt(9,2)-pt(3,2)
      st(10,2)=-(pt(9,1)-pt(3,1))     
      st(11,1)=pt(4,2)-pt(9,2)
      st(11,2)=-(pt(4,1)-pt(9,1))
      st(12,1)=pt(6,2)-pt(4,2)
      st(12,2)=-(pt(6,1)-pt(4,1))
      do 7 ntt=1,12
         rs(1)=xrs(ntt,1)
         rs(2)=xrs(ntt,2)
         call element(rs)
         call jacob(x,xj,xji,det)
         call bdpd(xji)
         do 39 ni=1,4
            htem(ntt,ni)=h(ni)
 39      continue
         call convel(nloc)
         xm(ntt)=fdensi*(st(ntt,1)*fvel(1)+&
             st(ntt,2)*fvel(2))
 7    continue
      tt=0.0d0
      do 70 i=1,12
         tt=tt+dabs(xm(i))
 70   continue
      if (nsmoo .eq. 5) then
         do 10 i=1,4
            do 11 j=1,4
               xcoe(i,j)=0.0d0
               ycoe(i,j)=0.0d0
 11         continue
 10      continue
         do 12 i=1,4
            if (xm(i) .lt. 0.0d0) then
               f=xm(nsip(i,1))/xm(i)
               if (f .gt. 1.0d0) then
                  f=1.0d0
               else
                  if (f .lt. 0.0d0) then
                     f=0.0d0
                  endif
               endif
               xcoe(i,i)=1.0d0
               xcoe(i,nsip(i,1))=-f
               ycoe(i,nsco(i,1))=(1.0d0-f)
            endif     
            if (xm(i) .gt. 0.0d0) then
               f=xm(nsip(i,2))/xm(i)
               if (f .gt. 1.0d0) then
                  f=1.0d0
               else
                  if (f .lt. 0.0d0) then
                     f=0.0d0
                  endif
               endif
               xcoe(i,i)=1.0d0
               xcoe(i,nsip(i,2))=-f
               ycoe(i,nsco(i,2))=(1.0d0-f)
            endif
            if (xm(i) .eq. 0.0d0) then
               xcoe(i,i)=1.0d0
               ycoe(i,i)=1.0d0
            endif
 12      continue
         call gaussj(xcoe,4,4,ycoe,4,4)
         do 91 j=1,2
            do 99 i=1,4
               utem(i,j)=0.0d0
               do 92 m=1,4
                  utem(i,j)=utem(i,j)+ycoe(i,m)*&
                      vel(j,nloc(m),1)
 92            continue
 99         continue
 91      continue
         do 391 j=1,2
            do 399 i=5,12
               utem(i,j)=0.0d0
               do 392 m=1,4
                  utem(i,j)=utem(i,j)+htem(i,m)*&
                     vel(j,nloc(m),1)
 392            continue
 399         continue
 391      continue
      else
         do 31 j=1,2
            do 32 i=1,12
               utem(i,j)=0.0d0
               do 35 m=1,4
                  utem(i,j)=utem(i,j)+htem(i,m)*&
                      vel(j,nloc(m),1)
 35            continue
 32         continue
 31      continue
      endif
      do 22 j=1,2
         nj=(j-1)*4+1
         xfutem(nj)=xfutem(nj)+&
             utem(1,j)*xm(1)-utem(4,j)*xm(4)+&
             utem(5,j)*xm(5)+utem(12,j)*xm(12)
 22   continue
      do 94 j=1,2
         nj=(j-1)*4+2
         xfutem(nj)=xfutem(nj)+&
             utem(2,j)*xm(2)-utem(1,j)*xm(1)+&
             utem(7,j)*xm(7)+utem(6,j)*xm(6)
 94   continue
      do 95 j=1,2
         nj=(j-1)*4+3
         xfutem(nj)=xfutem(nj)+&
             utem(3,j)*xm(3)-utem(2,j)*xm(2)+&
             utem(9,j)*xm(9)+utem(8,j)*xm(8)
 95   continue
      do 96 j=1,2
         nj=(j-1)*4+4
         xfutem(nj)=xfutem(nj)+&
             utem(4,j)*xm(4)-utem(3,j)*xm(3)+&
             utem(11,j)*xm(11)+utem(10,j)*xm(10)
 96   continue
      if (nsmoo .eq. 5) then
         do 13 m=1,4
            xfsktem(1,m)=xfsktem(1,m)+htem(5,m)*xm(5)
            xfsktem(1,m)=xfsktem(1,m)+ycoe(1,m)*xm(1)
            xfsktem(1,m)=xfsktem(1,m)-ycoe(4,m)*xm(4)
            xfsktem(1,m)=xfsktem(1,m)+htem(12,m)*xm(12)
            xfsktem(2,m)=xfsktem(2,m)+htem(7,m)*xm(7)
            xfsktem(2,m)=xfsktem(2,m)+ycoe(2,m)*xm(2)
            xfsktem(2,m)=xfsktem(2,m)-ycoe(1,m)*xm(1)
            xfsktem(2,m)=xfsktem(2,m)+htem(6,m)*xm(6)
            xfsktem(3,m)=xfsktem(3,m)+htem(9,m)*xm(9)
            xfsktem(3,m)=xfsktem(3,m)+ycoe(3,m)*xm(3)
            xfsktem(3,m)=xfsktem(3,m)-ycoe(2,m)*xm(2)
            xfsktem(3,m)=xfsktem(3,m)+htem(8,m)*xm(8)
            xfsktem(4,m)=xfsktem(4,m)+htem(11,m)*xm(11)
            xfsktem(4,m)=xfsktem(4,m)+ycoe(4,m)*xm(4)
            xfsktem(4,m)=xfsktem(4,m)-ycoe(3,m)*xm(3)
            xfsktem(4,m)=xfsktem(4,m)+htem(10,m)*xm(10)
 13      continue
      else
         do 33 m=1,4
            xfsktem(1,m)=xfsktem(1,m)+htem(5,m)*xm(5)
            xfsktem(1,m)=xfsktem(1,m)+htem(1,m)*xm(1)
            xfsktem(1,m)=xfsktem(1,m)-htem(4,m)*xm(4)
            xfsktem(1,m)=xfsktem(1,m)+htem(12,m)*xm(12)
            xfsktem(2,m)=xfsktem(2,m)+htem(7,m)*xm(7)
            xfsktem(2,m)=xfsktem(2,m)+htem(2,m)*xm(2)
            xfsktem(2,m)=xfsktem(2,m)-htem(1,m)*xm(1)
            xfsktem(2,m)=xfsktem(2,m)+htem(6,m)*xm(6)
            xfsktem(3,m)=xfsktem(3,m)+htem(9,m)*xm(9)
            xfsktem(3,m)=xfsktem(3,m)+htem(3,m)*xm(3)
            xfsktem(3,m)=xfsktem(3,m)-htem(2,m)*xm(2)
            xfsktem(3,m)=xfsktem(3,m)+htem(8,m)*xm(8)
            xfsktem(4,m)=xfsktem(4,m)+htem(11,m)*xm(11)
            xfsktem(4,m)=xfsktem(4,m)+htem(4,m)*xm(4)
            xfsktem(4,m)=xfsktem(4,m)-htem(3,m)*xm(3)
            xfsktem(4,m)=xfsktem(4,m)+htem(10,m)*xm(10)
 33      continue
      endif
      nis(ig)=nistem
      return
      end