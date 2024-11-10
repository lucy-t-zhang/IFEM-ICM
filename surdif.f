      subroutine surdif
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension nu(60),ni(60),nj(60)
      logical jjx     
      do i=1,nnds
         nu(i)=neqv+i
      enddo
      if (nsmoos .eq. 2) then
         do i=1,nnds
            jjx=(vel(1,nnn(i),1) .eq. 0.0d0) .or. 
     $           (i .eq. 1) .or. (i .eq. nnds)
            if (jjx) then
               sk(nu(i))=sk(nu(i))+1.0d0
               drf(nu(i))=drf(nu(i))-
     $              (velm(2,nnn(i),1)-vel(2,nnn(i),1))
            else
               dx1=coor(nnn(i+1),1)-coor(nnn(i),1)
               du1=dis(2,nnn(i+1),1)-dis(2,nnn(i),1)
               dx2=coor(nnn(i-1),1)-coor(nnn(i),1)
               du2=dis(2,nnn(i-1),1)-dis(2,nnn(i),1)
               dx3=coor(nnn(i+1),1)-coor(nnn(i-1),1)
               du3=dis(2,nnn(i+1),1)-dis(2,nnn(i-1),1)
               if (vel(1,nnn(i),1) .lt. 0.0d0) then
                  sk(nu(i))=sk(nu(i))+1.0d0-
     $                 alpha*dt/beta/dx2*vel(1,nnn(i),1)
                  drf(nu(i))=drf(nu(i))-(velm(2,nnn(i),1)-
     $                 vel(2,nnn(i),1)+vel(1,nnn(i),1)*du2/dx2)
               endif
               if (vel(1,nnn(i),1) .gt. 0.0d0) then
                  sk(nu(i))=sk(nu(i))+1.0d0-
     $                 alpha*dt/beta/dx1*vel(1,nnn(i),1)
                  drf(nu(i))=drf(nu(i))-(velm(2,nnn(i),1)-
     $                 vel(2,nnn(i),1)+vel(1,nnn(i),1)*du1/dx1)
               endif
            endif
         enddo
      endif
      if (nsmoos .eq. 3) then
         do i=1,nnds
            jjx=(vel(1,nnn(i),1) .eq. 0.0d0) .or. 
     $           (i .eq. 1) .or. (i .eq. nnds)
            if (jjx) then
               sk(nu(i))=sk(nu(i))+1.0d0
               drf(nu(i))=drf(nu(i))-
     $              (velm(2,nnn(i),1)-vel(2,nnn(i),1))
            else
		  if (vel(1,nnn(i),1) .lt. 0.0d0) then
                  m1=i-2
                  m2=i-1
                  m3=i
                  m4=i+1
                  if (m1 .le. 0) then
                     m1=1
                  endif
                  if (m2 .le. 0) then
                     m2=1
                  endif
               endif
               if (vel(1,nnn(i),1) .gt. 0.0d0) then
                  m1=i+2
                  m2=i+1
                  m3=i
                  m4=i-1
                  if (m1 .gt. nnds) then
                     m1=nnds
                  endif
                  if (m2 .gt. nnds) then
                     m2=nnds
                  endif
               endif
               dx5=2.0d0*coor(nnn(m4),1)+
     $              3.0d0*coor(nnn(m3),1)-6.0d0*coor(nnn(m2),1)+
     $              coor(nnn(m1),1)
               du5=2.0d0*dis(2,nnn(m4),1)+
     $              3.0d0*dis(2,nnn(m3),1)-6.0d0*dis(2,nnn(m2),1)+
     $              dis(2,nnn(m1),1)
               sk(nu(m3))=sk(nu(m3))+1.0d0+
     $              3.0d0*alpha*dt/beta/dx5*vel(1,nnn(m3),1)
               drf(nu(m3))=drf(nu(m3))-(velm(2,nnn(m3),1)-
     $              vel(2,nnn(m3),1)+vel(1,nnn(m3),1)*du5/dx5)
            endif
          enddo
      endif
      if (nsmoos .eq. 0) then
         do i=2,nnds-1
            dx3=coor(nnn(i+1),1)-coor(nnn(i-1),1)
            du3=dis(2,nnn(i+1),1)-dis(2,nnn(i-1),1)
            sk(nu(i))=sk(nu(i))+1.0d0
            drf(nu(i))=drf(nu(i))-(velm(2,nnn(i),1)-
     $           vel(2,nnn(i),1)+vel(1,nnn(i),1)*du3/dx3)
         enddo
      endif
      return
      end