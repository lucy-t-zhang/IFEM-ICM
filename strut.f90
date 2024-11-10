      subroutine strut(w,ne)
      implicit real*8 (a-h,o-z)     
      include 'common'
	ixmass=1
      if (nis(2) .eq. 2) then
            if (ixmass .eq. 2) then
               xfsk(3)=syoun*sarea/xl
               xfsk(4)=xfsk(3)
               xfsm(3)=sdensi*xl*sarea/3.0d0
               xfsm(4)=sdensi*xl*sarea/3.0d0
               nte1=nea(ne,1)
               nte2=nea(ne,2)
               xfu(3)=-xfsk(3)*dis(2,nte1,1)-&
                   xfsm(3)*acm(2,nte1,1)
               xfu(4)=-xfsk(4)*dis(2,nte2,1)-&
                   xfsm(3)*acm(2,nte1,1)
            endif
            if (ixmass .eq. 1) then
               xfsk(1)=syoun*sarea/xl
               xfsk(2)=xfsk(1)
               xfsm(1)=sdensi*xl*sarea/3.0d0
               xfsm(2)=sdensi*xl*sarea/3.0d0
               nte1=nea(ne,1)
               nte2=nea(ne,2)
               xfu(1)=-xfsk(1)*dis(1,nte1,1)-&
                   xfsm(1)*acm(1,nte1,1)
               xfu(2)=-xfsk(2)*dis(1,nte2,1)-&
                   xfsm(1)*acm(1,nte1,1) 
            endif
      else       
         if (nstr .ne. 3) then
            call material
            do ni=1,nis(2)
                  xfsm(ni)=xfsm(ni)+w*sdensi*h(ni)*h(ni)
                  xfsm(ni+nis(2))=xfsm(ni+nis(2))+&
                      w*sdensi*h(ni)*h(ni)
                  xfsk(ni)=xfsk(ni)+w*((bd(1,ni)*cmat(1,1)+&
                      bd(2,ni)*cmat(3,1))*bd(1,ni)+&
                      (bd(1,ni)*cmat(1,3)+bd(2,ni)*cmat(3,3))*bd(2,ni))
                  xfsk(ni+nis(2))=xfsk(ni+nis(2))+&
                      w*((bd(2,ni)*cmat(2,2)+&
                      bd(1,ni)*cmat(3,2))*bd(2,ni)+&
                      (bd(2,ni)*cmat(2,3)+bd(1,ni)*cmat(3,3))*bd(1,ni))
           enddo
         endif
      endif
      return
      end