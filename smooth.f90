      subroutine smooth
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension nen(9),ptem(4,2)
      if (nsmfun .eq. 1) then     
         do ne=1,numele(1)
            do m=1,9
               nen(m)=nea(ne,m)
            enddo     
            do j=1,2
               ptem(1,j)=velm(j,nen(1),1)+velm(j,nen(9),1)-&
                   velm(j,nen(5),1)-velm(j,nen(8),1)
               ptem(2,j)=velm(j,nen(5),1)+velm(j,nen(6),1)-&
                   velm(j,nen(2),1)-velm(j,nen(9),1)
               ptem(3,j)=velm(j,nen(8),1)+velm(j,nen(7),1)-&
                   velm(j,nen(9),1)-velm(j,nen(4),1)
               ptem(4,j)=velm(j,nen(3),1)+velm(j,nen(9),1)-&
                   velm(j,nen(6),1)-velm(j,nen(7),1)
               velm(j,nen(1),1)=velm(j,nen(1),1)-&
                   0.25d0*smoo*ptem(1,j)
               velm(j,nen(5),1)=velm(j,nen(5),1)+&
                   0.25d0*smoo*ptem(1,j)
               velm(j,nen(9),1)=velm(j,nen(9),1)-&
                   0.25d0*smoo*ptem(1,j)
               velm(j,nen(8),1)=velm(j,nen(8),1)+&
                   0.25d0*smoo*ptem(1,j)
               velm(j,nen(5),1)=velm(j,nen(5),1)-&
                   0.25d0*smoo*ptem(2,j)
               velm(j,nen(2),1)=velm(j,nen(2),1)+&
                  0.25d0*smoo*ptem(2,j)
               velm(j,nen(6),1)=velm(j,nen(6),1)-&
                   0.25d0*smoo*ptem(2,j)
               velm(j,nen(9),1)=velm(j,nen(9),1)+&
                   0.25d0*smoo*ptem(2,j)
               velm(j,nen(8),1)=velm(j,nen(8),1)-&
                   0.25d0*smoo*ptem(3,j)
               velm(j,nen(9),1)=velm(j,nen(9),1)+&
                   0.25d0*smoo*ptem(3,j)
               velm(j,nen(7),1)=velm(j,nen(7),1)-&
                   0.25d0*smoo*ptem(3,j)
               velm(j,nen(4),1)=velm(j,nen(4),1)+&
                   0.25d0*smoo*ptem(3,j)
               velm(j,nen(9),1)=velm(j,nen(9),1)-&
                   0.25d0*smoo*ptem(4,j)
               velm(j,nen(6),1)=velm(j,nen(6),1)+&
                   0.25d0*smoo*ptem(4,j)
               velm(j,nen(3),1)=velm(j,nen(3),1)-&
                   0.25d0*smoo*ptem(4,j)
               velm(j,nen(7),1)=velm(j,nen(7),1)+&
                   0.25d0*smoo*ptem(4,j)
            enddo
		enddo
          do i=1,numgb
            if (ndirgb(i) .eq. 111111) then
               do j=1,numdir(i)
                  nl=nodegb(i,j)
                  velm(1,nl,1)=0.0d0
                  velm(2,nl,1)=0.0d0
			 enddo
            endif
            if (ndirgb(i) .eq. 110111) then
               do j=1,numdir(i)
                  nl=nodegb(i,j)
                  velm(1,nl,1)=0.0d0
			 enddo
            endif
            if (ndirgb(i) .eq. 101111) then
               do j=1,numdir(i)
                  nl=nodegb(i,j)
                  velm(2,nl,1)=0.0d0
			 enddo
            endif
		enddo
      endif
      if (nsmfun .eq. 2) then
         do ne=1,numele(1)
            do m=1,9
               nen(m)=nea(ne,m)
		  enddo     
            do j=1,2
               ptem(1,j)=dis(j,nen(1),1)+dis(j,nen(9),1)-&
                   dis(j,nen(5),1)-dis(j,nen(8),1)
               ptem(2,j)=dis(j,nen(5),1)+dis(j,nen(6),1)-&
                   dis(j,nen(2),1)-dis(j,nen(9),1)
               ptem(3,j)=dis(j,nen(8),1)+dis(j,nen(7),1)-&
                  dis(j,nen(9),1)-dis(j,nen(4),1)
               ptem(4,j)=dis(j,nen(3),1)+dis(j,nen(9),1)-&
                   dis(j,nen(6),1)-dis(j,nen(7),1)
               dis(j,nen(1),1)=dis(j,nen(1),1)-&
                   0.25d0*smoo*ptem(1,j)
               dis(j,nen(5),1)=dis(j,nen(5),1)+&
                   0.25d0*smoo*ptem(1,j)
               dis(j,nen(9),1)=dis(j,nen(9),1)-&
                  0.25d0*smoo*ptem(1,j)
               dis(j,nen(8),1)=dis(j,nen(8),1)+&
                   0.25d0*smoo*ptem(1,j)
               dis(j,nen(5),1)=dis(j,nen(5),1)-&
                  0.25d0*smoo*ptem(2,j)
               dis(j,nen(2),1)=dis(j,nen(2),1)+&
                   0.25d0*smoo*ptem(2,j)
               dis(j,nen(6),1)=dis(j,nen(6),1)-&
                   0.25d0*smoo*ptem(2,j)
               dis(j,nen(9),1)=dis(j,nen(9),1)+&
                   0.25d0*smoo*ptem(2,j)
               dis(j,nen(8),1)=dis(j,nen(8),1)-&
                  0.25d0*smoo*ptem(3,j)
               dis(j,nen(9),1)=dis(j,nen(9),1)+&
                   0.25d0*smoo*ptem(3,j)
               dis(j,nen(7),1)=dis(j,nen(7),1)-&
                   0.25d0*smoo*ptem(3,j)
               dis(j,nen(4),1)=dis(j,nen(4),1)+&
                   0.25d0*smoo*ptem(3,j)
               dis(j,nen(9),1)=dis(j,nen(9),1)-&
                   0.25d0*smoo*ptem(4,j)
               dis(j,nen(6),1)=dis(j,nen(6),1)+&
                   0.25d0*smoo*ptem(4,j)
               dis(j,nen(3),1)=dis(j,nen(3),1)-&
                   0.25d0*smoo*ptem(4,j)
               dis(j,nen(7),1)=dis(j,nen(7),1)+&
                   0.25d0*smoo*ptem(4,j)
		enddo
	enddo
      endif
      return
      end