      subroutine inter
      implicit real*8 (a-h,o-z)
      include 'common'
      do i=1,numct
         nmastd=2*(nodema(i)-1)+imasdir(i)
         nmast=id(nmastd)
         nslavd=2*(nodesl(i)-1)+islavdir(i)
         nslav=id(nslavd)
         drf(nmast)=drf(nmast)+amct(i)*drf(nslav)
         sk(nmast)=sk(nmast)+amct(i)*sk(nslav)
         sk(nmast)=sk(nmast)+amct(i)*sk(nslav)
	enddo
      return
      end