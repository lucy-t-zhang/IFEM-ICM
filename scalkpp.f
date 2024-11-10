      subroutine scalkpp(xkpps,ocpp,k)
      implicit real*8 (a-h,o-z)
      include 'common'
      xkpps=ocpp*hp(k)*hp(k)
      return
      end