      subroutine cal_displace
      include 'moldy.h'
c
      integer*4 i, j, ii, ihe, ijk
      real*8 xij, yij, zij, rxij, ryij, rzij, r2
      real*8 xhea, yhea, zhea, rij
      integer*4 idya1, iikey

      do i = 1, nboost
         displace(i) = 0.0
      end do

      do i = 1, idya
         idi = idya_n(i)
         xi = x0(idi)
         yi = y0(idi)
         zi = z0(idi)
         xj = xjp(i)
         yj = yjp(i)
         zj = zjp(i)
         xij = (xi - xj)*b0(1,1)
         yij = (yi - yj)*b0(2,2)
         zij = (zi - zj)*b0(3,3)

         if(xij > b0(1,1)/2.0d0) xij = xij - b0(1,1)
         if(xij <-b0(1,1)/2.0d0) xij = xij + b0(1,1)
         if(yij > b0(2,2)/2.0d0) yij = yij - b0(2,2)
         if(yij <-b0(2,2)/2.0d0) yij = yij + b0(2,2)
         if(zij > b0(3,3)/2.0d0) zij = zij - b0(3,3)
         if(zij <-b0(3,3)/2.0d0) zij = zij + b0(3,3)

         dist = sqrt(xij*xij + yij*yij + zij*zij)
         displace(i) = dist
      end do

      return
      end subroutine 
