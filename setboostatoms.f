      subroutine setboostatoms
      include 'moldy.h'
c
      integer*4 i, j, ii, ihe, ijk
      real*8 xij, yij, zij, rxij, ryij, rzij, r2
      real*8 xhea, yhea, zhea, rij
      integer*4 idya1, iikey

      ihenum = 0
      do i = 1, nn1
         if(atpe(i) .gt. -3.89)then
           ihenum = ihenum + 1
           iheid(ihenum) = i
         endif
      end do
      write(*,*) 'total number of atoms with Ep>-3.89eV is ', ihenum
c determine the surrounding atoms of these accelerating atoms
      do ihe = 1,ihenum
         ijk = iheid(ihe)
         xhea = x0(ijk)
         yhea = y0(ijk)
         zhea = z0(ijk)
c idya: the number of found near-by atoms around current He position   
          idya1 = idya

c loop over all the atoms to find the atoms close to the current
c He position, distance less than 3 angstrom         
          do i = 1, nn1
            xij = (x0(i)-xhea)*b0(1,1)
            yij = (y0(i)-yhea)*b0(2,2)
            zij = (z0(i)-zhea)*b0(3,3)
            if(xij.ge. b0(1,1)/2.0d0) xij = xij - b0(1,1)
            if(xij.lt.-b0(1,1)/2.0d0) xij = xij + b0(1,1)
            if(yij.ge. b0(2,2)/2.0d0) yij = yij - b0(2,2)
            if(yij.lt.-b0(2,2)/2.0d0) yij = yij + b0(2,2)
            if(zij.ge. b0(3,3)/2.0d0) zij = zij - b0(3,3)
            if(zij.lt.-b0(3,3)/2.0d0) zij = zij + b0(3,3)
            rij = sqrt(yij*yij + zij*zij)
c for loops, we consider the distance around one B
           if(abs(xij).le.0.5*b0(1,1) .and. rij .lt. 3.5d0)then
              iikey = 1
c idya_n(i): store the id of atoms near by He atoms
              do j = 1, idya1
                 iii = idya_n(j)
                 if(iii .eq. i)then
                   iikey = 0
                   exit
                 endif
              end do

              if(iikey .eq. 1)then
                idya = idya + 1
                idya_n(idya) = i
                xjp(idya) = x0(i)
                yjp(idya) = y0(i)
                zjp(idya) = z0(i)
              endif
            endif
         end do
      end do

      write(*,*) 'total accelerating atoms = ', idya

      do i = 1, idya
         ii = idya_n(i)
         x00 = xjp(i)*b0(1,1)
         y00 = yjp(i)*b0(2,2)
         z00 = zjp(i)*b0(3,3)
         write(*,*) x00, y00, z00, id(ii), ii
      end do

      open(unit=21,file='mdyn_region.psi',position='append')
      write(21,'(a)') 'ITEM: TIMESTEP'
      write(21,'(i9)') nst
      write(21,'(a)') 'ITEM: NUMBER OF ATOMS'
      write(21,'(i6)') idya
      write(21,'(a)') 'ITEM: BOX BOUNDS'
      write(21,'(2(f18.9))') -b0(1,1)/2.0d0, b0(1,1)/2.0d0
      write(21,'(2(f18.9))') -b0(2,2)/2.0d0, b0(2,2)/2.0d0
      write(21,'(2(f18.9))') -b0(3,3)/2.0d0, b0(3,3)/2.0d0
      write(21,'(a)') 'ITEM: ATOMS id x y z itype pe ke pre shr'
      do 12211 j = 1, idya
         iij = idya_n(j)
         xiij = x0(iij)*b0(1,1)
         yiij = y0(iij)*b0(2,2)
         ziij = z0(iij)*b0(3,3)
         write(21,'(i6,2x,3(f12.5),2x,i2,4(f12.5))')
     &   iij,xiij,yiij,ziij,id(iij),atpe(iij),ekin(iij),atpress(iij),
     &        atshear(iij)
12211 continue
      close(21)

      return
      end subroutine 
