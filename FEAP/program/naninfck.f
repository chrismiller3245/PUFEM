!$Id:$
      logical function naninfck(y,sizey,isw)

!     * * F E A P * * A Finite Element Analysis Program

!.... Copyright (c) 1984-2016: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2006
!-----[--.----+----.----+----.-----------------------------------------]
!     Purpose: Check if double precision number is NAN or -INF or +INF
!              This is Fortran  replacement for C: isinf
!                                                  isnan
!       Inputs:
!          y(*)     : The array of values to be checked
!          sizey    : Length of list
!          isw      : Switch - 0: Initialize values; >0: Check values

!       Outputs:
!         naninfck  : A logical variable which is true or false
!-----[--.----+----.----+----.-----------------------------------------]
      use, intrinsic :: ieee_arithmetic

      implicit none

      real*8   y(*)
      integer  sizey,isw,ic

      save

!     Set values for checks

      if(isw.eq.0) then

        naninfck = .true.

!     Check values

      else

        naninfck = .false.
        do ic = 1,sizey

          if(ieee_is_finite(y(ic))) then ! Finite value check
            continue
          elseif ( ieee_is_nan(y(ic)) ) then  ! IEEE interface
            naninfck = .true.
            return
          else
            naninfck = .true.
            return
          end if

        end do ! ic

      endif ! isw

      end
