c$Id: umani4.f,v 1.1 2007/08/15 22:42:22 rlt Exp $
      subroutine umani4

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2007: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Dummy user mesh manipulation set routine

c      Inputs:
c         prtu   - Flag, output results if true

c      Outputs:
c         none   - Users are responsible for generating outputs
c                  through common blocks, etc.  See programmer
c                  manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include 'umac1.h'

      logical  pcomp,prtu

      if (pcomp(uct,'man4',4)) then
c         uct = 'name'
      else

      end if

      end
