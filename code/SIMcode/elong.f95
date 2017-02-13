! returns .true. on search if RP is outside cylinder
Subroutine elong(RP,RCylinder,LCylinder,search)

    IMPLICIT NONE
    DOUBLE PRECISION, intent(in) :: RP(3)
    DOUBLE PRECISION, intent(in) :: RCylinder
    DOUBLE PRECISION, intent(in) :: LCylinder
    Logical, intent(out) :: search
    DOUBLE PRECISION r(3)

    search=.False.
    if( (RP(2)-RCylinder)**2 + (RP(3)-RCylinder)**2 .GT. Rcylinder**2) then
        search=.true.
        return
    endif
    if( RP(1) .LT. RCylinder) then
        r(1) = RP(1) - RCylinder
        r(2) = RP(2) - RCylinder
        r(3) = RP(3) - RCylinder
        if (r(1)**2+r(2)**2+r(3)**2 .gt. RCylinder**2) then
            search=.true.
            return
        endif
    elseif ( RP(1) .GT. RCylinder+LCylinder) then
        r(1) = RP(1) - RCylinder - LCylinder
        r(2) = RP(2) - RCylinder
        r(3) = RP(3) - RCylinder
        if (r(1)**2+r(2)**2+r(3)**2 .gt. RCylinder**2) then
            search=.true.
            return
        endif
    endif
end subroutine
