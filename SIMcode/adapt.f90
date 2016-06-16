!
!       Code for MC move addaptations when you need to 
!       adapt both move size and number of beads.
!
!           Quinn MacPherson
!

Subroutine MCvar_adapt(mc,MCTYPE,ISTEP,rand_stat)
! Run this after say 1000 move in order to improve performance
    !use mt19937, only : grnd
    use simMod
    use mersenne_twister 
    IMPLICIT NONE
    INTEGER ISTEP    ! Step number, need to initialize if 1
    TYPE(MCvar) mc
    INTEGER MCTYPE   ! Type of move
    Double Precision floatWindow  !like window but a floating point
    real urand(1)
    type(random_stat) rand_stat  !for random number generator
    mc%PHIT(MCTYPE)=real(mc%SUCCESS(MCTYPE))/real(mc%NADAPT(MCTYPE))
    mc%SUCCESS(MCTYPE)=0

    ! If move type has no amplitude then only ajust window
    if (MCTYPE.eq.7) then
        if (mc%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
           mc%WINDOW(MCTYPE)=mc%WINDOW(MCTYPE)+1
        else
           mc%WINDOW(MCTYPE)=mc%WINDOW(MCTYPE)-1
        endif
        !window limits
        if (mc%WINDOW(MCTYPE).LT.1) then
           mc%WINDOW(MCTYPE)=1
        elseif (mc%WINDOW(MCTYPE).GT.mc%MAXWINDOW(MCTYPE)) then
           mc%WINDOW(MCTYPE)=mc%MAXWINDOW(MCTYPE)
        endif
        return
    endif

    ! Adjust Amplidtude
    if (mc%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MCAMP(MCTYPE)*1.05_dp
    else
       mc%MCAMP(MCTYPE)=mc%MCAMP(MCTYPE)*0.95_dp
    endif
    
    ! If move has no window it doesn't need to be ajusted
    if ((MCTYPE.eq.4) .or. & 
        (MCTYPE.eq.5) .or. &
        (MCTYPE.eq.6)) then
        ! amplitude limits
        if (mc%MCAMP(MCTYPE).GT.mc%MAXAMP(MCTYPE)) then
           mc%MCAMP(MCTYPE)=mc%MAXAMP(MCTYPE)
        elseif (mc%MCAMP(MCTYPE).LT.mc%MINAMP(MCTYPE)) then
           mc%MCAMP(MCTYPE)=mc%MINAMP(MCTYPE)
        endif
        return
    endif

    if (ISTEP.eq.1) then
        ! initialize
        mc%WA_ratio(MCTYPE)=mc%WINDOW(MCTYPE)/mc%MCAMP(MCTYPE)
    endif
    !drift to chosen ratio
    mc%WA_ratio(MCTYPE)=mc%WA_ratio(MCTYPE)*0.96_dp+0.04_dp*mc%moveSlope(MCTYPE)

    ! probabalisticaly round window
    floatWindow=mc%WA_ratio(MCTYPE)*mc%MCAMP(MCTYPE)
    call random_number(urand,rand_stat) 
    if (urand(1).lt.(floatWindow-floor(floatWindow))) then
        mc%Window(MCTYPE)=CEILING(mc%WA_ratio(MCTYPE)*mc%MCAMP(MCTYPE))
    else
        mc%Window(MCTYPE)=floor(mc%WA_ratio(MCTYPE)*mc%MCAMP(MCTYPE))
    endif

    ! amplitude limits
    if (mc%MCAMP(MCTYPE).GT.mc%MAXAMP(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MAXAMP(MCTYPE)
    elseif (mc%MCAMP(MCTYPE).LT.mc%MINAMP(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MINAMP(MCTYPE)
    endif
    
    !window limits
    if (mc%WINDOW(MCTYPE).LT.1) then
       mc%WINDOW(MCTYPE)=1
    elseif (mc%WINDOW(MCTYPE).GT.mc%MAXWINDOW(MCTYPE)) then
       mc%WINDOW(MCTYPE)=mc%MAXWINDOW(MCTYPE)
    endif

    mc%WA_ratio(MCTYPE)=dble(mc%WINDOW(MCTYPE))/mc%MCAMP(MCTYPE)   
end subroutine
