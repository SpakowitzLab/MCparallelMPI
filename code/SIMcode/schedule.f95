!---------------------------------------------------------------*
Subroutine strength_schedule(mc,inton)
    use setPrecision
    use simMod
    implicit none
    TYPE(MCvar), intent(inout) :: mc
    integer, intent(inout) :: inton

    if (mc%IND.LE.mc%NNOINT) then
        INTON=0
    else
        INTON=1
    endif
    if(mc%ind.lt.mc%N_KAP_ON) then
        mc%KAP_ON=0.0_dp
    else
        mc%KAP_ON=1.0_dp
    endif

    if(mc%ind.lt.mc%N_CHI_ON) then
!        PTON=.False.
        mc%CHI_ON=0.0_dp
    else
!        PTON=.True.
        mc%CHI_ON=1.0_dp
    endif
!    maximum=-28.0_dp    
!    if(ind.lt.101) then
!        Couple_ON=0.0_dp
!    elseif(ind.lt.111) then
!        Couple_ON=0.1_dp
!    elseif(ind.lt.121) then
!        Couple_ON=0.2_dp
!    elseif(ind.lt.131) then
!        Couple_ON=0.3_dp
!    elseif(ind.lt.141) then
!        Couple_ON=0.4_dp
!    elseif(ind.lt.151) then
!        Couple_ON=0.5_dp
!    elseif(ind.lt.161) then
!        Couple_ON=0.6_dp
!    elseif(ind.lt.171) then
!        Couple_ON=0.7_dp
!    elseif(ind.lt.181) then
!        Couple_ON=0.8_dp
!    elseif(ind.lt.191) then
!        Couple_ON=0.9_dp
!    elseif(ind.lt.201) then
!        Couple_ON=1.0_dp
!    elseif(ind.lt.211) then
!        Couple_ON=1.1_dp
!    elseif(ind.lt.221) then
!        Couple_ON=1.2_dp
!    elseif(ind.lt.231) then
!        Couple_ON=1.1_dp
!    elseif(ind.lt.241) then
!        Couple_ON=1.0_dp
!    elseif(ind.lt.251) then
!        Couple_ON=0.9_dp
!    elseif(ind.lt.261) then
!        Couple_ON=0.8_dp
!    elseif(ind.lt.271) then
!        Couple_ON=0.7_dp
!    endif
    return
end subroutine
!---------------------------------------------------------------*
