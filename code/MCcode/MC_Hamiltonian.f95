!--------------------------------------------------------------------
!
!
! This subroutine calculates the field Hamiltonian from the phi values.
!      
!      by Quinn MacPherson based on code from Shifan Mao
!       Made a separate function on 7/8/16
!
!   If initialize then calculate all bins.
!   Otherwise calcualte only specified bins.
!-------------------------------------------------------------------

function phi_function(phi_s,VV,mc) result(dx)
use simMod
    implicit none
    double precision, intent(in) :: phi_s! franction solvent
    double precision, intent(in) :: VV   !Fraction of a bin
    TYPE(MCvar), intent(in) :: mc
    double precision dx ! change in energy / KAP

    double precision ratio  !geometric ratio
    integer indexPhi
    integer indexRatio
    double precision abovePhi,aboveRatio, temp
    
    ratio=sqrt(4*mc%V/(mc%DEL**2*3.1416*mc%L0))  ! rod diameter / DEL
    !ratio=sqrt(mc%V/(mc%DEL**2*3.1416*mc%L0))  ! rod radious / DEL
    
    temp=(1.0-PHI_s-mc%phiMin)*(mc%nPhiValues-1)/(mc%phiMax-mc%phiMin) +1.0
    indexPhi=floor(temp)
    abovePhi=temp - dble(indexPhi)

    temp=(ratio-mc%dMin)*(mc%nDiameters-1)/(mc%dMax-mc%dMin) + 1.0
    indexRatio=floor(temp);
    aboveRatio=temp - dble(indexRatio)


    dx=0.0
    dx=dx+mc%exclusionMu(indexRatio,indexPHI)*(1-aboveRatio)*(1-abovePhi) 
    dx=dx+mc%exclusionMu(indexRatio,indexPHI+1)*(aboveRatio)*(1-abovePhi) 
    dx=dx+mc%exclusionMu(indexRatio+1,indexPHI)*(1-aboveRatio)*(abovePhi) 
    dx=dx+mc%exclusionMu(indexRatio+1,indexPHI+1)*(aboveRatio)*(abovePhi)


    dx=dx*VV*(1.0-PHI_s)/(ratio**2)  ! energy per length to energy

    !if (PHI_s.gt.0.05) then 
    !    dx = VV*( PHI_s*log(PHI_s) -  PHI_s + 1.0)
    !elseif (PHI_s.gt.0.0) then 
    !    dx = VV*( 1.0-4.0*PHI_s)
    !else
    !    dx = 10.0
    !endif
    
end function phi_function
subroutine hamiltonian(mc,md,initialize)
use simMod
use setPrecision
implicit none
TYPE(MCvar), intent(inout) :: mc   ! <---- Contains output
TYPE(MCData), intent(in) :: md  
logical, intent(in) :: initialize ! Need to do all beads
double precision PHIPoly ! fraction polymer
double precision phi_A ! demsotu of A
double precision phi_B ! density of B
double precision phi_h ! strength of field
double precision VV ! volume of bin
integer I,J ! for looping
double precision PHI_s
double precision phi_function

mc%dx_Chi=0.0_dp
mc%Dx_Couple=0.0_dp
mc%Dx_Kap=0.0_dp
mc%Dx_Field=0.0_dp


if (initialize) then  ! calculate absolute energy
    if (mc%simType.eq.0) then ! Melt Hamiltonian
        do I=1,mc%NBIN
            VV=md%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            mc%Dx_Chi=mc%Dx_Chi+(VV/mc%V)*(md%PHIA(I)*md%PHIB(I))
            mc%Dx_Kap=mc%dx_Kap+(VV/mc%V)*((md%PHIA(I)+md%PHIB(I)-1.0_dp)**2)
            mc%Dx_Field=mc%dx_Field-md%PHIH(I)*md%PHIA(I)
        enddo        
    elseif(mc%simType.eq.1) then ! Chromatin Hamiltonian
        do I=1,mc%NBIN
            VV=md%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            PHIPoly=md%PHIA(I)+md%PHIB(I)
            mc%Dx_Chi=mc%Dx_Chi+(VV/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
            mc%Dx_Couple=mc%Dx_Couple+VV*(md%PHIA(I))**2
            if(PHIPoly.GT.1.0_dp) then
               mc%Dx_Kap=mc%Dx_Kap+(VV/mc%V)*(PHIPoly-1.0_dp)**2
            endif
        enddo
    elseif(mc%simType.eq.2) then
        do I=1,mc%NBIN
            VV=md%Vol(I)
            if (VV.le.0.1_dp) CYCLE
            PHIPoly=md%PHIA(I)+md%PHIB(I)
            PHI_s=1.0_dp-PHIPoly
            mc%Dx_Kap = mc%Dx_Kap + phi_function(PHI_s,VV,mc)
        enddo
    else
        print*, "Error in MC_int, simType",mc%simType, &
                " notdefined"
    endif
else ! Calculate change in energy
    if (mc%simType.eq.0) then ! Melt Hamiltonian
        do I=1,mc%NPHI
            J=md%INDPHI(I)
            VV=md%Vol(J)
            if (VV.le.0.1_dp) CYCLE
            ! new
            phi_A=md%PHIA(J)+md%DPHIA(I)
            phi_B=md%PHIB(J)+md%DPHIB(I)
            phi_h=md%PHIH(J)
            mc%Dx_Chi=mc%Dx_Chi+(VV/mc%V)*phi_A*phi_B
            mc%Dx_Kap=mc%Dx_Kap+(VV/mc%V)*((phi_A+phi_B-1.0_dp)**2)
            mc%Dx_Field=mc%Dx_Field-phi_h*phi_A
            ! minus old
            mc%Dx_Chi=mc%Dx_Chi-(VV/mc%V)*(md%PHIA(J)*md%PHIB(J))
            mc%Dx_Kap=mc%Dx_Kap-(VV/mc%V)*((md%PHIA(J)+md%PHIB(J)-1.0_dp)**2)
            mc%Dx_Field=mc%Dx_Field+phi_h*md%PHIA(J)
        enddo
    elseif(mc%simType.eq.1) then ! Chromatin Hamiltonian
        do I=1,mc%NPHI
            J=md%INDPHI(I)
            VV=md%Vol(J)
            if (VV.le.0.1_dp) CYCLE
            ! new ...
            PHIPoly=md%PHIA(J)+md%DPHIA(I)+md%PHIB(J)+md%DPHIB(I)
            mc%Dx_Chi=mc%Dx_Chi+(VV/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
            mc%Dx_Couple=mc%Dx_Couple+VV*(md%PHIA(J)+md%DPHIA(I))**2
            if(PHIPoly.GT.1.0_dp) then
               mc%Dx_Kap=mc%Dx_Kap+(VV/mc%V)*(PHIPoly-1.0_dp)**2
            endif
            ! minus old
            PHIPoly=md%PHIA(J)+md%PHIB(J)
            mc%Dx_Chi=mc%Dx_Chi-(VV/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
            mc%Dx_Couple=mc%Dx_Couple-VV*(md%PHIA(J))**2
            if(PHIPoly.GT.1.0_dp) then
               mc%Dx_Kap=mc%Dx_Kap-(VV/mc%V)*(PHIPoly-1.0_dp)**2
            endif 
        enddo
    elseif(mc%simType.eq.2) then
        do I=1,mc%NPHI
            J=md%INDPHI(I)
            VV=md%Vol(J)

            !! new ...
            if (VV.le.0.1_dp) CYCLE
            PHIPoly=md%PHIA(J)+md%DPHIA(I)+md%PHIB(J)+md%DPHIB(I)
            PHI_s=1.0_dp-PHIPoly
            mc%Dx_Kap = mc%Dx_Kap + phi_function(PHI_s,VV,mc)

            ! minus old
            PHIPoly=md%PHIA(J)+md%PHIB(J)
            PHI_s=1.0_dp-PHIPoly
            
            mc%Dx_Kap = mc%Dx_Kap - phi_function(PHI_s,VV,mc)
        enddo
    endif
endif
mc%dx_chi=mc%dx_chi*mc%CHI_ON
mc%dx_couple=mc%dx_couple*mc%Couple_ON
mc%dx_Kap=mc%dx_Kap*mc%KAP_ON

mc%DEChi=mc%Chi*        mc%dx_chi
mc%DECouple=mc%HP1_Bind*mc%dx_couple
mc%DEKap=mc%Kap*        mc%dx_Kap
mc%DEField=mc%h_A*      mc%dx_Field
!print*, "NPHI",mc%NPHI," DEKap",mc%DEKap,"dx_Kap", mc%dx_Kap 
RETURN
END subroutine
!---------------------------------------------------------------!
