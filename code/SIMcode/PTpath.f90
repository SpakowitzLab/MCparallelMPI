function chi_path(s) result(chi)
    use setPrecision
    implicit none
    double precision, intent(in) :: s
    double precision chi
    double precision chi_max
    if (.false.) then
        chi_max=2.00_dp
        if (s.lt.0.5_dp) then
            chi=0.0_dp
        else
            chi=chi_max*2.0_dp*(s-0.5_dp)
        endif
    else
        chi=s
    endif
end function chi_path
function h_path(s) result(h)
    use setPrecision
    implicit none
    double precision, intent(in) :: s
    double precision h
    double precision h_max
    if (.false.) then
        h_max=10.0_dp
        if (s.lt.0.5_dp) then
             h=h_max*s*2.0_dp
        else
             h=h_max*(1.0_dp-s)*2.0_dp
        endif
    else
        h=0.0_dp
    endif
end function h_path
function mu_path(s) result(mu)
    use setPrecision
    implicit none
    double precision, intent(in) :: s
    double precision mu
    mu=s
end function mu_path
function kap_path(s) result(kap)
    use setPrecision
    implicit none
    double precision, intent(in) :: s
    double precision kap
    kap=s*10.0_dp
end function kap_path
function hp1_bind_path(s) result(hp1_bind)
    use setPrecision
    implicit none
    double precision, intent(in) :: s
    double precision hp1_bind
    hp1_bind=s
end function hp1_bind_path
