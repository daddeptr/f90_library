   program prova_rhooa_exp

     USE LambdaGeneral
     USE ModelParams, ONLY: nstep

     IMPLICIT NONE

     real(dl), dimension(nstep,2) :: rho

     CALL Init_Rhooa(rho)

     write(34) rho

   end program prova_rhooa_exp
