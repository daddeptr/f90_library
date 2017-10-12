modules.f90
type CAMBparams
!jv
         integer   :: Num_ISW_catalogs
         integer   :: Num_first_ISW
!ejv

cmbmain.f90
subroutine GetSourceMem
     
  if (CP%WantScalars) then
     if (CP%Dolensing) then
!jv
        Num_first_ISW = 4
        SourceNum = 3+Num_ISW_catalogs
!            SourceNum=3
            C_last = C_PhiTemp
         else
            Num_first_ISW = 3
            SourceNum = 2+Num_ISW_catalogs
!            SourceNum=2
            C_last = C_Cross
           end if
        else
            Num_first_ISW = 4
            SourceNum = 3+Num_ISW_catalogs
!           SourceNum=3 
   end if
!ejv


cmbmain.f90
subroutine ClTransferToCl
       if (CP%WantScalars) then
           lSamp = CTransS%ls
!jv
!           allocate(iCl_Scalar(CTransS%ls%l0,C_Temp:C_last,CP%InitPower%nn))
! For each ISW catalog we calculate C_l^{Tg} and C_l^{gg}
           allocate(iCl_Scalar(CTransS%ls%l0,C_Temp:C_last+2*Num_ISW_catalogs,CP%InitPower%nn))
!ejv
           iCl_scalar = 0

subroutine InterpolateCls

          if (CP%WantScalars) then
!jv
!               do i = C_Temp, C_last
               do i = C_Temp, C_last+2*Num_ISW_catalogs
!ejv              
                call InterpolateClArr(CTransS%ls,iCl_scalar(1,i,in),Cl_scalar(lmin, in, i),CTransS%ls%l0)
               end do
             end if

?

              
                !The unwrapped form is faster
     
!jv                 sums(1) = sums(1) + IV%Source_q(n,1)*J_l
!jv                 sums(2) = sums(2) + IV%Source_q(n,2)*J_l
!jv                 sums(3) = sums(3) + IV%Source_q(n,3)*J_l
!jv
                 sums(:) = sums(:) + IV%Source_q(n,:)*J_l
!ejv



!jv             if (.not. DebugEvolution .and. (EV%q*tauend > max_etak_scalar .and. tauend > taurend) &
!                  .and. .not. CP%Dolensing .and. &
!jv modified to ensure integration up to today!!!!!
      if (.not. DebugEvolution .and. (EV%q*tauend > max_etak_scalar .and. tauend > taurend) &
                 .and. .not. CP%Dolensing .and. Num_ISW_catalogs==0 &
                  (.not.CP%WantTransfer.or.tau > tautf(CP%Transfer%num_redshifts))) then
!ejv            



!jv Modified to ensure integration up to today
!                     if ((CP%Dolensing .or. IV%q*atau0(i) < max_etak_scalar) .and. xf > 1.e-8_dl) then
                     if ((CP%Dolensing .or. IV%q*atau0(i) < max_etak_scalar .or. Num_ISW_catalogs > 0) & 
                          .and. xf > 1.e-8_dl) then
!ejv

!jv Modified to ensure integration up to today
!      if ((num2*IntAccuracyBoost < dchisource .and. .not. CP%Dolensing) & !Oscillating fast 
!        .or. (nstart>IV%SourceSteps.and.nend>IV%SourceSteps)) then  
      if ((num2*IntAccuracyBoost < dchisource .and. .not. CP%Dolensing &
          .and. Num_ISW_catalogs==0) & !Oscillating fast 
        .or. (nstart>IV%SourceSteps.and.nend>IV%SourceSteps)) then  
!ejv





! jv The actual calculation of C_l in cmbmain.f90


          do q_ix = 1, CTrans%num_q_int 

             if (.not.(CP%closed.and.nint(CTrans%q_int(q_ix)*CP%r)<=CTrans%ls%l(j))) then 
               !cut off at nu = l + 1
             dlnk = dlnks(q_ix)
             apowers = pows(q_ix)

             iCl_scalar(j,C_Temp:C_E,pix) = iCl_scalar(j,C_Temp:C_E,pix) +  &
                          apowers*CTrans%Delta_p_l_k(1:2,j,q_ix)**2*dlnk
             iCl_scalar(j,C_Cross,pix) = iCl_scalar(j,C_Cross,pix) + &
                          apowers*CTrans%Delta_p_l_k(1,j,q_ix)*CTrans%Delta_p_l_k(2,j,q_ix)*dlnk
!jv
!             if (CTrans%NumSources>2) then
              if(C_last>C_Cross) then
!ejv
                        iCl_scalar(j,C_Phi,pix) = iCl_scalar(j,C_Phi,pix) +  &
                                                       apowers*CTrans%Delta_p_l_k(3,j,q_ix)**2*dlnk
                        iCl_scalar(j,C_PhiTemp,pix) = iCl_scalar(j,C_PhiTemp,pix) +  &
                                          apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(1,j,q_ix)*dlnk
             end if
!jv
! ISW-LSS correlation
             if (CP%Num_ISW_catalogs>0) then
!                do nISWcat=CP%Num_first_ISW, CP%Num_first_ISW+CP%Num_ISW_catalogs
                   ! gal-gal (gg)
                   iCl_scalar(j,C_last+1:C_last+CP%Num_ISW_catalogs,pix) = &
                          iCl_scalar(j,C_last+1:C_last+CP%Num_ISW_catalogs,pix) +  &
                   apowers*CTrans%Delta_p_l_k(CP%Num_first_ISW:CP%Num_first_ISW+CP%Num_ISW_catalogs,j,q_ix)**2*dlnk
                   ! gal-temp (gT)
                   iCl_scalar(j,C_last+CP%Num_ISW_catalogs+1:C_last+2*CP%Num_ISW_catalogs,pix) = &
                          iCl_scalar(j,C_last+CP%Num_ISW_catalogs+1:C_last+2*CP%Num_ISW_catalogs,pix) + &
                          apowers*CTrans%Delta_p_l_k(CP%Num_first_ISW:CP%Num_first_ISW+CP%Num_ISW_catalogs-1,j,q_ix)* &
                          CTrans%Delta_p_l_k(1,j,q_ix)*dlnk
!                end do
             end if
! end ISW-LSS correlation
!ejv             

             end if

           end do




equations.f90

in subroutine output
...
       else
         sources(3) = 0
       end if
      end if
 !        zz=1.0_dl/a - 1.0_dl
           zz=(1.0_dl-a)/a
           zbar = 0.50_dl
 !      zbar(1) = 0.1d0
 !      zbar(2) = 0.15d0
 !      zbar(3) = 0.3d0
 !      zbar(4) = 0.5d0
 !      zbar(5) = 0.9d0

              znot = zbar/1.41_dl
!jv
       if (CP%Num_ISW_catalogs>0) then
         do ISWcat=1, CP%Num_ISW_catalogs
! JPV 12/04/2006 This is now for matter C_l
!        if (CTransScal%NumSources > 2) then
       if (a.ge.1.0_dl) then
           sources(CP%Num_first_ISW + ISWcat - 1)=0.0_dl
        else
            if (ISWcat==2) sources(CP%Num_first_ISW + ISWcat - 1) = &
          adotoa/a * clxc * 1.5_dl*(1.41_dl**3)*(zz**2)/(zbar**3)*exp((-(1.41_dl*zz/zbar)**(1.5_dl)))
            if (ISWcat==1) sources(CP%Num_first_ISW + ISWcat - 1) = 0.0_dl
       end if
!              adotoa/a * clxc * ISWwindow(a,ISWcat)
                              !Tommaso writes this function and checks the minus sign etc.
       sources(CP%Num_first_ISW + ISWcat - 1) = &
            adotoa/a * clxc * 1.5_dl*(1.41_dl**3)*(zz**2)/(zbar**3)*exp((-(1.41_dl*zz/zbar)**(1.5_dl)))
         end do
       end if
!ejv          

! ------ :D
! --- Dec 2009
!!$            if (CTrans%NumSources>2) then
              if(C_last>C_Cross) then
! ---
                     iCl_scalar(j,C_Phi,pix)   =  &
                            iCl_scalar(j,C_Phi,pix)*fourpi*real(CTrans%ls%l(j)**2,dl)**2    
                     !The lensing power spectrum computed is l^4 C_l^{\phi\phi}
                     !We put pix extra factors of %l here to improve interpolation in CTrans%ls%l
                     iCl_scalar(j,C_PhiTemp,pix)   =  &
                            iCl_scalar(j,C_PhiTemp,pix)*fourpi*real(CTrans%ls%l(j)**2,dl)*CTrans%ls%l(j)
                      !Cross-correlation is CTrans%ls%l^3 C_l^{\phi T}
             end if
! ------ :D
! --- Dec 2009
! Normalization l(l+1)/2pi in both
             if (CP%Num_ISW_catalogs>0) then
!                do nISWcat=CP%Num_first_ISW, CP%Num_first_ISW+CP%Num_ISW_catalogs
                   ! gal-gal (gg)
                   iCl_scalar(j,C_last+1:C_last+CP%Num_ISW_catalogs,pix) = &
                          iCl_scalar(j,C_last+1:C_last+CP%Num_ISW_catalogs,pix) * dbletmp
                   ! gal-temp (gT)
                   iCl_scalar(j,C_last+CP%Num_ISW_catalogs+1:C_last+2*CP%Num_ISW_catalogs,pix) = &
                          iCl_scalar(j,C_last+CP%Num_ISW_catalogs+1:C_last+2*CP%Num_ISW_catalogs,pix) * dbletmp
!                end do
             end if
! ---
 
