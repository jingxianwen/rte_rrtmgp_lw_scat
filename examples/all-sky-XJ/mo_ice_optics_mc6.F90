module mo_ice_optics_mc6

   use mo_rte_kind, only: wp, wl
   use mo_mc6_table, only: absice4, extice4, ssaice4, asyice4, &
                        absice5, extice5, ssaice5, asyice5, &
                        absice6, extice6, ssaice6, asyice6   
   implicit none

   logical(wl), parameter :: flag_mc6  = .true.  ! T: use MC6 ice optics
                                                 ! F: use RRTMGP ice optics
   logical(wl), parameter :: flag_scat = .false. ! ONLY valid when flag_mc6=T
                                                 ! T: include scattering extinction
                                                 ! F: use absorption extinction only (ssa=g=0)
   public
contains

!==============================================================================
subroutine compute_all_from_mc6(ncol,nlay,nlwbands,icemsk,iciwp,reice, &
                tau,tau_w, tau_w_g)

  !-- input variables
   integer, intent(in)     :: ncol, nlay, nlwbands
   logical(wl), intent(in) :: icemsk(ncol,nlay)
   real(wp), intent(in)    :: iciwp(ncol,nlay), & !in cloud ice water path (kg/m2)
                              reice(ncol,nlay)    !ice effective radius (microns)
  !-- output variables
   real(wp),intent(out) :: tau    (ncol,nlay,nlwbands) ! extinction optical depth
   real(wp),intent(out) :: tau_w  (ncol,nlay,nlwbands) ! single scattering albedo * tau
   real(wp),intent(out) :: tau_w_g(ncol,nlay,nlwbands) ! assymetry parameter * tau * w

!--- Local ---!
   integer  :: ib,i,k
   real(wp) :: diaice               ! cloud ice effective diameter (microns)
   real(wp) :: diaice_max           !uplimit of diameter.
   real(wp) :: invrad
   real(wp) :: extcoice(nlwbands)   ! ice mass-extinction coefficients of TAMU scheme (m2/g)
   real(wp) :: ssacoice(nlwbands)   ! ice single scattering albedo of TAMU scheme (unitless)
   real(wp) :: asycoice(nlwbands)   ! ice asymmetric factor of TAMU scheme (unitless)
   real(wp) :: taucloud, ssacloud, asycloud, fpeak
   integer,parameter :: iceflag=3
       !  iceflag= 1 : MC6 ice cloud optical properties (MC6 ice crystal shape and 
       !                  0.1 variance gamma PSD, absorption and scattering)
       !           2 : THM ice cloud optical properties (Two-Habit model ice crystal 
       !                  shape and 0.1 variance gamma PSD, absorption and scattering)
       !           3 : TAMU MC6 ice cloud optical properties (MC6 ice crystal shape 
       !                  and Bryan Baum in-situ PSD, absorption and scattering)   
   real(wp), parameter :: gravit=9.80616_wp

   !----------------------------
   ! initialize return arrays
   !----------------------------   
   tau(:,:,:) = 0._wp
   tau_w(:,:,:) = 0._wp
   tau_w_g(:,:,:) = 0._wp

   !--------------------------------------------------------
   ! compute absorption coefficient and cloud optical depth
   !--------------------------------------------------------
   do i = 1,ncol
      do k = 1,nlay
        diaice = 2._wp*reice(i,k)
        extcoice(:) = 0._wp
        ssacoice(:) = 0._wp
        asycoice(:) = 0._wp

        !*** if there is no ice cloud ***
        !if( iciwp(i,k) < 1.e-80_wp .or. diaice .eq. 0._wp) then
        if(icemsk(i,k)) then

        !*** if there is an ice cloud layer ***    
        !-- set upper bound according to ice optical option -->
        !   iceflag=1. MC6 scattering + 0.1 variance Gamma PSD
        !           2. THM scattering + 0.1 variance Gamma PSD
        !           3. TAMU MC6 scattering + Bryan Baum in-situ PSD
          if (iceflag==1) then
             diaice_max=500._wp
          else if (iceflag==2) then
             diaice_max=500._wp
          else if (iceflag==3) then
             diaice_max=370._wp
          else
            print*,'ERROR in UMRad: current iceflag is not available.'
            stop !('ERROR in UMRad: current iceflag is not available.')
          end if

          diaice = min(diaice, diaice_max)

          if (diaice .ge. 3.0_wp .and. diaice .le. diaice_max) then
             if (iceflag==1) then
                invrad = 1.0/diaice - 0.05_wp
                do ib=1,nlwbands
                   if (diaice .ge. 20.0_wp) then
                      extcoice(ib) = ((extice4(1,1,ib) * invrad + &
                                       extice4(1,2,ib)) * invrad + &
                                       extice4(1,3,ib)) * invrad + &
                                       extice4(1,4,ib)
                      ssacoice(ib) = ((ssaice4(1,1,ib) * invrad + &
                                       ssaice4(1,2,ib)) * invrad + &
                                       ssaice4(1,3,ib)) * invrad + &
                                       ssaice4(1,4,ib)
                      asycoice(ib) = ((asyice4(1,1,ib) * invrad + &
                                       asyice4(1,2,ib)) * invrad + &
                                       asyice4(1,3,ib)) * invrad + &
                                       asyice4(1,4,ib)
                   else
                      extcoice(ib) = (((((extice4(2,1,ib) * invrad + &
                                           extice4(2,2,ib)) * invrad + &
                                           extice4(2,3,ib)) * invrad + &
                                           extice4(2,4,ib)) * invrad + &
                                           extice4(2,5,ib)) * invrad + &
                                           extice4(2,6,ib)) * invrad + &
                                           extice4(2,7,ib)
                      ssacoice(ib) = (((((ssaice4(2,1,ib) * invrad + &
                                           ssaice4(2,2,ib)) * invrad + &
                                           ssaice4(2,3,ib)) * invrad + &
                                           ssaice4(2,4,ib)) * invrad + &
                                           ssaice4(2,5,ib)) * invrad + &
                                           ssaice4(2,6,ib)) * invrad + &
                                           ssaice4(2,7,ib)
                      asycoice(ib) = (((((asyice4(2,1,ib) * invrad + &
                                           asyice4(2,2,ib)) * invrad + &
                                           asyice4(2,3,ib)) * invrad + &
                                           asyice4(2,4,ib)) * invrad + &
                                           asyice4(2,5,ib)) * invrad + &
                                           asyice4(2,6,ib)) * invrad + &
                                           asyice4(2,7,ib)
                   endif  ! end if of diaice .ge. 25.
                 enddo    ! end do of lwbands for cloud radiative coefficients
             else if (iceflag==2) then
                 invrad = 1.0/diaice - 0.05_wp
                 do ib=1,nlwbands
                   if (diaice .ge. 20.0_wp) then
                      extcoice(ib) = ((extice5(1,1,ib) * invrad + &
                                       extice5(1,2,ib)) * invrad + &
                                       extice5(1,3,ib)) * invrad + &
                                       extice5(1,4,ib)
                      ssacoice(ib) = ((ssaice5(1,1,ib) * invrad + &
                                       ssaice5(1,2,ib)) * invrad + &
                                       ssaice5(1,3,ib)) * invrad + &
                                       ssaice5(1,4,ib)
                      asycoice(ib) = ((asyice5(1,1,ib) * invrad + &
                                       asyice5(1,2,ib)) * invrad + &
                                       asyice5(1,3,ib)) * invrad + &
                                       asyice5(1,4,ib)

                   else
                      extcoice(ib) = (((((extice5(2,1,ib) * invrad + &
                                           extice5(2,2,ib)) * invrad + &
                                           extice5(2,3,ib)) * invrad + &
                                           extice5(2,4,ib)) * invrad + &
                                           extice5(2,5,ib)) * invrad + &
                                           extice5(2,6,ib)) * invrad + &
                                           extice5(2,7,ib)
                      ssacoice(ib) = (((((ssaice5(2,1,ib) * invrad + &
                                           ssaice5(2,2,ib)) * invrad + &
                                           ssaice5(2,3,ib)) * invrad + &
                                           ssaice5(2,4,ib)) * invrad + &
                                           ssaice5(2,5,ib)) * invrad + &
                                           ssaice5(2,6,ib)) * invrad + &
                                           ssaice5(2,7,ib)
                      asycoice(ib) = (((((asyice5(2,1,ib) * invrad + &
                                           asyice5(2,2,ib)) * invrad + &
                                           asyice5(2,3,ib)) * invrad + &
                                           asyice5(2,4,ib)) * invrad + &
                                           asyice5(2,5,ib)) * invrad + &
                                           asyice5(2,6,ib)) * invrad + &
                                           asyice5(2,7,ib)
                   endif  ! end if of diaice .ge. 20.
                 enddo    ! end do of lwbands for cloud radiative coefficients
             else if (iceflag==3) then
                 invrad = 1.0/diaice - 0.04_wp
                 do ib=1,nlwbands
                   if (diaice .ge. 25.0_wp) then
                      extcoice(ib) = ((extice6(1,1,ib) * invrad + &
                                       extice6(1,2,ib)) * invrad + &
                                       extice6(1,3,ib)) * invrad + &
                                       extice6(1,4,ib)
                      ssacoice(ib) = ((ssaice6(1,1,ib) * invrad + &
                                       ssaice6(1,2,ib)) * invrad + &
                                       ssaice6(1,3,ib)) * invrad + &
                                       ssaice6(1,4,ib)
                      asycoice(ib) = ((asyice6(1,1,ib) * invrad + &
                                       asyice6(1,2,ib)) * invrad + &
                                       asyice6(1,3,ib)) * invrad + &
                                       asyice6(1,4,ib)

                   else
                      extcoice(ib) = ((((((extice6(2,1,ib) * invrad + &
                                           extice6(2,2,ib)) * invrad + &
                                           extice6(2,3,ib)) * invrad + &
                                           extice6(2,4,ib)) * invrad + &
                                           extice6(2,5,ib)) * invrad + &
                                           extice6(2,6,ib)) * invrad + &
                                           extice6(2,7,ib)) * invrad + &
                                           extice6(2,8,ib)
                      ssacoice(ib) = ((((((ssaice6(2,1,ib) * invrad + &
                                           ssaice6(2,2,ib)) * invrad + &
                                           ssaice6(2,3,ib)) * invrad + &
                                           ssaice6(2,4,ib)) * invrad + &
                                           ssaice6(2,5,ib)) * invrad + &
                                           ssaice6(2,6,ib)) * invrad + &
                                           ssaice6(2,7,ib)) * invrad + &
                                           ssaice6(2,8,ib)
                      asycoice(ib) = ((((((asyice6(2,1,ib) * invrad + &
                                           asyice6(2,2,ib)) * invrad + &
                                           asyice6(2,3,ib)) * invrad + &
                                           asyice6(2,4,ib)) * invrad + &
                                           asyice6(2,5,ib)) * invrad + &
                                           asyice6(2,6,ib)) * invrad + &
                                           asyice6(2,7,ib)) * invrad + &
                                           asyice6(2,8,ib)
                   endif  ! end if of diaice .ge. 25.
                 enddo    ! end do of lwbands for cloud radiative coefficients
             else
                 print*,'ERROR in UMRad: current iceflag is not available.'
                 stop
             end if  ! iceflag=1

             do ib=1,nlwbands
                taucloud = iciwp(i,k) * extcoice(ib)  ! iwp units: g/m2
                ssacloud = ssacoice(ib)
                ! delta-scaling (Liou, 2002, 313p)
                asycloud = asycoice(ib)
                fpeak = asycloud * asycloud

                taucloud = (1._wp-ssacloud*fpeak) * taucloud   ! delta-scaling technique, ref: Joseph, Wiscombe and Weinman (1976, JAS)
                ssacloud = (1._wp-fpeak)*ssacloud / &
                                (1._wp-ssacloud*fpeak)
                asycloud = (asycloud-fpeak) / (1._wp-fpeak)
                if (flag_scat) then
                  tau(i,k,ib) = taucloud
                  tau_w(i,k,ib) = taucloud*ssacloud
                  tau_w_g(i,k,ib) = taucloud*ssacloud*asycloud
                else 
                  tau(i,k,ib) = iciwp(i,k) * extcoice(ib) * (1._wp-ssacoice(ib))
                  tau_w(i,k,ib) = 0._wp
                  tau_w_g(i,k,ib) = 0._wp
                end if
             enddo    ! end do of lwbands for cloud radiative coefficients

          else
            print*,'ice effective diameter is out of range in mo_ice_optics_mc6'
            stop
          endif ! end if of ice effective diameter range
        endif   ! end if of ice cloud layer
      enddo     ! end do of k
   enddo        ! end do of i

   return
end subroutine compute_all_from_mc6

end module mo_ice_optics_mc6
