! This code is part of Radiative Transfer for Energetics (RTE) and
!   RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! This module reads and writes to netCDF files (nominally for the Garand atmospheres) using
!   serial I/O. The files follow arbitrary conventions adopted by RTE+RRTMGP developers.
!   It may be useful as an example for other formats/conventions.
!
! Reading routines use "allocation on assignment," a feature of Fortran 2003 that may require
!   particular compilation flags.
!
module mo_garand_atmos_io
  !
  ! RTE+RRTMGP modules
  !
  use mo_rte_kind,           only: wp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_rrtmgp_util_reorder,only: reorder123x312
  use mo_optical_props,      only: ty_optical_props
  !
  ! NetCDF I/O routines, shared with other RTE+RRTMGP examples
  !
  use mo_simple_netcdf,      only: read_field, read_string, var_exists, get_dim_size, &
                                   write_field, create_dim, create_var

  use netcdf
  implicit none
  private

  public :: read_atmos, write_lw_fluxes, write_sw_fluxes
!! U-MICH add -->
  public :: write_lw_optics, write_sw_optics, write_cloud_props
!! <--

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read profiles for all columns  -- T, p, and gas concentrations
  !   Allocation occurs on assignments (says the F2003 standard)
  !
  subroutine read_atmos(fileName,                          &
                        p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry)
    character(len=*),   intent(in   ) :: fileName
    real(wp), dimension(:,:), allocatable,                 &
                        intent(inout) :: p_lay, t_lay, p_lev, t_lev, col_dry
    type(ty_gas_concs), intent(inout) :: gas_concs
    ! -------------------
    integer :: igas
    integer :: ncid, ncol, nlay, nlev
    integer, parameter :: ngas = 8
    character(len=3), dimension(ngas) &
                       :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']
    character(len=7) :: vmr_name
    ! -------------------

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_atmos: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    nlev = get_dim_size(ncid, 'lev')
    if(nlev /= nlay+1) call stop_on_err("read_atmos: nlev should be nlay+1")

    !
    ! These lines assume that compilers follow the Fortran 2003 standard for
    !   allocating on assignment. This may require explicit compiler support
    !   e.g. -assume realloc_lhs flag for Intel
    !
    p_lay = read_field(ncid, 'p_lay', ncol, nlay)
    t_lay = read_field(ncid, 't_lay', ncol, nlay)
    p_lev = read_field(ncid, 'p_lev', ncol, nlev)
    t_lev = read_field(ncid, 't_lev', ncol, nlev)

    call stop_on_err(gas_concs%init(gas_names))
    do igas = 1, ngas
      vmr_name = 'vmr_' // trim(gas_names(igas))
      if(.not. var_exists(ncid, trim(vmr_name))) &
        call stop_on_err("read_atmos: can't read concentration of " // trim(gas_names(igas)))
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(igas)), read_field(ncid, trim(vmr_name), ncol, nlay)))
    end do

    ! col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
    if(var_exists(ncid, 'col_dry')) col_dry = read_field(ncid, 'col_dry', ncol, nlay)

    ncid = nf90_close(ncid)

  end subroutine read_atmos
  !--------------------------------------------------------------------------------------------------------------------
  subroutine write_sw_fluxes(fileName, flux_up, flux_dn, flux_dir)
    character(len=*),         intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: flux_up, flux_dn, flux_dir
    ! -------------------
    integer :: ncid, ncol, nlev
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = size(flux_up, dim=1)
    nlev  = get_dim_size(ncid, 'lev')
    call create_dim(ncid, "col_flx", ncol)

    call create_var(ncid, "sw_flux_up",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "sw_flux_dn",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "sw_flux_dir", ["col_flx",  "lev    "], [ncol, nlev])

    call stop_on_err(write_field(ncid, "sw_flux_up",  flux_up ))
    call stop_on_err(write_field(ncid, "sw_flux_dn",  flux_dn ))
    call stop_on_err(write_field(ncid, "sw_flux_dir", flux_dir))

    ncid = nf90_close(ncid)
  end subroutine write_sw_fluxes
  !--------------------------------------------------------------------------------------------------------------------
  subroutine write_lw_fluxes(fileName, flux_up, flux_dn)
    character(len=*),         intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: flux_up, flux_dn
    ! -------------------
    integer :: ncid, ncol, nlev
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))

    ncol  = size(flux_up, dim=1)
    nlev  = get_dim_size(ncid, 'lev')
    call create_dim(ncid, "col_flx", ncol)

    call create_var(ncid, "lw_flux_up",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, "lw_flux_dn",  ["col_flx",  "lev    "], [ncol, nlev])

    call stop_on_err(write_field(ncid, "lw_flux_up",  flux_up ))
    call stop_on_err(write_field(ncid, "lw_flux_dn",  flux_dn ))

    ncid = nf90_close(ncid)
  end subroutine write_lw_fluxes
  !--------------------------------------------------------------------------------------------------------------------
!! U-MICH add -->
  subroutine write_lw_optics(fileName, tau, ssa, g)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp), optional, dimension(:,:,:), intent(in) ::  ssa, g
    real(wp), dimension(:,:), allocatable :: tau_tot, ssa_tot, g_tot, taussa_tot, taussag_tot
    ! -------------------
    integer :: ncid, ncol, nlay, nbnd
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))

    ncol  = size(tau, dim=1)
    nlay  = get_dim_size(ncid, 'lay')
    nbnd  = size(tau, dim=3)
    call create_dim(ncid, "col_flx", ncol)
    call create_dim(ncid, "bnd_lw" , nbnd)
  ! create profile variables
    call create_var(ncid, "lw_tau",  ["col_flx",  "lay    ", "bnd_lw "], [ncol, nlay, nbnd])
    if (present(ssa)) &
    call create_var(ncid, "lw_ssa",  ["col_flx",  "lay    ", "bnd_lw "], [ncol, nlay, nbnd])
    if (present(g)) &
    call create_var(ncid, "lw_g"  ,  ["col_flx",  "lay    ", "bnd_lw "], [ncol, nlay, nbnd])

  ! create accumulated variables
    call create_var(ncid, "lw_tau_tot",  ["col_flx", "bnd_lw "], [ncol, nbnd])
    if (present(ssa)) &
    call create_var(ncid, "lw_ssa_tot",  ["col_flx", "bnd_lw "], [ncol, nbnd])
    if (present(g)) &
    call create_var(ncid, "lw_g_tot",  ["col_flx", "bnd_lw "], [ncol, nbnd])

  ! calculate vertical accumulation 
    allocate(tau_tot(ncol,nbnd))
    tau_tot(:,:)=sum(tau(:,:,:),dim=2)

    if (present(ssa)) then
      allocate(taussa_tot(ncol,nbnd))
      allocate(ssa_tot(ncol,nbnd))
      ssa_tot(:,:)=0._wp
      taussa_tot(:,:)=sum(tau(:,:,:)*ssa(:,:,:),dim=2)
      where(tau_tot .gt. epsilon(0._wp)) ssa_tot=taussa_tot/tau_tot
      deallocate(taussa_tot)
    end if

    if (present(g)) then
      allocate(taussag_tot(ncol,nbnd))
      allocate(taussa_tot(ncol,nbnd))
      allocate(g_tot(ncol,nbnd))
      g_tot(:,:)=0._wp
      taussag_tot(:,:)=sum(tau(:,:,:)*ssa(:,:,:)*g(:,:,:),dim=2)
      taussa_tot(:,:)=sum(tau(:,:,:)*ssa(:,:,:),dim=2)
      where(taussa_tot .gt. epsilon(0._wp)) g_tot=taussag_tot/taussa_tot
      deallocate(taussag_tot)
      deallocate(taussa_tot)
    end if

  ! write data to file
    call stop_on_err(write_field(ncid, "lw_tau"    ,      tau ))
    call stop_on_err(write_field(ncid, "lw_tau_tot",  tau_tot ))
    if (present(ssa)) then
      call stop_on_err(write_field(ncid, "lw_ssa"    ,      ssa ))
      call stop_on_err(write_field(ncid, "lw_ssa_tot",  ssa_tot ))
    end if
    if (present(g)) then
      call stop_on_err(write_field(ncid, "lw_g"      ,      g   ))
      call stop_on_err(write_field(ncid, "lw_g_tot"  ,  g_tot   ))
    end if
    ncid = nf90_close(ncid)
  end subroutine write_lw_optics
!! <--
  !--------------------------------------------------------------------------------------------------------------------
!! U-MICH add -->
  subroutine write_sw_optics(fileName, tau, ssa, g)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp),optional, dimension(:,:,:), intent(in) :: ssa, g
    ! -------------------
    integer :: ncid, ncol, nlay, nbnd
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))

    ncol  = size(tau, dim=1)
    nlay  = get_dim_size(ncid, 'lay')
    nbnd  = size(tau, dim=3)
    call create_dim(ncid, "col_flx", ncol)
    call create_dim(ncid, "bnd_sw" , nbnd)

    call create_var(ncid, "sw_tau",  ["col_flx",  "lay    ", "bnd_sw "], [ncol, nlay, nbnd])
    if (present(ssa)) &
    call create_var(ncid, "sw_ssa",  ["col_flx",  "lay    ", "bnd_sw "], [ncol, nlay, nbnd])
    if (present(g)) &
    call create_var(ncid, "sw_g"  ,  ["col_flx",  "lay    ", "bnd_sw "], [ncol, nlay, nbnd])

    call stop_on_err(write_field(ncid, "sw_tau",  tau ))
    if (present(ssa)) &
    call stop_on_err(write_field(ncid, "sw_ssa",  ssa ))
    if (present(g)) &
    call stop_on_err(write_field(ncid, "sw_g"  ,  g   ))

    ncid = nf90_close(ncid)
  end subroutine write_sw_optics
!! <--
  !--------------------------------------------------------------------------------------------------------------------
!! U-MICH add -->
  subroutine write_cloud_props(fileName, lwp, iwp, rel, rei)
    character(len=*),         intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: lwp, iwp, rel, rei
    ! -------------------
    integer :: ncid, ncol, nlay
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))

    ncol  = size(lwp, dim=1)
    nlay  = get_dim_size(ncid, 'lay')
    call create_dim(ncid, "col_flx", ncol)

    call create_var(ncid, "lwp",  ["col_flx",  "lay    "], [ncol, nlay])
    call create_var(ncid, "iwp",  ["col_flx",  "lay    "], [ncol, nlay])
    call create_var(ncid, "rel",  ["col_flx",  "lay    "], [ncol, nlay])
    call create_var(ncid, "rei",  ["col_flx",  "lay    "], [ncol, nlay])

    call stop_on_err(write_field(ncid, "lwp",  lwp ))
    call stop_on_err(write_field(ncid, "iwp",  iwp ))
    call stop_on_err(write_field(ncid, "rel",  rel ))
    call stop_on_err(write_field(ncid, "rei",  rei ))

    ncid = nf90_close(ncid)
  end subroutine write_cloud_props
!! <--
  !--------------------------------------------------------------------------------------------------------------------
  subroutine stop_on_err(msg)
    !
    ! Print error meiwpge and stop
    !
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg
    if(len_trim(msg) > 0) then
      write(error_unit,*) trim(msg)
      stop
    end if
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------
end module mo_garand_atmos_io
