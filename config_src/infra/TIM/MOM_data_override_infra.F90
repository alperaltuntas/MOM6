!> These interfaces allow for ocean or sea-ice variables to be replaced with data.
!!
!! TIM infrastructure: data override is NOT supported. In CESM coupled mode it
!! is dormant (impose_data_init is only reached when ALLOW_FLUX_ADJUSTMENTS or
!! LIQUID_RUNOFF_FROM_DATA is enabled, both default False), so these stubs
!! fail loudly if a configuration actually requests the capability rather
!! than silently skipping the override. Configurations that need data
!! override should use the FMS1/FMS2 infrastructure.
module MOM_data_override_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domain_infra,  only : MOM_domain_type, domain2d
use MOM_time_manager,  only : time_type
use MOM_error_infra,   only : MOM_err, FATAL

implicit none ; private

public :: impose_data_init, impose_data, impose_data_unset_domains

!> Potentially override the values of a field in the model with values from a dataset.
interface impose_data
  module procedure data_override_MD, data_override_2d
end interface

contains

!> Data override is unavailable under the TIM infrastructure; reaching this
!! call means a configuration enabled an override-dependent option.
subroutine impose_data_init(MOM_domain_in, Ocean_domain_in, Ice_domain_in)
  type (MOM_domain_type), intent(in), optional :: MOM_domain_in !< Ocean model domain
  type (domain2d),        intent(in), optional :: Ocean_domain_in !< Ocean model domain
  type (domain2d),        intent(in), optional :: Ice_domain_in !< Sea-ice model domain

  call MOM_err(FATAL, "MOM_data_override_infra: data override is not "//&
      "supported under the TIM infrastructure (enable it only with FMS1/FMS2).")
end subroutine impose_data_init

!> Data override is unavailable under the TIM infrastructure (fatal).
subroutine data_override_MD(domain, fieldname, data_2D, time, scale, override, is_ice)
  type(MOM_domain_type), intent(in)   :: domain   !< MOM domain from which to extract information
  character(len=*),     intent(in)    :: fieldname !< Name of the field to override
  real, dimension(:,:), intent(inout) :: data_2D  !< Data that may be modified by this call.
  type(time_type),      intent(in)    :: time     !< The model time, and the time for the data
  real,       optional, intent(in)    :: scale    !< A scaling factor for overridden fields
  logical,    optional, intent(out)   :: override !< True if the field has been overridden successfully
  logical,    optional, intent(in)    :: is_ice   !< If present and true, use the ice domain.

  if (present(override)) override = .false.
  call MOM_err(FATAL, "MOM_data_override_infra: data override of '"//&
      trim(fieldname)//"' is not supported under the TIM infrastructure.")
end subroutine data_override_MD

!> Data override is unavailable under the TIM infrastructure (fatal).
subroutine data_override_2d(gridname, fieldname, data_2D, time, override)
  character(len=3),     intent(in)    :: gridname !< Model component ('OCN' or 'ICE')
  character(len=*),     intent(in)    :: fieldname !< Name of the field to override
  real, dimension(:,:), intent(inout) :: data_2D  !< Data that may be modified by this call
  type(time_type),      intent(in)    :: time     !< The model time, and the time for the data
  logical,    optional, intent(out)   :: override !< True if the field has been overridden successfully

  if (present(override)) override = .false.
  call MOM_err(FATAL, "MOM_data_override_infra: data override of '"//&
      trim(fieldname)//"' is not supported under the TIM infrastructure.")
end subroutine data_override_2d

!> Unset domains that had previously been set for use by data_override.
!! A no-op under TIM (nothing can have been set).
subroutine impose_data_unset_domains(unset_Ocean, unset_Ice, must_be_set)
  logical, intent(in), optional :: unset_Ocean !< If present and true, unset the ocean domain for overrides
  logical, intent(in), optional :: unset_Ice   !< If present and true, unset the sea-ice domain for overrides
  logical, intent(in), optional :: must_be_set !< If present and true, it is a fatal error to unset
                                               !! a domain that is not set.
end subroutine impose_data_unset_domains

end module MOM_data_override_infra

!> \namespace MOM_data_override_infra
!!
!! Under the TIM infrastructure the FMS data_override capability is stubbed
!! out: dormant in CESM coupled mode, fatal if actually requested.
