!> A wrapper for the diag_manager routines. This module is MOM6's only
!! import point for diagnostic-manager infrastructure. It dispatches at
!! runtime between the FMS diag_manager (default) and the TIM prototype
!! diagnostics manager (config key tim.diag / env TIM_DIAG) — a WHOLE-RUN
!! switch: exactly one backend owns every history file in a run.
module MOM_diag_manager_infra

! This file is part of MOM6. See LICENSE.md for the license.

use, intrinsic :: iso_fortran_env, only : real64
use, intrinsic :: iso_c_binding,   only : c_double, c_char, c_null_char, &
                                          c_long_long, c_ptr, c_null_ptr, c_loc
! Position constants: the same mpp_domains values MOM_domain_infra re-exports
! as EAST_FACE/NORTH_FACE (FMS diag_axis used the identical constants).
use MOM_domain_infra, only : EAST => EAST_FACE, NORTH => NORTH_FACE
use MOM_time_manager, only : time_type, get_time, get_calendar_type
use MOM_domain_infra, only : MOM_domain_type
use MOM_error_infra,  only : MOM_err, FATAL, WARNING
use MOM_io_infra,     only : tim_get_domain_handle
use tim_diag_interface, only : tim_diag_init, tim_diag_active, tim_diag_axis_init
use tim_diag_interface, only : tim_diag_axis_name, tim_diag_register_field
use tim_diag_interface, only : tim_diag_field_id, tim_diag_attr_text
use tim_diag_interface, only : tim_diag_attr_ints, tim_diag_attr_reals
use tim_diag_interface, only : tim_diag_post, tim_diag_send_complete
use tim_diag_interface, only : tim_diag_set_time_end, tim_diag_end
use tim_diag_interface, only : tim_diag_save_state, tim_diag_restore_state
use tim_diag_interface, only : cstr

implicit none ; private

!> Special axis id returned for scalar diagnostics (TIM null axis).
integer, parameter :: null_axis_id = 0
!> Returned by the register calls when the diag_table requests no output.
integer, parameter :: DIAG_FIELD_NOT_FOUND = -1

!> transmit data for diagnostic output
interface register_diag_field_infra
  module procedure register_diag_field_infra_scalar
  module procedure register_diag_field_infra_array
end interface register_diag_field_infra

!> transmit data for diagnostic output
interface send_data_infra
  module procedure send_data_infra_0d, send_data_infra_1d
  module procedure send_data_infra_2d, send_data_infra_3d
#ifdef OVERLOAD_R8
  module procedure send_data_infra_2d_r8, send_data_infra_3d_r8
#endif
end interface send_data_infra

!> Add an attribute to a diagnostic field
interface MOM_diag_field_add_attribute
  module procedure MOM_diag_field_add_attribute_scalar_r
  module procedure MOM_diag_field_add_attribute_scalar_i
  module procedure MOM_diag_field_add_attribute_scalar_c
  module procedure MOM_diag_field_add_attribute_r1d
  module procedure MOM_diag_field_add_attribute_i1d
end interface MOM_diag_field_add_attribute


! Public interfaces
public MOM_diag_axis_init
public get_MOM_diag_axis_name
public MOM_diag_manager_init
public MOM_diag_manager_end
public send_data_infra
public diag_send_complete_infra
public diag_manager_set_time_end_infra
public MOM_diag_field_add_attribute
public register_diag_field_infra
public register_static_field_infra
public get_MOM_diag_field_id
! Restart-spanning averaging (TIM capability; no-ops under FMS)
public MOM_diag_save_state, MOM_diag_restore_state
! Public data
public null_axis_id
public DIAG_FIELD_NOT_FOUND
public EAST, NORTH

contains


!> Warn (once) that coarsened (downsampled) diag axes are not supported.
subroutine coarse_axis_warning()
  logical, save :: warned = .false.
  if (.not. warned) then
    call MOM_err(WARNING, "TIM diag: coarsened (downsampled) diagnostics are "//&
                          "not supported by the prototype; such fields will not be output.")
    warned = .true.
  endif
end subroutine coarse_axis_warning

!> Warn (once) about send_data features the TIM prototype ignores.
subroutine tim_unsupported_mask_warning()
  logical, save :: warned = .false.
  if (.not. warned) then
    call MOM_err(WARNING, "TIM diag: logical mask argument to send_data is "//&
                          "ignored (MOM6 uses rmask; logical masks are not bridged)")
    warned = .true.
  endif
end subroutine tim_unsupported_mask_warning

!> Initialize a diagnostic axis
integer function MOM_diag_axis_init(name, data, units, cart_name, long_name, MOM_domain, position, &
          & direction, edges, set_name, coarsen, null_axis)
  character(len=*),   intent(in) :: name      !< The name of this axis
  real, dimension(:), intent(in) :: data      !< The array of coordinate values
  character(len=*),   intent(in) :: units     !< The units for the axis data
  character(len=*),   intent(in) :: cart_name !< Cartesian axis ("X", "Y", "Z", "T", or "N" for none)
  character(len=*), &
            optional, intent(in) :: long_name !< The long name of this axis
  type(MOM_domain_type), &
            optional, intent(in) :: MOM_Domain !< A MOM_Domain that describes the decomposition
  integer,  optional, intent(in) :: position  !< This indicates the relative position of this
                                              !! axis.  The default is CENTER, but EAST and NORTH
                                              !! are common options.
  integer,  optional, intent(in) :: direction !< This indicates the direction along which this
                                              !! axis increases: 1 for upward, -1 for downward, or
                                              !! 0 for non-vertical axes (the default)
  integer,  optional, intent(in) :: edges     !< The axis_id of the complementary axis that
                                              !! describes the edges of this axis
  character(len=*), &
            optional, intent(in) :: set_name  !< A name to use for this set of axes.
  integer,  optional, intent(in) :: coarsen   !< An optional degree of coarsening for the grid, 1
                                              !! by default.
  logical,  optional, intent(in) :: null_axis !< If present and true, return the special null axis
                                              !! id for use with scalars.

  integer :: coarsening ! The degree of grid coarsening
  integer :: dom_handle, staggered, dir, edge_id
  character(len=256) :: lname, sname
  real(kind=c_double) :: cdata(size(data))

  if (present(null_axis)) then ; if (null_axis) then
    MOM_diag_axis_init = 0  ! the TIM null axis id
    return
  endif ; endif
  if (present(coarsen)) then ; if (coarsen /= 1) then
    ! MOM's diag mediator registers downsampled axes unconditionally; only
    ! fields actually requesting downsampled output ever use them. Hand
    ! back a sentinel id so such fields fail registration cleanly.
    call coarse_axis_warning()
    MOM_diag_axis_init = -2
    return
  endif ; endif
  dom_handle = -1
  if (present(MOM_domain)) dom_handle = tim_get_domain_handle(MOM_domain)
  staggered = 0
  if (present(position)) then
    if (position == EAST .or. position == NORTH) staggered = 1
  endif
  dir = 0 ; if (present(direction)) dir = direction
  edge_id = 0 ; if (present(edges)) edge_id = edges
  lname = "" ; if (present(long_name)) lname = long_name
  sname = "" ; if (present(set_name)) sname = set_name
  cdata(:) = real(data(:), kind=c_double)
  MOM_diag_axis_init = tim_diag_axis_init(cstr(name), cdata, size(data), &
      cstr(units), cstr(cart_name), cstr(lname), dom_handle, staggered, &
      dir, edge_id, cstr(sname))

end function MOM_diag_axis_init

!> Returns the short name of the axis
subroutine get_MOM_diag_axis_name(id, name)
  integer,          intent(in)  :: id   !< The axis numeric id
  character(len=*), intent(out) :: name !< The short name of the axis

  character(kind=c_char) :: cbuf(256)
  integer :: k

  call tim_diag_axis_name(id, cbuf, size(cbuf))
  name = ""
  do k = 1, min(size(cbuf), len(name))
    if (cbuf(k) == c_null_char) exit
    name(k:k) = cbuf(k)
  enddo

end subroutine get_MOM_diag_axis_name

!> Return a unique numeric ID field a module/field name combination.
integer function get_MOM_diag_field_id(module_name, field_name)
  character(len=*), intent(in) :: module_name !< A module name string to query.
  character(len=*), intent(in) :: field_name  !< A field name string to query.

  get_MOM_diag_field_id = tim_diag_field_id(cstr(module_name), cstr(field_name))

end function get_MOM_diag_field_id

!> Initializes the diagnostic manager
subroutine MOM_diag_manager_init(diag_model_subset, time_init, err_msg)
  integer,               optional, intent(in) :: diag_model_subset !< An optional diagnostic subset
  integer, dimension(6), optional, intent(in) :: time_init !< An optional reference time for diagnostics
                                                           !! The default uses the value contained in the
                                                           !! diag_table. Format is Y-M-D-H-M-S
  character(len=*),     optional, intent(out) :: err_msg   !< Error message.

  integer :: y, mo, d, h, mi, s, rc

  y = -1 ; mo = 1 ; d = 1 ; h = 0 ; mi = 0 ; s = 0
  if (present(time_init)) then
    y = time_init(1) ; mo = time_init(2) ; d = time_init(3)
    h = time_init(4) ; mi = time_init(5) ; s = time_init(6)
  endif
  rc = tim_diag_init(get_calendar_type(), y, mo, d, h, mi, s)
  if (rc /= 0) call MOM_err(FATAL, "MOM_diag_manager_init: tim_diag_init failed")
  if (present(err_msg)) err_msg = ""

end subroutine MOM_diag_manager_init

!> Close the diagnostic manager
subroutine MOM_diag_manager_end(time)
  type(time_type), intent(in) :: time !< Model time at call to close.

  integer :: days, secs

  call get_time(time, secs, days)
  if (tim_diag_end(days, secs) /= 0) &
    call MOM_err(WARNING, "MOM_diag_manager_end: tim_diag_end reported errors")

end subroutine MOM_diag_manager_end

!> Save every in-progress diagnostic averaging window to a restart file
!! (TIM restart-spanning averaging; a no-op under the FMS backend, which
!! cannot carry partial windows across restarts).
subroutine MOM_diag_save_state(filename)
  character(len=*), intent(in) :: filename !< Path of the diag state file

  if (tim_diag_save_state(cstr(filename)) /= 0) &
    call MOM_err(WARNING, "MOM_diag_save_state: save failed for "//trim(filename))
end subroutine MOM_diag_save_state

!> Restore in-progress diagnostic averaging windows from a state file.
!! Returns true when state was restored; false = cold start (no file) or FMS.
logical function MOM_diag_restore_state(filename)
  character(len=*), intent(in) :: filename !< Path of the diag state file

  MOM_diag_restore_state = (tim_diag_restore_state(cstr(filename)) == 0)
end function MOM_diag_restore_state

!> Register a MOM diagnostic field for scalars
integer function register_diag_field_infra_scalar(module_name, field_name, init_time, &
                        long_name, units, missing_value, range, standard_name, do_not_log, &
                        err_msg, area, volume)
  character(len=*),              intent(in) :: module_name !< The name of the associated module
  character(len=*),              intent(in) :: field_name !< The name of the field
  type(time_type),     optional, intent(in) :: init_time !< The registration time
  character(len=*),    optional, intent(in) :: long_name !< A long name for the field
  character(len=*),    optional, intent(in) :: units     !< Field units
  character(len=*),    optional, intent(in) :: standard_name !< A standard name for the field
  real,                optional, intent(in) :: missing_value !< Missing value attribute
  real,  dimension(2), optional, intent(in) :: range     !< A valid range of the field
  logical,             optional, intent(in) :: do_not_log !< if TRUE, field information is not logged
  character(len=*),    optional, intent(out):: err_msg   !< An error message to return
  integer,             optional, intent(in) :: area      !< Diagnostic ID of the field containing the area attribute
  integer,             optional, intent(in) :: volume    !< Diagnostic ID of the field containing the volume attribute

  integer :: axes(1)

  axes(1) = 0
  register_diag_field_infra_scalar = tim_register_bridge(module_name, field_name, &
      axes, 0, init_time, long_name, units, standard_name, missing_value=missing_value, &
      range=range, is_static=.false., area=area, volume=volume)
  if (present(err_msg)) err_msg = ""

end function register_diag_field_infra_scalar

!> Register a MOM diagnostic field for scalars
integer function register_diag_field_infra_array(module_name, field_name, axes, init_time, &
                        long_name, units, missing_value, range, mask_variant, standard_name, verbose, &
                        do_not_log, err_msg, interp_method, tile_count, area, volume)
  character(len=*),             intent(in) :: module_name !< The name of the associated module
  character(len=*),             intent(in) :: field_name !< The name of the field
  integer, dimension(:),        intent(in) :: axes      !< Diagnostic IDs of axis attributes for the field
  type(time_type),    optional, intent(in) :: init_time !< The registration time
  character(len=*),   optional, intent(in) :: long_name !< A long name for the field
  character(len=*),   optional, intent(in) :: units     !< Units of the field
  real,               optional, intent(in) :: missing_value !< Missing value attribute
  real, dimension(2), optional, intent(in) :: range     !< A valid range of the field
  logical,            optional, intent(in) :: mask_variant !< If true, the field mask is varying in time
  character(len=*),   optional, intent(in) :: standard_name !< A standard name for the field
  logical,            optional, intent(in) :: verbose    !< If true, provide additional log information
  logical,            optional, intent(in) :: do_not_log !< if TRUE, field information is not logged
  character(len=*),   optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,            optional, intent(in) :: tile_count !< The tile number for the current PE
  character(len=*),   optional, intent(out):: err_msg   !< An error message to return
  integer,            optional, intent(in) :: area      !< Diagnostic ID of the field containing the area attribute
  integer,            optional, intent(in) :: volume    !< Diagnostic ID of the field containing the volume attribute

  register_diag_field_infra_array = tim_register_bridge(module_name, field_name, &
      axes, size(axes), init_time, long_name, units, standard_name, &
      interp_method=interp_method, missing_value=missing_value, range=range, &
      mask_variant=mask_variant, is_static=.false., area=area, volume=volume)
  if (present(err_msg)) err_msg = ""

end function register_diag_field_infra_array


integer function register_static_field_infra(module_name, field_name, axes, long_name, units, &
                        missing_value, range, mask_variant, standard_name, do_not_log, interp_method, &
                        tile_count, area, volume)
  character(len=*),             intent(in) :: module_name !< The name of the associated module
  character(len=*),             intent(in) :: field_name !< The name of the field
  integer, dimension(:),        intent(in) :: axes      !< Diagnostic IDs of axis attributes for the field
  character(len=*),   optional, intent(in) :: long_name !< A long name for the field
  character(len=*),   optional, intent(in) :: units     !< Units of the field
  real,               optional, intent(in) :: missing_value !< Missing value attribute
  real, dimension(2), optional, intent(in) :: range     !< A valid range of the field
  logical,            optional, intent(in) :: mask_variant !< If true, the field mask is varying in time
  character(len=*),   optional, intent(in) :: standard_name !< A standard name for the field
  logical,            optional, intent(in) :: do_not_log !< if TRUE, field information is not logged
  character(len=*),   optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,            optional, intent(in) :: tile_count !< The tile number for the current PE
  integer,            optional, intent(in) :: area      !< Diagnostic ID of the field containing the area attribute
  integer,            optional, intent(in) :: volume    !< Diagnostic ID of the field containing the volume attribute

  register_static_field_infra = tim_register_bridge(module_name, field_name, &
      axes, size(axes), long_name=long_name, units=units, standard_name=standard_name, &
      interp_method=interp_method, missing_value=missing_value, range=range, &
      mask_variant=mask_variant, is_static=.true., area=area, volume=volume)

end function register_static_field_infra

!> Shared TIM registration bridge: optional-argument resolution in one place.
integer function tim_register_bridge(module_name, field_name, axes, naxes, init_time, &
                        long_name, units, standard_name, interp_method, missing_value, &
                        range, mask_variant, is_static, area, volume)
  character(len=*),             intent(in) :: module_name !< The name of the associated module
  character(len=*),             intent(in) :: field_name  !< The name of the field
  integer, dimension(:),        intent(in) :: axes        !< TIM axis ids (0 = null axis)
  integer,                      intent(in) :: naxes       !< Number of valid axes (0 = scalar)
  type(time_type),    optional, intent(in) :: init_time   !< The registration time
  character(len=*),   optional, intent(in) :: long_name   !< A long name for the field
  character(len=*),   optional, intent(in) :: units       !< Units of the field
  character(len=*),   optional, intent(in) :: standard_name !< A standard name for the field
  character(len=*),   optional, intent(in) :: interp_method !< Scalar-interpolation hint
  real,               optional, intent(in) :: missing_value !< Missing value attribute
  real, dimension(2), optional, intent(in) :: range       !< A valid range of the field
  logical,            optional, intent(in) :: mask_variant !< Time-varying mask flag
  logical,                      intent(in) :: is_static   !< Static (no time axis) field
  integer,            optional, intent(in) :: area        !< Diag ID of the cell_measures area field
  integer,            optional, intent(in) :: volume      !< Diag ID of the cell_measures volume field

  character(len=256) :: lname, uname, sname, imeth
  integer :: init_days, init_secs, has_miss, has_range, maskv, statf, aid, vid
  real(kind=c_double) :: missv, rlo, rhi

  lname = "" ; if (present(long_name)) lname = long_name
  uname = "" ; if (present(units)) uname = units
  sname = "" ; if (present(standard_name)) sname = standard_name
  imeth = "" ; if (present(interp_method)) imeth = interp_method
  init_days = -1 ; init_secs = 0
  if (present(init_time)) call get_time(init_time, init_secs, init_days)
  has_miss = 0 ; missv = 0.0_c_double
  if (present(missing_value)) then
    has_miss = 1 ; missv = real(missing_value, kind=c_double)
  endif
  has_range = 0 ; rlo = 0.0_c_double ; rhi = 0.0_c_double
  if (present(range)) then
    has_range = 1 ; rlo = real(range(1), kind=c_double) ; rhi = real(range(2), kind=c_double)
  endif
  maskv = 0 ; if (present(mask_variant)) then ; if (mask_variant) maskv = 1 ; endif
  statf = 0 ; if (is_static) statf = 1
  aid = -1 ; if (present(area)) aid = area
  vid = -1 ; if (present(volume)) vid = volume

  tim_register_bridge = tim_diag_register_field(cstr(module_name), cstr(field_name), &
      axes, naxes, init_days, init_secs, cstr(lname), cstr(uname), cstr(sname), &
      cstr(imeth), has_miss, missv, has_range, rlo, rhi, maskv, statf, aid, vid)

end function tim_register_bridge

!> Returns true if the argument data are successfully passed to a diagnostic manager
!! with the indicated unique reference id, false otherwise.
logical function send_data_infra_0d(diag_field_id, field, time, err_msg)
  integer,                    intent(in)  :: diag_field_id !< The diagnostic manager identifier for this field
  real,                       intent(in)  :: field   !< The value being recorded
  TYPE(time_type),  optional, intent(in)  :: time    !< The time for the current record
  CHARACTER(len=*), optional, intent(out) :: err_msg !< An optional error message

  integer :: days, secs
  real(kind=c_double) :: dbuf(1)

  days = -1 ; secs = 0
  if (present(time)) call get_time(time, secs, days)
  dbuf(1) = real(field, kind=c_double)
  send_data_infra_0d = (tim_diag_post(diag_field_id, days, secs, dbuf, &
      1_c_long_long, c_null_ptr, 1.0_c_double) /= 0)
  if (present(err_msg)) err_msg = ""

end function send_data_infra_0d

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_1d(diag_field_id, field, is_in, ie_in, time, mask, rmask, weight, err_msg)
  integer,                         intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real, dimension(:),              intent(in) :: field !< A 1-d array of values being recorded
  integer,               optional, intent(in) :: is_in !< The starting index for the data being recorded
  integer,               optional, intent(in) :: ie_in !< The end index for the data being recorded
  type(time_type),       optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:), optional, intent(in) :: mask  !< An optional rank 1 logical mask
  real, dimension(:),    optional, intent(in) :: rmask !< An optional rank 1 mask array
  real,                  optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                       !! record if there is averaging in time
  character(len=*),      optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                       !! returning to the calling routine

  integer :: isv, iev, days, secs, i, n
  real(kind=c_double), allocatable :: dbuf(:)
  real(kind=c_double), allocatable, target :: rbuf(:)
  type(c_ptr) :: rptr
  real(kind=c_double) :: wt

  if (present(mask)) call tim_unsupported_mask_warning()
  isv = 1 ; if (present(is_in)) isv = is_in
  iev = isv + size(field) - 1 ; if (present(ie_in)) iev = ie_in
  allocate(dbuf(iev-isv+1))
  n = 0
  do i=isv,iev ; n = n+1 ; dbuf(n) = real(field(i), kind=c_double) ; enddo
  rptr = c_null_ptr
  if (present(rmask)) then
    allocate(rbuf(iev-isv+1))
    n = 0
    do i=isv,iev ; n = n+1 ; rbuf(n) = real(rmask(i), kind=c_double) ; enddo
    rptr = c_loc(rbuf(1))
  endif
  days = -1 ; secs = 0
  if (present(time)) call get_time(time, secs, days)
  wt = 1.0_c_double ; if (present(weight)) wt = real(weight, kind=c_double)
  send_data_infra_1d = (tim_diag_post(diag_field_id, days, secs, dbuf, &
      int(iev-isv+1, kind=c_long_long), rptr, wt) /= 0)
  if (present(err_msg)) err_msg = ""

end function send_data_infra_1d

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_2d(diag_field_id, field, is_in, ie_in, js_in, je_in, &
                                    time, mask, rmask, weight, err_msg)
  integer,                           intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real, dimension(:,:),              intent(in) :: field !< A 2-d array of values being recorded
  integer,                 optional, intent(in) :: is_in !< The starting i-index for the data being recorded
  integer,                 optional, intent(in) :: ie_in !< The end i-index for the data being recorded
  integer,                 optional, intent(in) :: js_in !< The starting j-index for the data being recorded
  integer,                 optional, intent(in) :: je_in !< The end j-index for the data being recorded
  type(time_type),         optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:,:), optional, intent(in) :: mask  !< An optional 2-d logical mask
  real, dimension(:,:),    optional, intent(in) :: rmask !< An optional 2-d mask array
  real,                    optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                         !! record if there is averaging in time
  character(len=*),        optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                         !! returning to the calling routine

  if (present(mask)) call tim_unsupported_mask_warning()
  if (present(rmask)) then
    send_data_infra_2d = tim_post_bridge(diag_field_id, &
        reshape(field, (/size(field,1), size(field,2), 1/)), &
        is_in, ie_in, js_in, je_in, 1, 1, time, &
        rmask3d=reshape(rmask, (/size(rmask,1), size(rmask,2), 1/)), &
        weight=weight)
  else
    send_data_infra_2d = tim_post_bridge(diag_field_id, &
        reshape(field, (/size(field,1), size(field,2), 1/)), &
        is_in, ie_in, js_in, je_in, 1, 1, time, weight=weight)
  endif
  if (present(err_msg)) err_msg = ""

end function send_data_infra_2d

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_3d(diag_field_id, field, is_in, ie_in, js_in, je_in, ks_in, ke_in, &
                                    time, mask, rmask, weight, err_msg)
  integer,                             intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real, dimension(:,:,:),              intent(in) :: field !< A rank 1 array of floating point values being recorded
  integer,                   optional, intent(in) :: is_in !< The starting i-index for the data being recorded
  integer,                   optional, intent(in) :: ie_in !< The end i-index for the data being recorded
  integer,                   optional, intent(in) :: js_in !< The starting j-index for the data being recorded
  integer,                   optional, intent(in) :: je_in !< The end j-index for the data being recorded
  integer,                   optional, intent(in) :: ks_in !< The starting k-index for the data being recorded
  integer,                   optional, intent(in) :: ke_in !< The end k-index for the data being recorded
  type(time_type),           optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:,:,:), optional, intent(in) :: mask  !< An optional 3-d logical mask
  real, dimension(:,:,:),    optional, intent(in) :: rmask !< An optional 3-d mask array
  real,                      optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                           !! record if there is averaging in time
  character(len=*),          optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                           !! returning to the calling routine

  if (present(mask)) call tim_unsupported_mask_warning()
  send_data_infra_3d = tim_post_bridge(diag_field_id, field, &
      is_in, ie_in, js_in, je_in, ks_in, ke_in, time, rmask3d=rmask, &
      weight=weight)
  if (present(err_msg)) err_msg = ""

end function send_data_infra_3d

!> Shared TIM post bridge: FMS is/ie window defaulting, contiguous r8 copy,
!! one tim_diag_post call. The field is the (possibly larger) local array;
!! the is/ie...ks/ke window selects the points that go to the manager.
logical function tim_post_bridge(diag_field_id, field, is_in, ie_in, js_in, je_in, &
                                 ks_in, ke_in, time, rmask3d, weight)
  integer,                intent(in) :: diag_field_id !< TIM diag field id
  real, dimension(:,:,:), intent(in) :: field !< Local data (2-d passed as single layer)
  integer,      optional, intent(in) :: is_in !< Start i-index of the window
  integer,      optional, intent(in) :: ie_in !< End i-index of the window
  integer,      optional, intent(in) :: js_in !< Start j-index of the window
  integer,      optional, intent(in) :: je_in !< End j-index of the window
  integer,      optional, intent(in) :: ks_in !< Start k-index of the window
  integer,      optional, intent(in) :: ke_in !< End k-index of the window
  type(time_type), optional, intent(in) :: time !< The time for the current record
  real, dimension(:,:,:), optional, intent(in) :: rmask3d !< Real mask, field-shaped
  real,         optional, intent(in) :: weight !< Weight for time averaging

  integer :: isv, iev, jsv, jev, ksv, kev, days, secs, i, j, k, n
  real(kind=c_double), allocatable :: dbuf(:)
  real(kind=c_double), allocatable, target :: rbuf(:)
  type(c_ptr) :: rptr
  real(kind=c_double) :: wt

  isv = 1 ; if (present(is_in)) isv = is_in
  iev = isv + size(field,1) - 1 ; if (present(ie_in)) iev = ie_in
  jsv = 1 ; if (present(js_in)) jsv = js_in
  jev = jsv + size(field,2) - 1 ; if (present(je_in)) jev = je_in
  ksv = 1 ; if (present(ks_in)) ksv = ks_in
  kev = ksv + size(field,3) - 1 ; if (present(ke_in)) kev = ke_in

  allocate(dbuf((iev-isv+1)*(jev-jsv+1)*(kev-ksv+1)))
  n = 0
  do k=ksv,kev ; do j=jsv,jev ; do i=isv,iev
    n = n+1 ; dbuf(n) = real(field(i,j,k), kind=c_double)
  enddo ; enddo ; enddo

  rptr = c_null_ptr
  if (present(rmask3d)) then
    allocate(rbuf(size(dbuf)))
    n = 0
    do k=ksv,kev ; do j=jsv,jev ; do i=isv,iev
      n = n+1 ; rbuf(n) = real(rmask3d(i,j,k), kind=c_double)
    enddo ; enddo ; enddo
    rptr = c_loc(rbuf(1))
  endif

  days = -1 ; secs = 0
  if (present(time)) call get_time(time, secs, days)
  wt = 1.0_c_double ; if (present(weight)) wt = real(weight, kind=c_double)

  tim_post_bridge = (tim_diag_post(diag_field_id, days, secs, dbuf, &
      int(size(dbuf), kind=c_long_long), rptr, wt) /= 0)

end function tim_post_bridge


#ifdef OVERLOAD_R8
!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_2d_r8(diag_field_id, field, is_in, ie_in, js_in, je_in, &
                                       time, mask, rmask, weight, err_msg)
  integer,                           intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real(kind=real64), dimension(:,:), intent(in) :: field !< A 2-d array of values being recorded
  integer,                 optional, intent(in) :: is_in !< The starting i-index for the data being recorded
  integer,                 optional, intent(in) :: ie_in !< The end i-index for the data being recorded
  integer,                 optional, intent(in) :: js_in !< The starting j-index for the data being recorded
  integer,                 optional, intent(in) :: je_in !< The end j-index for the data being recorded
  type(time_type),         optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:,:), optional, intent(in) :: mask  !< An optional 2-d logical mask
  real, dimension(:,:),    optional, intent(in) :: rmask !< An optional 2-d mask array
  real,                    optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                         !! record if there is averaging in time
  character(len=*),        optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                         !! returning to the calling routine

  if (present(mask)) call tim_unsupported_mask_warning()
  if (present(rmask)) then
    send_data_infra_2d_r8 = tim_post_bridge(diag_field_id, &
        real(reshape(field, (/size(field,1), size(field,2), 1/))), &
        is_in, ie_in, js_in, je_in, 1, 1, time, &
        rmask3d=reshape(rmask, (/size(rmask,1), size(rmask,2), 1/)), &
        weight=weight)
  else
    send_data_infra_2d_r8 = tim_post_bridge(diag_field_id, &
        real(reshape(field, (/size(field,1), size(field,2), 1/))), &
        is_in, ie_in, js_in, je_in, 1, 1, time, weight=weight)
  endif
  if (present(err_msg)) err_msg = ""

end function send_data_infra_2d_r8

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_3d_r8(diag_field_id, field, is_in, ie_in, js_in, je_in, ks_in, ke_in, &
                                    time, mask, rmask, weight, err_msg)
  integer,                             intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real(kind=real64), dimension(:,:,:), intent(in) :: field !< A rank 1 array of floating point values being recorded
  integer,                   optional, intent(in) :: is_in !< The starting i-index for the data being recorded
  integer,                   optional, intent(in) :: ie_in !< The end i-index for the data being recorded
  integer,                   optional, intent(in) :: js_in !< The starting j-index for the data being recorded
  integer,                   optional, intent(in) :: je_in !< The end j-index for the data being recorded
  integer,                   optional, intent(in) :: ks_in !< The starting k-index for the data being recorded
  integer,                   optional, intent(in) :: ke_in !< The end k-index for the data being recorded
  type(time_type),           optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:,:,:), optional, intent(in) :: mask  !< An optional 3-d logical mask
  real, dimension(:,:,:),    optional, intent(in) :: rmask !< An optional 3-d mask array
  real,                      optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                           !! record if there is averaging in time
  character(len=*),          optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                           !! returning to the calling routine

  if (present(mask)) call tim_unsupported_mask_warning()
  send_data_infra_3d_r8 = tim_post_bridge(diag_field_id, real(field), &
      is_in, ie_in, js_in, je_in, ks_in, ke_in, time, rmask3d=rmask, &
      weight=weight)
  if (present(err_msg)) err_msg = ""

end function send_data_infra_3d_r8
#endif

!> Add a real scalar attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_scalar_r(diag_field_id, att_name, att_value)
  integer,          intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*), intent(in) :: att_name  !< The name of the attribute
  real,             intent(in) :: att_value !< A real scalar value

  real(kind=c_double) :: v(1)

  v(1) = real(att_value, kind=c_double)
  call tim_diag_attr_reals(diag_field_id, cstr(att_name), v, 1)

end subroutine MOM_diag_field_add_attribute_scalar_r

!> Add an integer attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_scalar_i(diag_field_id, att_name, att_value)
  integer,          intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*), intent(in) :: att_name  !< The name of the attribute
  integer,          intent(in) :: att_value !< An integer scalar value

  integer :: v(1)

  v(1) = att_value
  call tim_diag_attr_ints(diag_field_id, cstr(att_name), v, 1)

end subroutine MOM_diag_field_add_attribute_scalar_i

!> Add a character string attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_scalar_c(diag_field_id, att_name, att_value)
  integer,          intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*), intent(in) :: att_name  !< The name of the attribute
  character(len=*), intent(in) :: att_value !< A character string value

  call tim_diag_attr_text(diag_field_id, cstr(att_name), cstr(att_value))

end subroutine MOM_diag_field_add_attribute_scalar_c

!> Add a real list of attributes attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_r1d(diag_field_id, att_name, att_value)
  integer,            intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*),   intent(in) :: att_name  !< The name of the attribute
  real, dimension(:), intent(in) :: att_value !< An array of real values

  real(kind=c_double) :: v(size(att_value))

  v(:) = real(att_value(:), kind=c_double)
  call tim_diag_attr_reals(diag_field_id, cstr(att_name), v, size(v))

end subroutine MOM_diag_field_add_attribute_r1d

!> Add a integer list of attributes attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_i1d(diag_field_id, att_name, att_value)
  integer,               intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*),      intent(in) :: att_name  !< The name of the attribute
  integer, dimension(:), intent(in) :: att_value !< An array of integer values

  call tim_diag_attr_ints(diag_field_id, cstr(att_name), att_value, size(att_value))

end subroutine MOM_diag_field_add_attribute_i1d

!> Finishes the diag manager reduction methods as needed for the time_step
subroutine diag_send_complete_infra ()
  !! The time_step in the diag_send_complete call is a dummy argument, needed for backwards compatibility
  !! It won't be used at all when diag_manager_nml::use_modern_diag=.true.
  !! It won't have any impact when diag_manager_nml::use_modern_diag=.false.
  call tim_diag_send_complete()

end subroutine diag_send_complete_infra

!> Sets the time that the simulation ends in the diag manager
subroutine diag_manager_set_time_end_infra(time)
  type(time_type),           optional, intent(in) :: time  !< The time the simulation ends

  integer :: days, secs

  if (present(time)) then
    call get_time(time, secs, days)
    call tim_diag_set_time_end(days, secs)
  endif

end subroutine diag_manager_set_time_end_infra

end module MOM_diag_manager_infra
