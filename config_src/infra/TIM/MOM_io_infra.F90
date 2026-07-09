!> This module contains a thin interface to the TIM (PIO-based) I/O code
module MOM_io_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms_helpers,      only : PE_here, root_PE, num_PEs, is_root_pe
use MOM_domain_infra,     only : MOM_domain_type, rescale_comp_data, AGRID, BGRID_NE, CGRID_NE
use MOM_domain_infra,     only : domain2d, domain1d, CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_error_infra,      only : MOM_err, NOTE, FATAL, WARNING

use mpp_domains_mod,      only : mpp_get_compute_domain, mpp_get_global_domain
use mpp_domains_mod,      only : mpp_get_data_domain
use mpp_mod,              only : stdout_if_root=>stdout
use mpp_mod,              only : mpp_get_current_pelist_name
use mpp_mod,              only : mpp_max
use iso_fortran_env,      only : int64

use tim_io_interface, only : tim_io_register_domain, tim_io_read_decomposed, cstr
use tim_io_interface, only : tim_io_read_plain
use tim_io_interface, only : tim_io_createfile, tim_io_def_axis, tim_io_def_var
use tim_io_interface, only : tim_io_put_global_att, tim_io_write_axis, tim_io_var_stagger
use tim_io_interface, only : tim_io_write_decomposed, tim_io_write_plain, tim_io_closefile
use tim_io_interface, only : tim_io_file_num_times, tim_io_file_time
use tim_io_interface, only : tim_io_finalize
use tim_io_interface, only : tim_io_file_exists, tim_io_file_info, tim_io_file_times
use tim_io_interface, only : tim_io_var_exists
use tim_io_interface, only : tim_io_file_var_name, tim_io_var_att, tim_io_var_sizes
use tim_io_interface, only : tim_io_read_slab
use, intrinsic :: iso_c_binding, only : c_double, c_int

implicit none ; private

! Duplication of FMS1 parameter values
! NOTE: Only kept to emulate FMS1 behavior, and may be removed in the future.
integer, parameter :: WRITEONLY_FILE = 100
integer, parameter :: READONLY_FILE = 101
integer, parameter :: APPEND_FILE = 102
integer, parameter :: OVERWRITE_FILE = 103
integer, parameter :: ASCII_FILE = 200
integer, parameter :: NETCDF_FILE = 203
integer, parameter :: SINGLE_FILE = 400
integer, parameter :: MULTIPLE = 401

! TIM PIO decomposition registry state
integer, parameter :: MAX_TIM_DOMAINS = 64 !< Max distinct decompositions memoized
integer :: n_tim_domains = 0          !< Number of registered decompositions
integer :: tim_dom_sig(7, MAX_TIM_DOMAINS) = 0 !< Signatures of registered decomps
integer :: tim_dom_handle(MAX_TIM_DOMAINS) = -1 !< TIM-side handles
character(len=64) :: tim_filename_suffix = "" !< Filename appendix (ensemble suffix)

! PROTOTYPE: seam-level read timing (timer brackets the whole
! read_field/read_vector body).
real(kind=8) :: seam_read_secs = 0.0 !< Accumulated read seconds on this rank
integer :: seam_read_count = 0       !< Number of timed read calls on this rank
real(kind=8) :: seam_write_secs = 0.0 !< Accumulated write seconds on this rank
integer :: seam_write_count = 0      !< Number of timed write calls on this rank

! These interfaces are actually implemented or have explicit interfaces in this file.
public :: open_file, open_ASCII_file, file_is_open, close_file, flush_file, file_exists
public :: get_file_info, get_file_fields, get_file_times, get_filename_suffix, set_filename_suffix
public :: read_field, read_vector, write_metadata, write_field
public :: field_exists, get_field_atts, get_field_size, read_field_chksum
public :: get_axis_data, set_axis_data
public :: io_infra_init, io_infra_end, MOM_namelist_file, check_namelist_error, write_version
public :: stdout_if_root
! TIM prototype: expose the decomposition handle for the diag seam.
public :: tim_get_domain_handle, tim_get_domain2d_handle
! These types act as containers for information about files, fields and axes, respectively,
! and may also wrap opaque types from the underlying infrastructure.
public :: file_type, fieldtype, axistype
! These are encoding constant parmeters.
public :: ASCII_FILE, NETCDF_FILE, SINGLE_FILE, MULTIPLE
public :: APPEND_FILE, READONLY_FILE, OVERWRITE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

!> Indicate whether a file exists, perhaps with domain decomposition
interface file_exists
  module procedure FMS_file_exists
  module procedure MOM_file_exists
end interface

!> Read a data field from a file
interface read_field
  module procedure read_field_4d
  module procedure read_field_3d, read_field_3d_region
  module procedure read_field_2d, read_field_2d_region
  module procedure read_field_1d, read_field_1d_int
  module procedure read_field_0d, read_field_0d_int
end interface

!> Write a registered field to an output file
interface write_field
  module procedure write_field_4d
  module procedure write_field_3d
  module procedure write_field_2d
  module procedure write_field_1d
  module procedure write_field_0d
  module procedure MOM_write_axis
end interface write_field

!> Read a pair of data fields representing the two components of a vector from a file
interface read_vector
  module procedure read_vector_3d
  module procedure read_vector_2d
end interface read_vector

!> Write metadata about a variable or axis to a file and store it for later reuse
interface write_metadata
  module procedure write_metadata_axis, write_metadata_field, write_metadata_global
end interface write_metadata

!> Close a file (or fileset).  If the file handle does not point to an open file,
!! close_file simply returns without doing anything.
interface close_file
  module procedure close_file_type, close_file_unit
end interface close_file

!> Type for holding a handle to an open file and related information
type :: file_type ; private
  integer :: unit = -1 !< The framework identfier or netCDF unit number of an output file
  character(len=:), allocatable :: filename !< The path to this file, if it is open
  logical :: open_to_read  = .false. !< If true, this file or fileset can be read
  logical :: open_to_write = .false. !< If true, this file or fileset can be written to
  integer :: num_times !< The number of time levels in this file
  real    :: file_time !< The time of the latest entry in the file.
  integer :: tim_fh = -1 !< TIM write-file handle (>=0 when TIM owns this file)
  logical :: tim_read = .false. !< True when TIM owns this read handle
end type file_type

!> This type is a container for information about a variable in a file.
type :: fieldtype ; private
  character(len=256)  :: name !< The name of this field in the files.
  character(len=:), allocatable :: longname !< The long name for this field
  character(len=:), allocatable :: units    !< The units for this field
  integer(kind=int64) :: chksum_read !< A checksum that has been read from a file
  logical :: valid_chksum !< If true, this field has a valid checksum value.
end type fieldtype

!> This type is a container for information about an axis in a file.
type :: axistype ; private
  character(len=256) :: name !< The name of this axis in the files.
  real, allocatable, dimension(:) :: ax_data !< The values of the data on the axis.
  logical :: domain_decomposed = .false.  !< True if axis is domain-decomposed
end type axistype

contains

!> Reads the checksum value for a field that was recorded in a file, along with a flag indicating
!! whether the file contained a valid checksum for this field.
subroutine read_field_chksum(field, chksum, valid_chksum)
  type(fieldtype),     intent(in)  :: field !< The field whose checksum attribute is to be read.
  integer(kind=int64), intent(out) :: chksum !< The checksum for the field.
  logical,             intent(out) :: valid_chksum  !< If true, chksum has been successfully read.

  chksum = -1
  valid_chksum = field%valid_chksum
  if (valid_chksum) chksum = field%chksum_read

end subroutine read_field_chksum

!> Returns true if the named file or its domain-decomposed variant exists.
logical function MOM_file_exists(filename, MOM_Domain)
  character(len=*),       intent(in) :: filename   !< The name of the file being inquired about
  type(MOM_domain_type),  intent(in) :: MOM_Domain !< The MOM_Domain that describes the decomposition

  ! Plain existence check (works for any file type; netCDF callers may omit .nc)
  inquire(file=trim(filename), exist=MOM_file_exists)
  if (.not. MOM_file_exists) &
    inquire(file=tim_norm_path(filename), exist=MOM_file_exists)
end function MOM_file_exists

!> Returns true if the named file or its domain-decomposed variant exists.
logical function FMS_file_exists(filename)
  character(len=*),         intent(in) :: filename  !< The name of the file being inquired about
  ! Determines whether a named file (or its decomposed variant) exists.

  ! Plain existence check: this is used for ASCII files (input.nml!) as well
  ! as netCDF files named with or without their .nc suffix.
  inquire(file=trim(filename), exist=FMS_file_exists)
  if (.not. FMS_file_exists) &
    inquire(file=tim_norm_path(filename), exist=FMS_file_exists)
end function FMS_file_exists

!> indicates whether an I/O handle is attached to an open file
logical function file_is_open(IO_handle)
  type(file_type), intent(in) :: IO_handle !< Handle to a file to inquire about

  file_is_open = ((IO_handle%unit >= 0) .or. IO_handle%tim_read .or. (IO_handle%tim_fh >= 0))
end function file_is_open

!> closes a file (or fileset).  If the file handle does not point to an open file,
!! close_file_type simply returns without doing anything.
subroutine close_file_type(IO_handle)
  type(file_type), intent(inout) :: IO_handle   !< The I/O handle for the file to be closed
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer (close = flush; count it as write time)

  t0_seam = seam_tic()
  if (IO_handle%tim_read) then
    ! The context cache keeps the file open (cheap reopens); just release.
    IO_handle%tim_read = .false.
  elseif (IO_handle%tim_fh >= 0) then
    if (tim_io_closefile(IO_handle%tim_fh) /= 0) &
      call MOM_err(FATAL, "TIM: error closing "//trim(IO_handle%filename))
    IO_handle%tim_fh = -1
  endif
  if (allocated(IO_handle%filename)) deallocate(IO_handle%filename)
  IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .false.
  IO_handle%num_times = 0 ; IO_handle%file_time = 0.0
  call seam_toc_w(t0_seam)
end subroutine close_file_type

! TODO: close_file_unit is only used for ASCII files, which are opened outside
! of the framework, so this could probably be removed, and those calls could
! just be replaced with close(unit).

!> closes a file.  If the unit does not point to an open file,
!! close_file_unit simply returns without doing anything.
subroutine close_file_unit(iounit)
  integer, intent(inout) :: iounit   !< The I/O unit for the file to be closed

  logical :: unit_is_open

  inquire(iounit, opened=unit_is_open)
  if (unit_is_open) close(iounit)
end subroutine close_file_unit

!> Ensure that the output stream associated with a file handle is fully sent to disk.
subroutine flush_file(IO_handle)
  type(file_type), intent(in) :: IO_handle    !< The I/O handle for the file to flush

  ! TIM/PIO does not expose an explicit flush at this seam; data are
  ! flushed when the file is closed.
end subroutine flush_file

!> Initialize the underlying I/O infrastructure
subroutine io_infra_init(maxunits)
  integer,   optional, intent(in) :: maxunits !< An optional maximum number of file
                                              !! unit numbers that can be used.

  ! TIM's I/O context is created by tim_io_init (from MOM_infra_init), so
  ! this is a null function.
end subroutine io_infra_init

!> Gracefully close out and terminate the underlying I/O infrastructure
subroutine io_infra_end()
  ! PROTOTYPE: close TIM's cached
  ! files and finalize the PIO iosystem BEFORE MPI_Finalize, and report
  ! seam read timing.
  real(kind=8) :: max_secs
  character(len=160) :: mesg
  call tim_io_finalize()
  max_secs = seam_read_secs
  call mpp_max(max_secs)
  if (is_root_pe() .and. (seam_read_count > 0)) then
    write(mesg, '("MOM_io_infra seam reads: ",I6," calls, root ",F10.3," s, max ",F10.3," s")') &
      seam_read_count, seam_read_secs, max_secs
    call MOM_err(NOTE, trim(mesg))
  endif
  max_secs = seam_write_secs
  call mpp_max(max_secs)
  if (is_root_pe() .and. (seam_write_count > 0)) then
    write(mesg, '("MOM_io_infra seam writes:",I6," calls, root ",F10.3," s, max ",F10.3," s")') &
      seam_write_count, seam_write_secs, max_secs
    call MOM_err(NOTE, trim(mesg))
  endif
end subroutine io_infra_end

!> Start a seam read timer (returns current time in seconds).
function seam_tic() result(t0)
  real(kind=8) :: t0
  integer(kind=8) :: c, cr
  call system_clock(c, cr)
  t0 = real(c, kind=8) / real(cr, kind=8)
end function seam_tic

!> Stop a seam write timer and accumulate.
subroutine seam_toc_w(t0)
  real(kind=8), intent(in) :: t0
  integer(kind=8) :: c, cr
  call system_clock(c, cr)
  seam_write_secs = seam_write_secs + (real(c, kind=8) / real(cr, kind=8) - t0)
  seam_write_count = seam_write_count + 1
end subroutine seam_toc_w

!> Stop a seam read timer and accumulate.
subroutine seam_toc(t0)
  real(kind=8), intent(in) :: t0
  integer(kind=8) :: c, cr
  call system_clock(c, cr)
  seam_read_secs = seam_read_secs + (real(c, kind=8) / real(cr, kind=8) - t0)
  seam_read_count = seam_read_count + 1
end subroutine seam_toc

!> Open a single namelist file that is potentially readable by all PEs.
function MOM_namelist_file(filepath) result(iounit)
  character(len=*), optional, intent(in) :: filepath
    !< The file to open, by default "input.nml".
  integer                                :: iounit
    !< The opened unit number of the namelist file

  character(len=:), allocatable :: nmlpath
    ! Namelist path
  character(len=:), allocatable :: nmlpath_pe
    ! Hypothetical namelist path exclusive to the current PE list

  if (present(filepath)) then
    nmlpath = trim(filepath)
  else
    ! FMS1 first checks for a namelist unique to the PE list, `input_{}.nml`.
    ! If not found, it defaults to `input.nml`.
    nmlpath_pe = 'input_' // trim(mpp_get_current_pelist_name()) // '.nml'
    if (file_exists(nmlpath_pe)) then
      nmlpath = nmlpath_pe
    else
      nmlpath = 'input.nml'
    endif
  endif
  call open_ASCII_file(iounit, nmlpath, action=READONLY_FILE)
end function MOM_namelist_file

!> Checks the iostat argument that is returned after reading a namelist variable and writes a
!! message if there is an error.
subroutine check_namelist_error(IOstat, nml_name)
  integer,          intent(in) :: IOstat   !< An I/O status field from a namelist read call
  character(len=*), intent(in) :: nml_name !< The name of the namelist
  ! Native replacement for FMS check_nml_error. A negative status means the
  ! group was not found before end-of-file, which FMS tolerates (optional
  ! namelists); a positive status is a real read/syntax error and is fatal.
  if (IOstat > 0) call MOM_err(FATAL, "Error while reading namelist "//trim(nml_name))
end subroutine check_namelist_error

!> Write a file version number to the log file or other output file
subroutine write_version(version, tag, unit)
  character(len=*),           intent(in) :: version !< A string that contains the routine name and version
  character(len=*), optional, intent(in) :: tag  !< A tag name to add to the message
  integer,          optional, intent(in) :: unit !< An alternate unit number for output

  ! Native replacement for FMS write_version_number: a one-line version
  ! banner from the root PE (to stdout, or `unit` when given).
  integer :: out
  out = stdout_if_root()
  if (present(unit)) out = unit
  if (out >= 0) then
    if (present(tag)) then
      write(out,'(a)') trim(version)//" "//trim(tag)
    else
      write(out,'(a)') trim(version)
    endif
  endif
end subroutine write_version

!> open_file opens a file for parallel or single-file I/O.
subroutine open_file(IO_handle, filename, action, MOM_domain, threading, fileset)
  type(file_type),          intent(inout) :: IO_handle !< The handle for the opened file
  character(len=*),         intent(in)    :: filename !< The path name of the file being opened
  integer,        optional, intent(in)    :: action !< A flag indicating whether the file can be read
                                                    !! or written to and how to handle existing files.
                                                    !! The default is WRITE_ONLY.
  type(MOM_domain_type), &
                  optional, intent(in)    :: MOM_Domain !< A MOM_Domain that describes the decomposition
  integer,        optional, intent(in)    :: threading !< A flag indicating whether one (SINGLE_FILE)
                                                    !! or multiple PEs (MULTIPLE) participate in I/O.
                                                    !! With the default, the root PE does I/O.
  integer,        optional, intent(in)    :: fileset !< A flag indicating whether multiple PEs doing I/O due
                                                    !! to threading=MULTIPLE write to the same file (SINGLE_FILE)
                                                    !! or to one file per PE (MULTIPLE, the default).

  ! Local variables
  integer :: file_mode       ! An integer that encodes whether the file is to be opened for
                             ! reading, writing or appending
  character(len=:), allocatable :: filename_tmp  ! A copy of filename with .nc appended if necessary.
  integer :: index_nc

  if (IO_handle%open_to_write) then
    call MOM_err(WARNING, "open_file called for file "//trim(filename)//&
        " with an IO_handle that is already open to to write.")
    return
  endif
  if (IO_handle%open_to_read) then
    call MOM_err(FATAL, "open_file called for file "//trim(filename)//&
        " with an IO_handle that is already open to to read.")
  endif

  file_mode = WRITEONLY_FILE ; if (present(action)) file_mode = action

  ! Domains are currently required, matching the historical behavior of this interface.
  if (.not. present(MOM_Domain)) &
    call MOM_err(FATAL, 'open_file: a domain input is required.')

  ! The FMS1 interface automatically appended .nc if necessary; retain that behavior.
  index_nc = index(trim(filename), ".nc")
  if (index_nc > 0) then
    filename_tmp = trim(filename)
  else
    filename_tmp = trim(filename)//".nc"
    if (is_root_PE()) call MOM_err(WARNING, "Open_file is appending .nc to the filename "//trim(filename))
  endif

  if (file_mode == READONLY_FILE) then
    ! TIM PIO read handle (metadata/query use; data reads are served by the
    ! stateless read calls against the same cached file)
    block
      integer :: nd_, nv_, nt_
      if (tim_io_file_exists(cstr(filename_tmp)) == 0) &
        call MOM_err(FATAL, "TIM: unable to open for read: "//trim(filename_tmp))
      if (tim_io_file_info(cstr(filename_tmp), nd_, nv_, nt_) /= 0) &
        call MOM_err(FATAL, "TIM: file_info failed: "//trim(filename_tmp))
      IO_handle%tim_read = .true.
      IO_handle%filename = trim(filename_tmp)
      IO_handle%open_to_read = .true. ; IO_handle%open_to_write = .false.
      IO_handle%num_times = nt_
      IO_handle%file_time = 0.0
    end block
  elseif ((file_mode == WRITEONLY_FILE) .or. (file_mode == OVERWRITE_FILE) .or. &
          (file_mode == APPEND_FILE)) then
    ! TIM PIO write path
    block
      integer :: tim_mode
      tim_mode = 0
      if (file_mode == OVERWRITE_FILE) tim_mode = 1
      if (file_mode == APPEND_FILE) tim_mode = 2
      IO_handle%tim_fh = tim_io_createfile(cstr(filename_tmp), &
                                           tim_get_domain_handle(MOM_Domain), tim_mode)
      if (IO_handle%tim_fh < 0) &
        call MOM_err(FATAL, "TIM: unable to open for write: "//trim(filename_tmp))
      IO_handle%filename = trim(filename)
      IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .true.
      IO_handle%num_times = tim_io_file_num_times(IO_handle%tim_fh)
      IO_handle%file_time = tim_io_file_time(IO_handle%tim_fh)
    end block
  else
    call MOM_err(FATAL, "open_file called with unrecognized action.")
  endif

end subroutine open_file

!> open_file opens an ascii file for parallel or single-file I/O using Fortran read and write calls.
subroutine open_ASCII_file(unit, file, action, threading, fileset)
  integer,                  intent(out) :: unit   !< The I/O unit for the opened file
  character(len=*),         intent(in)  :: file   !< The name of the file being opened
  integer,        optional, intent(in)  :: action !< A flag indicating whether the file can be read
                                                  !! or written to and how to handle existing files.
  integer,        optional, intent(in)  :: threading !< A flag indicating whether one (SINGLE_FILE)
                                                  !! or multiple PEs (MULTIPLE) participate in I/O.
                                                  !! With the default, the root PE does I/O.
  integer,        optional, intent(in)  :: fileset !< A flag indicating whether multiple PEs doing I/O due
                                                  !! to threading=MULTIPLE write to the same file (SINGLE_FILE)
                                                  !! or to one file per PE (MULTIPLE, the default).

  integer :: action_flag
  integer :: threading_flag
  integer :: fileset_flag
  logical :: exists
  logical :: is_open
  character(len=6) :: action_arg, position_arg
  character(len=:), allocatable :: filename

  ! NOTE: This function is written to emulate the original behavior of mpp_open
  !   from the FMS1 library, on which the MOM API is still based.  Much of this
  !   can be removed if we choose to drop this compatibility, but for now we
  !   try to retain as much as possible.

  ! NOTE: Default FMS1 I/O settings are summarized below.
  !
  !   access: Fortran and mpp_open default to SEQUENTIAL.
  !   form:   The Fortran and mpp_open default (for MPP_ASCII) is FORMATTED.
  !   recl:   mpp_open uses Fortran defaults when unset, so can be ignored.
  !   ios:    FMS1 allowed this to be caught, but we do not support it.
  !   action/position:  In mpp_open, these are inferred from `action`.
  !
  !     MOM flag        FMS1 flag     action    position
  !     --------        --------      ------    --------
  !     READONLY_FILE   MPP_RDONLY    READ      REWIND
  !     WRITEONLY_FILE  MPP_WRONLY    WRITE     REWIND
  !     OVERWRITE_FILE  MPP_OVERWR    WRITE     REWIND
  !     APPEND_FILE     MPP_APPEND    WRITE     APPEND
  !
  ! From this, we can omit `access`, `form`, and `recl`, and can construct
  !   `action` and `position` from the input arguments.

  ! I/O configuration

  action_flag = WRITEONLY_FILE
  if (present(action)) action_flag = action

  action_arg = 'write'
  if (action_flag == READONLY_FILE) action_arg = 'read'

  position_arg = 'rewind'
  if (action_flag == APPEND_FILE) position_arg = 'append'

  ! Threading configuration

  threading_flag = SINGLE_FILE
  if (present(threading)) threading_flag = threading

  fileset_flag = MULTIPLE
  if (present(fileset)) fileset_flag = fileset

  ! Force fileset to be consistent with threading (as in FMS1)
  if (threading_flag == SINGLE_FILE) fileset_flag = SINGLE_FILE

  ! Construct the distributed filename, if needed
  filename = file
  if (fileset_flag == MULTIPLE) then
    if (num_PEs() > 10000) then
      write(filename, '(a,".",i6.6)') trim(filename), PE_here() - root_PE()
    else
      write(filename, '(a,".",i4.4)') trim(filename), PE_here() - root_PE()
    endif
  endif

  inquire(file=filename, exist=exists)
  if (exists .and. action_flag == WRITEONLY_FILE) &
    call MOM_err(WARNING, 'open_ASCII_file: File ' // trim(filename) // &
                            ' opened WRITEONLY already exists!')

  open(newunit=unit, file=filename, action=trim(action_arg), &
       position=trim(position_arg))

  ! This checks if open() failed but did not raise a runtime error.
  inquire(unit, opened=is_open)
  if (.not. is_open) &
    call MOM_err(FATAL, &
        'open_ASCII_file: File "' // trim(filename) // '" failed to open.')

  ! NOTE: There are two possible mpp_write_meta functions in FMS1:
  ! - call mpp_write_meta( unit, 'filename', cval=mpp_file(unit)%name)
  ! - call mpp_write_meta( unit, 'NumFilesInSet', ival=nfiles)
  ! I'm not convinced we actually want these, but note them here in case.
end subroutine open_ASCII_file


!> Provide a string to append to filenames, to differentiate ensemble members, for example.
subroutine get_filename_suffix(suffix)
  character(len=*), intent(out) :: suffix !< A string to append to filenames

  suffix = trim(tim_filename_suffix)
end subroutine get_filename_suffix

!> Store a string to append to filenames (e.g. an ensemble-member suffix);
!! native replacement for the fms2_io filename appendix.
subroutine set_filename_suffix(suffix)
  character(len=*), intent(in) :: suffix !< A string to append to filenames

  tim_filename_suffix = trim(suffix)
end subroutine set_filename_suffix


!> Get information about the number of dimensions, variables and time levels
!! in the file associated with an open file unit
subroutine get_file_info(IO_handle, ndim, nvar, ntime)
  type(file_type),    intent(in)  :: IO_handle !< Handle for a file that is open for I/O
  integer,  optional, intent(out) :: ndim  !< The number of dimensions in the file
  integer,  optional, intent(out) :: nvar  !< The number of variables in the file
  integer,  optional, intent(out) :: ntime !< The number of time levels in the file

  ! Local variables
  integer :: nd_, nv_, nt_

  if (.not. IO_handle%tim_read) &
    call MOM_err(FATAL, "get_file_info called for a file that is not open for reading.")

  if (tim_io_file_info(cstr(IO_handle%filename), nd_, nv_, nt_) /= 0) &
    call MOM_err(FATAL, "TIM: file_info failed: "//trim(IO_handle%filename))
  if (present(ndim)) ndim = nd_
  if (present(nvar)) nvar = nv_
  if (present(ntime)) ntime = nt_
end subroutine get_file_info


!> Get the times of records from a file
subroutine get_file_times(IO_handle, time_values, ntime)
  type(file_type),                 intent(in)    :: IO_handle !< Handle for a file that is open for I/O
  real, allocatable, dimension(:), intent(inout) :: time_values !< The real times for the records in file.
  integer,               optional, intent(out)   :: ntime !< The number of time levels in the file

  integer :: ntimes  ! The number of time levels in the file
  real(kind=c_double), allocatable :: tbuf(:)

  !### Modify this routine to optionally convert to time_type, using information about the dimensions?

  if (allocated(time_values)) deallocate(time_values)
  call get_file_info(IO_handle, ntime=ntimes)
  if (present(ntime)) ntime = ntimes
  if (ntimes <= 0) return
  allocate(time_values(ntimes))
  allocate(tbuf(ntimes))
  if (tim_io_file_times(cstr(IO_handle%filename), tbuf, ntimes) /= 0) &
    call MOM_err(FATAL, "TIM: file_times failed: "//trim(IO_handle%filename))
  time_values(:) = tbuf(:)
  deallocate(tbuf)
end subroutine get_file_times

!> Set up the field information (e.g., names and metadata) for all of the variables in a file.  The
!! argument fields must be allocated with a size that matches the number of variables in a file.
subroutine get_file_fields(IO_handle, fields)
  type(file_type),               intent(in)    :: IO_handle !< Handle for a file that is open for I/O
  type(fieldtype), dimension(:), intent(inout) :: fields !< Field-type descriptions of all of
                                                         !! the variables in a file.
  character(len=256),  dimension(size(fields)) :: var_names ! The names of all variables
  character(len=256)  :: units    ! The units of a variable as recorded in the file
  character(len=2048) :: longname ! The long-name of a variable as recorded in the file
  character(len=64)   :: checksum_char ! The hexadecimal checksum read from the file
  integer(kind=int64), dimension(3) :: checksum_file ! The checksums for a variable in the file
  integer :: nvar  ! The number of variables in the file
  integer :: i

  nvar = size(fields)

  if (.not. IO_handle%tim_read) &
    call MOM_err(FATAL, "get_file_fields called for a file that is not open for reading.")

  do i=1,nvar
    var_names(i) = ""
    if (tim_io_file_var_name(cstr(IO_handle%filename), i, var_names(i), 256) /= 0) &
      call MOM_err(FATAL, "TIM: var_name failed: "//trim(IO_handle%filename))
    call tim_cstr_to_f(var_names(i))
    fields(i)%name = trim(var_names(i))
    longname = ""
    if (tim_io_var_att(cstr(IO_handle%filename), cstr(var_names(i)), cstr("long_name"), &
                       longname, 2048) == 0) call tim_cstr_to_f(longname)
    fields(i)%longname = trim(longname)
    units = ""
    if (tim_io_var_att(cstr(IO_handle%filename), cstr(var_names(i)), cstr("units"), &
                       units, 256) == 0) call tim_cstr_to_f(units)
    fields(i)%units = trim(units)
    checksum_char = ""
    fields(i)%valid_chksum = (tim_io_var_att(cstr(IO_handle%filename), cstr(var_names(i)), &
                              cstr("checksum"), checksum_char, 64) == 0)
    fields(i)%chksum_read = -1
    if (fields(i)%valid_chksum) then
      call tim_cstr_to_f(checksum_char)
      read (checksum_char(1:16), '(Z16)') checksum_file(1)
      fields(i)%chksum_read = checksum_file(1)
    endif
  enddo
end subroutine get_file_fields

!> Extract information from a field type, as stored or as found in a file
subroutine get_field_atts(field, name, units, longname, checksum)
  type(fieldtype),            intent(in)  :: field !< The field to extract information from
  character(len=*), optional, intent(out) :: name  !< The variable name
  character(len=*), optional, intent(out) :: units !< The units of the variable
  character(len=*), optional, intent(out) :: longname  !< The long name of the variable
  integer(kind=int64),  dimension(:), &
                    optional, intent(out) :: checksum !< The checksums of the variable in a file

  if (present(name)) name = trim(field%name)
  if (present(units)) units = trim(field%units)
  if (present(longname)) longname = trim(field%longname)
  if (present(checksum)) checksum = field%chksum_read

end subroutine get_field_atts

!> Field_exists returns true if the field indicated by field_name is present in the
!! file file_name.  If file_name does not exist, it returns false.
function field_exists(filename, field_name, domain, no_domain, MOM_domain)
  character(len=*),                 intent(in) :: filename   !< The name of the file being inquired about
  character(len=*),                 intent(in) :: field_name !< The name of the field being sought
  type(domain2d), target, optional, intent(in) :: domain     !< A domain2d type that describes the decomposition
  logical,                optional, intent(in) :: no_domain  !< This file does not use domain decomposition
  type(MOM_domain_type),  optional, intent(in) :: MOM_Domain !< A MOM_Domain that describes the decomposition
  logical                                      :: field_exists !< True if filename exists and field_name is in filename

  ! Local variables
  logical :: domainless      ! If true, this file does not use a domain-decomposed file.

  domainless = .not.(present(MOM_domain) .or. present(domain))
  if (present(no_domain)) then
    if (domainless .and. .not.no_domain) call MOM_err(FATAL, &
        "field_exists: When no_domain is present and false, a domain must be supplied in query about "//&
        trim(field_name)//" in file "//trim(filename))
  endif

  field_exists = .false.
  if (tim_io_file_exists(cstr(tim_norm_path(filename))) /= 0) &
    field_exists = (tim_io_var_exists(cstr(tim_norm_path(filename)), cstr(field_name)) /= 0)
end function field_exists

!> Given filename and fieldname, this subroutine returns the size of the field in the file
subroutine get_field_size(filename, fieldname, sizes, field_found, no_domain)
  character(len=*),      intent(in)    :: filename  !< The name of the file to read
  character(len=*),      intent(in)    :: fieldname !< The name of the variable whose sizes are returned
  integer, dimension(:), intent(inout) :: sizes     !< The sizes of the variable in each dimension
  logical,     optional, intent(out)   :: field_found !< This indicates whether the field was found in
                                                    !! the input file.  Without this argument, there
                                                    !! is a fatal error if the field is not found.
  logical,     optional, intent(in)    :: no_domain !< If present and true, do not check for file
                                                    !! names with an appended tile number
  ! Local variables
  logical :: field_exists    ! True if filename exists and field_name is in filename
  integer :: i
  integer :: sizes4(4), nd_

  field_exists = .false.
  sizes4(:) = 1
  if (tim_io_file_exists(cstr(tim_norm_path(filename))) /= 0) then
    nd_ = tim_io_var_sizes(cstr(tim_norm_path(filename)), cstr(fieldname), sizes4)
    if (nd_ >= 0) then
      field_exists = .true.
      sizes(:) = 1
      do i=1,min(nd_, size(sizes))
        sizes(i) = sizes4(i)
      enddo
    endif
  endif
  if (present(field_found)) field_found = field_exists
end subroutine get_field_size


!> Extracts and returns the axis data stored in an axistype.
subroutine get_axis_data(axis, axis_name, axis_data)
  type(axistype), intent(in) :: axis   !< Infra axis
  character(len=256), intent(out) :: axis_name   !< Axis name
  real, dimension(:), intent(out) :: axis_data   !< Axis points

  integer :: i

  if (allocated(axis%ax_data)) then
    if (size(axis%ax_data) > size(axis_data)) &
      call MOM_err(FATAL, "get_axis_data called with too small of an " &
          // "output data array for " // trim(axis%name) // ".")
    do i=1,size(axis%ax_data)
      axis_data(i) = axis%ax_data(i)
    enddo
  endif

  axis_name = axis%name

end subroutine get_axis_data

!> Return a new axistype based on axis specs
subroutine set_axis_data(axis, axis_name, axis_data)
  type(axistype), intent(inout) :: axis
    !< Target axis
  character(len=256), intent(in) :: axis_name
    !< Target axis name
  real, intent(in) :: axis_data(:)
    !< Target axis values

  axis%name = axis_name

  if (allocated(axis%ax_data)) deallocate(axis%ax_data)
  allocate(axis%ax_data(size(axis_data)))

  axis%ax_data(:) = axis_data(:)

  ! NOTE: We do not yet consider domain-decomposed axes.
  axis%domain_decomposed = .false.
end subroutine set_axis_data


!> This routine uses the TIM PIO path to read a scalar named
!! "fieldname" from a single or domain-decomposed file "filename".
subroutine read_field_0d(filename, fieldname, data, timelevel, scale, MOM_Domain, &
                         global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real,                   intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays, but
                                                     !! with the TIM I/O interfaces this does not matter.

  ! Local variables
  real(kind=c_double) :: buf0(1)   ! TIM read buffer

  call tim_read_plain(filename, fieldname, 1, buf0, timelevel)
  data = buf0(1)
  if (present(scale)) then ; if (scale /= 1.0) data = scale*data ; endif

end subroutine read_field_0d

!> This routine uses the TIM PIO path to read a 1-D data field named
!! "fieldname" from a single or domain-decomposed file "filename".
subroutine read_field_1d(filename, fieldname, data, timelevel, scale, MOM_Domain, &
                         global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before they are returned.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays, but
                                                     !! with the TIM I/O interfaces this does not matter.

  ! Local variables
  real(kind=c_double), allocatable :: buf1(:) ! TIM read buffer

  allocate(buf1(size(data)))
  call tim_read_plain(filename, fieldname, size(data), buf1, timelevel)
  data(:) = buf1(:)
  deallocate(buf1)
  if (present(scale)) then ; if (scale /= 1.0) data(:) = scale*data(:) ; endif

end subroutine read_field_1d

!> This routine uses the TIM PIO path to read a distributed
!! 2-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine read_field_2d(filename, fieldname, data, MOM_Domain, &
                         timelevel, position, scale, global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data      !< The 2-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays, but
                                                     !! with the TIM I/O interfaces this does not matter.

  ! Local variables
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer

  t0_seam = seam_tic()
  call tim_read_field_dd(filename, fieldname, data2d=data, MOM_Domain=MOM_Domain, &
                         timelevel=timelevel, position=position, scale=scale)
  call seam_toc(t0_seam)

end subroutine read_field_2d

!> This routine uses the TIM PIO path to read a region from a distributed or
!! global 2-D data field named "fieldname" from file "filename".
subroutine read_field_2d_region(filename, fieldname, data, start, nread, MOM_domain, &
                                no_domain, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data      !< The 2-dimensional array into which the data
                                                     !! should be read
  integer, dimension(:),  intent(in)    :: start     !< The starting index to read in each of 4
                                                     !! dimensions.  For this 2-d read, the 3rd
                                                     !! and 4th values are always 1.
  integer, dimension(:),  intent(in)    :: nread     !< The number of points to read in each of 4
                                                     !! dimensions.  For this 2-d read, the 3rd
                                                     !! and 4th values are always 1.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: no_domain !< If present and true, this variable does not
                                                     !! use domain decomposion.
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

  ! Local variables
  real(kind=c_double), allocatable :: sbuf(:)
  integer :: st4(4), nr4(4), k_

  st4(:) = 1 ; nr4(:) = 1
  do k_=1,min(4,size(start)) ; st4(k_) = start(k_) ; enddo
  do k_=1,min(4,size(nread)) ; nr4(k_) = nread(k_) ; enddo
  allocate(sbuf(size(data)))
  if (tim_io_read_slab(cstr(tim_norm_path(filename)), cstr(fieldname), st4, nr4, sbuf) /= 0) &
    call MOM_err(FATAL, "TIM: region read failed: "//trim(fieldname)//" from "//trim(filename))
  data = reshape(sbuf, shape(data))
  deallocate(sbuf)
  if (present(scale)) then ; if (scale /= 1.0) then
    if (present(MOM_Domain)) then
      call rescale_comp_data(MOM_Domain, data, scale)
    else
      ! Dangerously rescale the whole array
      data(:,:) = scale*data(:,:)
    endif
  endif ; endif

end subroutine read_field_2d_region

!> This routine uses the TIM PIO path to read a distributed
!! 3-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine read_field_3d(filename, fieldname, data, MOM_Domain, &
                         timelevel, position, scale, global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data      !< The 3-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays, but
                                                     !! with the TIM I/O interfaces this does not matter.

  ! Local variables
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer

  t0_seam = seam_tic()
  call tim_read_field_dd(filename, fieldname, data3d=data, MOM_Domain=MOM_Domain, &
                         timelevel=timelevel, position=position, scale=scale)
  call seam_toc(t0_seam)

end subroutine read_field_3d

!> This routine uses the TIM PIO path to read a region from a distributed or
!! global 3-D data field named "fieldname" from file "filename".
subroutine read_field_3d_region(filename, fieldname, data, start, nread, MOM_domain, &
                                no_domain, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:),   intent(inout) :: data      !< The 3-dimensional array into which the data
                                                     !! should be read
  integer, dimension(:),  intent(in)    :: start     !< The starting index to read in each of 3
                                                     !! dimensions.  For this 3-d read, the
                                                     !! 4th value is always 1.
  integer, dimension(:),  intent(in)    :: nread     !< The number of points to read in each of 4
                                                     !! dimensions.  For this 3-d read, the
                                                     !! 4th values are always 1.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: no_domain !< If present and true, this variable does not
                                                     !! use domain decomposion.
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

  ! Local variables
  real(kind=c_double), allocatable :: sbuf(:)
  integer :: st4(4), nr4(4), k_

  st4(:) = 1 ; nr4(:) = 1
  do k_=1,min(4,size(start)) ; st4(k_) = start(k_) ; enddo
  do k_=1,min(4,size(nread)) ; nr4(k_) = nread(k_) ; enddo
  allocate(sbuf(size(data)))
  if (tim_io_read_slab(cstr(tim_norm_path(filename)), cstr(fieldname), st4, nr4, sbuf) /= 0) &
    call MOM_err(FATAL, "TIM: region read failed: "//trim(fieldname)//" from "//trim(filename))
  data = reshape(sbuf, shape(data))
  deallocate(sbuf)
  if (present(scale)) then ; if (scale /= 1.0) then
    if (present(MOM_Domain)) then
      call rescale_comp_data(MOM_Domain, data, scale)
    else
      ! Dangerously rescale the whole array
      data(:,:,:) = scale*data(:,:,:)
    endif
  endif ; endif

end subroutine read_field_3d_region

!> This routine uses the TIM PIO path to read a distributed
!! 4-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine read_field_4d(filename, fieldname, data, MOM_Domain, &
                         timelevel, position, scale, global_file)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(inout) :: data    !< The 4-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file


  ! Local variables
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer

  t0_seam = seam_tic()
  call tim_read_field_dd(filename, fieldname, data4d=data, MOM_Domain=MOM_Domain, &
                         timelevel=timelevel, position=position, scale=scale)
  call seam_toc(t0_seam)

end subroutine read_field_4d

!> This routine uses the TIM PIO path to read a scalar integer
!! data field named "fieldname" from file "filename".
subroutine read_field_0d_int(filename, fieldname, data, timelevel)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  integer,                intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read

  ! Local variables
  real(kind=c_double) :: buf0(1)   ! TIM read buffer

  ! This routine might not be needed for MOM6.

  call tim_read_plain(filename, fieldname, 1, buf0, timelevel)
  data = nint(buf0(1))
end subroutine read_field_0d_int

!> This routine uses the TIM PIO path to read a 1-D integer
!! data field named "fieldname" from file "filename".
subroutine read_field_1d_int(filename, fieldname, data, timelevel)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  integer, dimension(:),  intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read

  ! Local variables
  real(kind=c_double), allocatable :: buf1(:) ! TIM read buffer

  ! This routine might not be needed for MOM6.

  allocate(buf1(size(data)))
  call tim_read_plain(filename, fieldname, size(data), buf1, timelevel)
  data(:) = nint(buf1(:))
  deallocate(buf1)
end subroutine read_field_1d_int


!> This routine uses the TIM PIO path to read a pair of distributed
!! 2-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine read_vector_2d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                          timelevel, stagger, scalar_pair, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:),   intent(inout) :: u_data    !< The 2 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:),   intent(inout) :: v_data    !< The 2 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  logical,      optional, intent(in)    :: scalar_pair !< If true, a pair of scalars are to be read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  ! Local variables
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer
  integer :: u_pos, v_pos           ! Flags indicating the positions of the u- and v- components.

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  t0_seam = seam_tic()
  call tim_read_field_dd(filename, u_fieldname, data2d=u_data, MOM_Domain=MOM_Domain, &
                         timelevel=timelevel, position=u_pos, scale=scale)
  call tim_read_field_dd(filename, v_fieldname, data2d=v_data, MOM_Domain=MOM_Domain, &
                         timelevel=timelevel, position=v_pos, scale=scale)
  call seam_toc(t0_seam)

end subroutine read_vector_2d

!> This routine uses the TIM PIO path to read a pair of distributed
!! 3-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine read_vector_3d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                          timelevel, stagger, scalar_pair, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:,:), intent(inout) :: u_data    !< The 3 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:,:), intent(inout) :: v_data    !< The 3 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  logical,      optional, intent(in)    :: scalar_pair !< If true, a pair of scalars are to be read.
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.

  ! Local variables
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer
  integer :: u_pos, v_pos           ! Flags indicating the positions of the u- and v- components.

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  t0_seam = seam_tic()
  call tim_read_field_dd(filename, u_fieldname, data3d=u_data, MOM_Domain=MOM_Domain, &
                         timelevel=timelevel, position=u_pos, scale=scale)
  call tim_read_field_dd(filename, v_fieldname, data3d=v_data, MOM_Domain=MOM_Domain, &
                         timelevel=timelevel, position=v_pos, scale=scale)
  call seam_toc(t0_seam)

end subroutine read_vector_3d


!> Write a 4d field to an output file.
subroutine write_field_4d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),          intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),          intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),    intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:,:), intent(inout) :: field      !< Field to write
  real,           optional, intent(in)    :: tstamp     !< Model time of this field
  integer,        optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,           optional, intent(in)    :: fill_value !< Missing data fill value


  ! Local variables
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "write_field_4d: file is not open for writing through TIM.")

  t0_seam = seam_tic()
  call tim_write_field_dd(IO_handle, field_md%name, MOM_domain, data4d=field, tstamp=tstamp)
  call seam_toc_w(t0_seam)
end subroutine write_field_4d

!> Write a 3d field to an output file.
subroutine write_field_3d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:), intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

  ! Local variables
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "write_field_3d: file is not open for writing through TIM.")

  t0_seam = seam_tic()
  call tim_write_field_dd(IO_handle, field_md%name, MOM_domain, data3d=field, tstamp=tstamp)
  call seam_toc_w(t0_seam)
end subroutine write_field_3d

!> Write a 2d field to an output file.
subroutine write_field_2d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:),   intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

  ! Local variables
  real(kind=8) :: t0_seam ! PROTOTYPE seam timer

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "write_field_2d: file is not open for writing through TIM.")

  t0_seam = seam_tic()
  call tim_write_field_dd(IO_handle, field_md%name, MOM_domain, data2d=field, tstamp=tstamp)
  call seam_toc_w(t0_seam)
end subroutine write_field_2d

!> Write a 1d field to an output file.
subroutine write_field_1d(IO_handle, field_md, field, tstamp)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real, dimension(:),     intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field

  ! Local variables
  real(kind=c_double), allocatable :: buf1(:)

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "write_field_1d: file is not open for writing through TIM.")

  allocate(buf1(size(field))) ; buf1(:) = field(:)
  call tim_write_plain_wrap(IO_handle, field_md%name, buf1, size(field), tstamp)
  deallocate(buf1)
end subroutine write_field_1d

!> Write a 0d field to an output file.
subroutine write_field_0d(IO_handle, field_md, field, tstamp)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real,                   intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field

  ! Local variables
  real(kind=c_double) :: buf0(1)

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "write_field_0d: file is not open for writing through TIM.")

  buf0(1) = field
  call tim_write_plain_wrap(IO_handle, field_md%name, buf0, 1, tstamp)
end subroutine write_field_0d

!> Write the data for an axis
subroutine MOM_write_axis(IO_handle, axis)
  type(file_type), intent(in) :: IO_handle  !< Handle for a file that is open for writing
  type(axistype),  intent(in) :: axis       !< An axis type variable with information to write

  real(kind=c_double), allocatable :: axbuf(:)

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "MOM_write_axis: file is not open for writing through TIM.")

  allocate(axbuf(size(axis%ax_data))) ; axbuf(:) = axis%ax_data(:)
  if (tim_io_write_axis(IO_handle%tim_fh, cstr(axis%name), axbuf, size(axbuf)) /= 0) &
    call MOM_err(FATAL, "TIM: axis write failed: "//trim(axis%name))
  deallocate(axbuf)
end subroutine MOM_write_axis

!> Store information about an axis in a previously defined axistype and write this
!! information to the file indicated by unit.
subroutine write_metadata_axis(IO_handle, axis, name, units, longname, cartesian, sense, domain, data, &
                               edge_axis, calendar)
  type(file_type),            intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(axistype),             intent(inout) :: axis  !< The axistype where this information is stored.
  character(len=*),           intent(in)    :: name  !< The name in the file of this axis
  character(len=*),           intent(in)    :: units !< The units of this axis
  character(len=*),           intent(in)    :: longname !< The long description of this axis
  character(len=*), optional, intent(in)    :: cartesian !< A variable indicating which direction
                                                     !! this axis corresponds with. Valid values
                                                     !! include 'X', 'Y', 'Z', 'T', and 'N' for none.
  integer,          optional, intent(in)    :: sense !< This is 1 for axes whose values increase upward, or
                                                     !! -1 if they increase downward.
  type(domain1D),   optional, intent(in)    :: domain !< The domain decomposion for this axis
  real, dimension(:), optional, intent(in)  :: data   !< The coordinate values of the points on this axis
  logical,          optional, intent(in)    :: edge_axis !< If true, this axis marks an edge of the tracer cells
  character(len=*), optional, intent(in)    :: calendar !< The name of the calendar used with a time axis

  logical :: is_x, is_y, is_t  ! If true, this is a domain-decomposed axis in one of the directions.
  integer :: position    ! A flag indicating the axis staggering position.
  integer :: kind_, n_, sense_, has_sense_, rc_ ! TIM axis definition args
  character(len=8) :: cart_str

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "write_metadata_axis: file is not open for writing through TIM.")

  axis%name = trim(name)
  is_x = .false. ; is_y = .false. ; is_t = .false.
  cart_str = " "
  if (present(cartesian)) then
    cart_str = trim(adjustl(cartesian))
    if ((index(cart_str, "X") == 1) .or. (index(cart_str, "x") == 1)) is_x = .true.
    if ((index(cart_str, "Y") == 1) .or. (index(cart_str, "y") == 1)) is_y = .true.
    if ((index(cart_str, "T") == 1) .or. (index(cart_str, "t") == 1)) is_t = .true.
  endif
  position = CENTER
  if (present(edge_axis)) then ; if (edge_axis) then
    if (is_x) position = EAST_FACE
    if (is_y) position = NORTH_FACE
  endif ; endif
  kind_ = 3 ; n_ = 0
  if (is_x) then ; kind_ = 0
  elseif (is_y) then ; kind_ = 1
  elseif (is_t .and. .not.present(data)) then ; kind_ = 2
  else
    if (.not.present(data)) call MOM_err(FATAL, "TIM write_metadata_axis: "//&
                      "a data argument is required to register the axis "//trim(name))
    n_ = size(data)
  endif
  if (is_x .or. is_y) axis%domain_decomposed = .true.
  sense_ = 0 ; has_sense_ = 0
  if (present(sense)) then ; sense_ = sense ; has_sense_ = 1 ; endif
  rc_ = tim_io_def_axis(IO_handle%tim_fh, cstr(name), kind_, tim_stag_code(position), n_, &
                        cstr(units), cstr(longname), cstr(cart_str), sense_, has_sense_)
  if (rc_ /= 0) call MOM_err(FATAL, "TIM: def_axis failed for "//trim(name))
  if (present(data)) then
    allocate(axis%ax_data(size(data))) ; axis%ax_data(:) = data(:)
  endif
end subroutine write_metadata_axis

!> Store information about an output variable in a previously defined fieldtype and write this
!! information to the file indicated by unit.
subroutine write_metadata_field(IO_handle, field, axes, name, units, longname, &
                                pack, standard_name, checksum)
  type(file_type),            intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),            intent(inout) :: field !< The fieldtype where this information is stored
  type(axistype), dimension(:), intent(in)  :: axes  !< Handles for the axis used for this variable
  character(len=*),           intent(in)    :: name  !< The name in the file of this variable
  character(len=*),           intent(in)    :: units !< The units of this variable
  character(len=*),           intent(in)    :: longname !< The long description of this variable
  integer,          optional, intent(in)    :: pack  !< A precision reduction factor with which the
                                                     !! variable.  The default, 1, has no reduction,
                                                     !! but 2 is not uncommon.
  character(len=*), optional, intent(in)    :: standard_name !< The standard (e.g., CMOR) name for this variable
  integer(kind=int64), dimension(:), &
                    optional, intent(in)    :: checksum !< Checksum values that can be used to verify reads.

  ! Local variables
  character(len=256), dimension(size(axes)) :: dim_names ! The names of the dimensions
  character(len=2048) :: joined
  character(len=64) :: cks
  character(len=256) :: sname
  integer :: p_, rc_
  integer :: i, ndims

  ndims = size(axes)
  do i=1,ndims ; dim_names(i) = trim(axes(i)%name) ; enddo

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "write_metadata_field: file is not open for writing through TIM.")

  joined = trim(dim_names(1))
  do i=2,ndims ; joined = trim(joined)//char(10)//trim(dim_names(i)) ; enddo
  cks = " "
  if (present(checksum)) write (cks,'(Z16)') checksum(1)
  sname = " " ; if (present(standard_name)) sname = standard_name
  p_ = 1 ; if (present(pack)) p_ = pack
  rc_ = tim_io_def_var(IO_handle%tim_fh, cstr(name), cstr(joined), cstr(units), &
                       cstr(longname), cstr(sname), p_, cstr(cks))
  if (rc_ /= 0) call MOM_err(FATAL, "TIM: def_var failed for "//trim(name))
  field%name = trim(name)
  field%longname = trim(longname)
  field%units = trim(units)
  field%chksum_read = -1
  field%valid_chksum = .false.

end subroutine write_metadata_field

!> Write a global text attribute to a file.
subroutine write_metadata_global(IO_handle, name, attribute)
  type(file_type),            intent(in)    :: IO_handle !< Handle for a file that is open for writing
  character(len=*),           intent(in)    :: name      !< The name in the file of this global attribute
  character(len=*),           intent(in)    :: attribute !< The value of this attribute

  if (IO_handle%tim_fh < 0) &
    call MOM_err(FATAL, "write_metadata_global: file is not open for writing through TIM.")

  if (tim_io_put_global_att(IO_handle%tim_fh, cstr(name), cstr(attribute)) /= 0) &
    call MOM_err(FATAL, "TIM: global attribute failed: "//trim(name))
end subroutine write_metadata_global

! ---------------------------------------------------------------------------
! TIM PIO I/O path helpers.
! ---------------------------------------------------------------------------

!> Returns the TIM-side handle for this domain's decomposition, registering it
!! on first use (memoized by decomposition signature).
integer function tim_get_domain_handle(MOM_Domain)
  type(MOM_domain_type), intent(in) :: MOM_Domain !< Decomposition to register
  integer :: isc, iec, jsc, jec, sym, k
  integer :: sig(7)
  call mpp_get_compute_domain(MOM_Domain%mpp_domain, isc, iec, jsc, jec)
  sym = 0 ; if (MOM_Domain%symmetric) sym = 1
  sig = (/ MOM_Domain%niglobal, MOM_Domain%njglobal, isc, iec, jsc, jec, sym /)
  do k=1,n_tim_domains
    if (all(tim_dom_sig(:,k) == sig)) then
      tim_get_domain_handle = tim_dom_handle(k) ; return
    endif
  enddo
  if (n_tim_domains >= MAX_TIM_DOMAINS) &
    call MOM_err(FATAL, "tim_get_domain_handle: too many distinct domains")
  n_tim_domains = n_tim_domains + 1
  tim_dom_sig(:,n_tim_domains) = sig
  tim_dom_handle(n_tim_domains) = tim_io_register_domain(sig(1), sig(2), &
      sig(3), sig(4), sig(5), sig(6), sig(7))
  tim_get_domain_handle = tim_dom_handle(n_tim_domains)
end function tim_get_domain_handle

!> Returns the TIM decomposition handle for a raw domain2d (used by the
!! external-field seam, which receives mpp domains directly). External
!! fields are cell-centered and on-grid, so the symmetric flag is moot and
!! recorded as false; the center window is identical either way.
integer function tim_get_domain2d_handle(mpp_domain)
  type(domain2d), intent(in) :: mpp_domain !< Decomposition to register
  integer :: isc, iec, jsc, jec, isg, ieg, jsg, jeg, k
  integer :: sig(7)
  call mpp_get_compute_domain(mpp_domain, isc, iec, jsc, jec)
  call mpp_get_global_domain(mpp_domain, isg, ieg, jsg, jeg)
  sig = (/ ieg-isg+1, jeg-jsg+1, isc, iec, jsc, jec, 0 /)
  do k=1,n_tim_domains
    if (all(tim_dom_sig(:,k) == sig)) then
      tim_get_domain2d_handle = tim_dom_handle(k) ; return
    endif
  enddo
  if (n_tim_domains >= MAX_TIM_DOMAINS) &
    call MOM_err(FATAL, "tim_get_domain2d_handle: too many distinct domains")
  n_tim_domains = n_tim_domains + 1
  tim_dom_sig(:,n_tim_domains) = sig
  tim_dom_handle(n_tim_domains) = tim_io_register_domain(sig(1), sig(2), &
      sig(3), sig(4), sig(5), sig(6), sig(7))
  tim_get_domain2d_handle = tim_dom_handle(n_tim_domains)
end function tim_get_domain2d_handle

!> Reads a domain-decomposed 2-d, 3-d or 4-d field via the TIM prototype PIO path.
!! Exactly one of data2d/data3d/data4d must be supplied.
subroutine tim_read_field_dd(filename, fieldname, MOM_Domain, data2d, data3d, data4d, &
                             timelevel, position, scale)
  character(len=*),       intent(in)    :: filename  !< File to read (with or without .nc)
  character(len=*),       intent(in)    :: fieldname !< Variable to read (case-insensitive)
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< Decomposition of the data
  real, dimension(:,:),   optional, intent(inout) :: data2d !< 2-d target array
  real, dimension(:,:,:), optional, intent(inout) :: data3d !< 3-d target array
  real, dimension(:,:,:,:), optional, intent(inout) :: data4d !< 4-d target array
  integer,      optional, intent(in)    :: timelevel !< Record number to read (1-based)
  integer,      optional, intent(in)    :: position  !< Staggering flag (CENTER etc.)
  real,         optional, intent(in)    :: scale     !< Scaling factor applied after read

  real(kind=c_double), allocatable :: buf(:)  ! contiguous compute window, x fastest
  character(len=len(filename)+8) :: fpath
  integer :: handle, stag, tl, rc, nk, nk2
  integer(kind=c_int) :: fsx, fsy  ! file staggering actually present
  integer :: i0, j0                ! first filled window indices
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  integer :: is, ie, js, je, ni, nj, di, dj, i, j, k, k2
  integer :: csz_x, csz_y, dsz_x, dsz_y, sx, sy
  logical :: sym

  fpath = trim(filename)
  if (len_trim(fpath) < 3) call MOM_err(FATAL, "tim_read: bad filename "//trim(filename))
  if (fpath(len_trim(fpath)-2:len_trim(fpath)) /= ".nc") fpath = trim(fpath)//".nc"

  handle = tim_get_domain_handle(MOM_Domain)
  stag = 0
  if (present(position)) then
    if (position == EAST_FACE) stag = 1
    if (position == NORTH_FACE) stag = 2
    if (position == CORNER) stag = 3
  endif
  sym = MOM_Domain%symmetric
  sx = 0 ; if (sym .and. (stag==1 .or. stag==3)) sx = 1
  sy = 0 ; if (sym .and. (stag==2 .or. stag==3)) sy = 1

  ! This rank's read window (must mirror the C++ staggeredExtents rule):
  ! staggered windows extend one point east/north on EVERY rank, matching
  ! mpp_get_compute_domain's position shift; windows overlap at shared edges,
  ! which is legal for reads and fills each rank's full staggered compute window.
  call mpp_get_compute_domain(MOM_Domain%mpp_domain, isc, iec, jsc, jec)
  is = isc ; ie = iec ; if (sx==1) ie = iec+1
  js = jsc ; je = jec ; if (sy==1) je = jec+1
  ni = ie-is+1 ; nj = je-js+1
  nk = 1 ; nk2 = 1
  if (present(data3d)) nk = size(data3d,3)
  if (present(data4d)) then ; nk = size(data4d,3)*size(data4d,4) ; nk2 = size(data4d,4) ; endif

  allocate(buf(ni*nj*nk))
  tl = 0 ; if (present(timelevel)) tl = timelevel
  rc = tim_io_read_decomposed(cstr(fpath), cstr(fieldname), handle, stag, tl, nk, nk2, buf, fsx, fsy)
  if (rc /= 0) call MOM_err(FATAL, "tim_read: failed reading "//trim(fieldname)// &
                            " from "//trim(fpath))

  ! Place the window in the caller's array: offset 0 for compute-sized arrays,
  ! isc-isd for halo (data-domain) arrays; uniform for centered and staggered.
  call mpp_get_data_domain(MOM_Domain%mpp_domain, isd, ied, jsd, jed)
  csz_x = (iec-isc+1) + sx ; dsz_x = (ied-isd+1) + sx
  csz_y = (jec-jsc+1) + sy ; dsz_y = (jed-jsd+1) + sy

  ! Files may lack the symmetric low-edge staggered points (axis size nig
  ! instead of nig+1); those window positions were not filled — skip them, the
  ! model fills the low edges via halo/edge updates (as with FMS).
  i0 = 1 + (sx - fsx) ; j0 = 1 + (sy - fsy)

  if (present(data2d)) then
    di = tim_target_offset(size(data2d,1), csz_x, dsz_x, isc-isd, fieldname)
    dj = tim_target_offset(size(data2d,2), csz_y, dsz_y, jsc-jsd, fieldname)
    do j=j0,nj ; do i=i0,ni
      data2d(di+i, dj+j) = buf(i + (j-1)*ni)
    enddo ; enddo
    if (present(scale)) then ; if (scale /= 1.0) then
      call rescale_comp_data(MOM_Domain, data2d, scale)
    endif ; endif
  elseif (present(data3d)) then
    di = tim_target_offset(size(data3d,1), csz_x, dsz_x, isc-isd, fieldname)
    dj = tim_target_offset(size(data3d,2), csz_y, dsz_y, jsc-jsd, fieldname)
    do k=1,nk ; do j=j0,nj ; do i=i0,ni
      data3d(di+i, dj+j, k) = buf(i + (j-1)*ni + (k-1)*ni*nj)
    enddo ; enddo ; enddo
    if (present(scale)) then ; if (scale /= 1.0) then
      call rescale_comp_data(MOM_Domain, data3d, scale)
    endif ; endif
  else
    di = tim_target_offset(size(data4d,1), csz_x, dsz_x, isc-isd, fieldname)
    dj = tim_target_offset(size(data4d,2), csz_y, dsz_y, jsc-jsd, fieldname)
    do k2=1,nk2 ; do k=1,nk/nk2 ; do j=j0,nj ; do i=i0,ni
      data4d(di+i, dj+j, k, k2) = buf(i + (j-1)*ni + (k-1)*ni*nj + (k2-1)*ni*nj*(nk/nk2))
    enddo ; enddo ; enddo ; enddo
    if (present(scale)) then ; if (scale /= 1.0) then
      call rescale_comp_data(MOM_Domain, data4d, scale)
    endif ; endif
  endif
  deallocate(buf)
end subroutine tim_read_field_dd

!> Reads a whole (replicated, non-decomposed) 0-d or 1-d real variable via the
!! TIM prototype PIO path; every rank receives the full data.
subroutine tim_read_plain(filename, fieldname, n, buf, timelevel)
  character(len=*),    intent(in)    :: filename  !< File to read (with or without .nc)
  character(len=*),    intent(in)    :: fieldname !< Variable to read (case-insensitive)
  integer,             intent(in)    :: n         !< Number of values to read
  real(kind=c_double), intent(inout) :: buf(n)    !< The values read
  integer,   optional, intent(in)    :: timelevel !< Record number to read (1-based)

  character(len=len(filename)+8) :: fpath
  integer :: tl, rc

  fpath = trim(filename)
  if (len_trim(fpath) < 3) call MOM_err(FATAL, "tim_read: bad filename "//trim(filename))
  if (fpath(len_trim(fpath)-2:len_trim(fpath)) /= ".nc") fpath = trim(fpath)//".nc"
  tl = 0 ; if (present(timelevel)) tl = timelevel
  rc = tim_io_read_plain(cstr(fpath), cstr(fieldname), tl, n, buf)
  if (rc /= 0) call MOM_err(FATAL, "tim_read_plain: failed reading "//trim(fieldname)// &
                            " from "//trim(fpath))
end subroutine tim_read_plain

!> Returns the index offset into a caller array for the compute window, judging
!! from its extent whether it is compute-domain or data-domain (halo) sized.
integer function tim_target_offset(sz, compute_sz, data_sz, halo_off, fieldname)
  integer, intent(in) :: sz         !< Actual array extent in this direction
  integer, intent(in) :: compute_sz !< Compute-domain extent (incl. stagger)
  integer, intent(in) :: data_sz    !< Data-domain extent (incl. stagger)
  integer, intent(in) :: halo_off   !< Offset if halo-sized (isc-isd or jsc-jsd)
  character(len=*), intent(in) :: fieldname !< For the error message
  if (sz == compute_sz) then
    tim_target_offset = 0
  elseif (sz == data_sz) then
    tim_target_offset = halo_off
  else
    call MOM_err(FATAL, "tim_read: unexpected array extent for "//trim(fieldname))
    tim_target_offset = 0
  endif
end function tim_target_offset

!> Returns the filename with ".nc" appended when missing (TIM path form).
function tim_norm_path(filename) result(fpath)
  character(len=*), intent(in) :: filename !< A file path, with or without .nc
  character(len=:), allocatable :: fpath
  if (len_trim(filename) >= 3) then
    if (filename(len_trim(filename)-2:len_trim(filename)) == ".nc") then
      fpath = trim(filename)
      return
    endif
  endif
  fpath = trim(filename)//".nc"
end function tim_norm_path

!> Replaces a C string terminator (and everything after it) with blanks.
subroutine tim_cstr_to_f(str)
  character(len=*), intent(inout) :: str !< String possibly holding a c_null_char
  integer :: i
  i = index(str, char(0))
  if (i > 0) str(i:) = " "
end subroutine tim_cstr_to_f

!> Maps an mpp position flag to the TIM stagger code (0..3).
integer function tim_stag_code(position)
  integer, intent(in) :: position !< An mpp position flag (CENTER, EAST_FACE, ...)
  tim_stag_code = 0
  if (position == EAST_FACE) tim_stag_code = 1
  if (position == NORTH_FACE) tim_stag_code = 2
  if (position == CORNER) tim_stag_code = 3
end function tim_stag_code

!> Writes a domain-decomposed 2-d, 3-d or 4-d field via the TIM prototype PIO
!! path. The window buffer layout mirrors the read path; the variable's
!! staggering is looked up from what was registered at define time.
subroutine tim_write_field_dd(IO_handle, name, MOM_Domain, data2d, data3d, data4d, tstamp)
  type(file_type),        intent(inout) :: IO_handle !< Open TIM write file
  character(len=*),       intent(in)    :: name      !< Variable name
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< Decomposition of the data
  real, dimension(:,:),   optional, intent(in) :: data2d !< 2-d data
  real, dimension(:,:,:), optional, intent(in) :: data3d !< 3-d data
  real, dimension(:,:,:,:), optional, intent(in) :: data4d !< 4-d data
  real,         optional, intent(in)    :: tstamp    !< Model time of this field

  real(kind=c_double), allocatable :: buf(:)
  real(kind=c_double) :: ts
  integer :: stag, has_ts, rc, nk
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  integer :: is, ie, js, je, ni, nj, di, dj, i, j, k, k2, nk2
  integer :: csz_x, csz_y, dsz_x, dsz_y, sx, sy
  logical :: sym

  stag = tim_io_var_stagger(IO_handle%tim_fh, cstr(name))
  if (stag < 0) call MOM_err(FATAL, "tim_write: unregistered variable "//trim(name))
  sym = MOM_Domain%symmetric
  sx = 0 ; if (sym .and. (stag==1 .or. stag==3)) sx = 1
  sy = 0 ; if (sym .and. (stag==2 .or. stag==3)) sy = 1

  call mpp_get_compute_domain(MOM_Domain%mpp_domain, isc, iec, jsc, jec)
  is = isc ; ie = iec ; if (sx==1) ie = iec+1
  js = jsc ; je = jec ; if (sy==1) je = jec+1
  ni = ie-is+1 ; nj = je-js+1
  nk = 1 ; nk2 = 1
  if (present(data3d)) nk = size(data3d,3)
  if (present(data4d)) then ; nk = size(data4d,3)*size(data4d,4) ; nk2 = size(data4d,4) ; endif

  call mpp_get_data_domain(MOM_Domain%mpp_domain, isd, ied, jsd, jed)
  csz_x = (iec-isc+1) + sx ; dsz_x = (ied-isd+1) + sx
  csz_y = (jec-jsc+1) + sy ; dsz_y = (jed-jsd+1) + sy

  allocate(buf(ni*nj*nk))
  buf(:) = 0.0
  if (present(data2d)) then
    di = tim_target_offset(size(data2d,1), csz_x, dsz_x, isc-isd, name)
    dj = tim_target_offset(size(data2d,2), csz_y, dsz_y, jsc-jsd, name)
    do j=1,nj ; do i=1,ni
      buf(i + (j-1)*ni) = data2d(di+i, dj+j)
    enddo ; enddo
  elseif (present(data3d)) then
    di = tim_target_offset(size(data3d,1), csz_x, dsz_x, isc-isd, name)
    dj = tim_target_offset(size(data3d,2), csz_y, dsz_y, jsc-jsd, name)
    do k=1,nk ; do j=1,nj ; do i=1,ni
      buf(i + (j-1)*ni + (k-1)*ni*nj) = data3d(di+i, dj+j, k)
    enddo ; enddo ; enddo
  else
    di = tim_target_offset(size(data4d,1), csz_x, dsz_x, isc-isd, name)
    dj = tim_target_offset(size(data4d,2), csz_y, dsz_y, jsc-jsd, name)
    do k2=1,nk2 ; do k=1,nk/nk2 ; do j=1,nj ; do i=1,ni
      buf(i + (j-1)*ni + (k-1)*ni*nj + (k2-1)*ni*nj*(nk/nk2)) = data4d(di+i, dj+j, k, k2)
    enddo ; enddo ; enddo ; enddo
  endif

  ts = 0.0 ; has_ts = 0
  if (present(tstamp)) then ; ts = tstamp ; has_ts = 1 ; endif
  rc = tim_io_write_decomposed(IO_handle%tim_fh, cstr(name), buf, ts, has_ts)
  if (rc /= 0) call MOM_err(FATAL, "tim_write: failed writing "//trim(name))
  deallocate(buf)
  IO_handle%num_times = tim_io_file_num_times(IO_handle%tim_fh)
  IO_handle%file_time = tim_io_file_time(IO_handle%tim_fh)
end subroutine tim_write_field_dd

!> Writes a non-decomposed 0-d/1-d variable via the TIM prototype PIO path.
subroutine tim_write_plain_wrap(IO_handle, name, buf, n, tstamp)
  type(file_type),     intent(inout) :: IO_handle !< Open TIM write file
  character(len=*),    intent(in)    :: name      !< Variable name
  integer,             intent(in)    :: n         !< Number of values
  real(kind=c_double), intent(in)    :: buf(n)    !< The values
  real,      optional, intent(in)    :: tstamp    !< Model time of this field

  real(kind=c_double) :: ts
  integer :: has_ts, rc
  ts = 0.0 ; has_ts = 0
  if (present(tstamp)) then ; ts = tstamp ; has_ts = 1 ; endif
  rc = tim_io_write_plain(IO_handle%tim_fh, cstr(name), buf, n, ts, has_ts)
  if (rc /= 0) call MOM_err(FATAL, "tim_write: failed writing "//trim(name))
  IO_handle%num_times = tim_io_file_num_times(IO_handle%tim_fh)
  IO_handle%file_time = tim_io_file_time(IO_handle%tim_fh)
end subroutine tim_write_plain_wrap

end module MOM_io_infra
