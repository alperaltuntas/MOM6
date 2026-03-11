!> Handy functions for manipulating strings
module MOM_string_infra

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public lowercase

contains

!> Return a string in which all uppercase letters have been replaced by
!! their lowercase counterparts.
function lowercase(input_string)
  character(len=*),     intent(in) :: input_string !< The string to modify
  character(len=len(input_string)) :: lowercase !< The modified output string
!   This function returns a string in which all uppercase letters have been
! replaced by their lowercase counterparts.  It is loosely based on the
! lowercase function in mpp_util.F90.
  integer, parameter :: co=iachar('a')-iachar('A') ! case offset
  integer :: k

  lowercase = input_string
  do k=1, len_trim(input_string)
    if (lowercase(k:k) >= 'A' .and. lowercase(k:k) <= 'Z') &
        lowercase(k:k) = achar(ichar(lowercase(k:k))+co)
  enddo
end function lowercase

end module MOM_string_infra
