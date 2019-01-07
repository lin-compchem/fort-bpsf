!-----------------------------------------------------------------------
!    A tool to convert binary training set files to an ASCII format
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2016 Nongnuch Artrith and Alexander Urban
!+
!+ This program is free software: you can redistribute it and/or modify
!+ it under the terms of the GNU General Public License as published by
!+ the Free Software Foundation, either version 3 of the License, or
!+ (at your option) any later version.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!+ General Public License for more details.
!+
!+ You should have received a copy of the GNU General Public License
!+ along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------
! 2014-09-29 Alexander Urban (AU) and Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
program trnset2ASCII

  use trainset, only: open_TrnSet,       &
                      close_TrnSet,      &
                      save_TrnSet_ASCII, &
                      load_TrnSet_ASCII, &
                      TrnSet

  implicit none

  character(len=1024) :: infile, outfile
  logical             :: to_bin
  type(TrnSet)        :: ts

  call initialize(infile, outfile, to_bin)

  if (to_bin) then
     write(*,*) 'Converting ASCII format to binary format.'
     ts = load_TrnSet_ASCII(outfile, file=infile)
     call close_TrnSet(ts)
     write(*,*) 'Done.'
  else
     write(*,*) 'Converting binary format to ASCII text format.'
     ts = open_TrnSet(infile)
     call save_TrnSet_ASCII(ts, file=outfile)
     write(*,*) 'Done.'
  end if

contains

  subroutine initialize(infile, outfile, to_bin)

    implicit none

    character(len=*), intent(out) :: infile, outfile
    logical,          intent(out) :: to_bin

    integer :: iarg, nargs
    character(len=100) :: arg

    nargs = command_argument_count()
    if (nargs < 1) then
       write(0,*) "Error: No input file provided."
       call print_usage()
       call finalize()
       stop
    end if

    infile = ' '
    outfile = ' '
    to_bin = .false.

    iarg = 1
    do while(iarg <= nargs)
       call get_command_argument(iarg, value=arg)
       select case(trim(arg))
       case('--to-binary')
          to_bin = .true.
       case default
          if (len_trim(infile) == 0) then
             infile = trim(arg)
          else if (len_trim(outfile) == 0) then
             outfile = trim(arg)
          else
             write(0,*) 'Error: Unknown argument: ', trim(arg)
             call finalize()
             stop
          end if
       end select
       iarg = iarg + 1
    end do

    if ((len(infile) == 0) .or. (len(outfile) == 0))then
       write(0,*) 'Error: No input file specified.'
       call finalize()
       stop
    end if

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none

  end subroutine finalize

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

  end subroutine print_usage

end program trnset2ASCII
