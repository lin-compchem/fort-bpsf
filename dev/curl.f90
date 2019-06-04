module curl_interface
use, intrinsic :: iso_c_binding
#include /usr/include/x86_64-linux-gnu/curl/easy.h
implicit none

type, bind(C) :: CURL
end type CURL

type, bind(C) :: CURLcode
end type CURLcode

type, bind(C) :: CURLOPT_URL
end type CURLOPT_URL

interface
 function curl_easy_init()  bind(C, name="curl_easy_init")
 end function curl_easy_init
end interface

interface
 subroutine curl_easy_setopt()  bind(C, name="curl_easy_setopt")
 end subroutine curl_easy_setopt
end interface

interface
 subroutine curl_easy_perform()  bind(C, name="curl_easy_perform")
 end subroutine curl_easy_perform
end interface

interface
 subroutine curl_easy_cleanup()  bind(C, name="curl_easy_cleanup")
 end subroutine curl_easy_cleanup
end interface

end module curl_interface

program call_tfserving
!    #include <curl/curl.h>
!    #include /usr/include/x86_64-linux-gnu/curl/curl.h
    use iso_c_binding
    use curl_interface
!    include '/usr/include/x86_64-linux-gnu/curl/curl.h'
    implicit none
    character(len=*), parameter :: send_str = '{"instances": [1.0, 2.0, 5.0]}'
    integer handle
    integer(c_int) :: handle2 = 10

    print *, send_str
    print *, handle2

    handle = curl_easy_init()
 !   handle = curl_easy_cleanup()
end program call_tfserving
