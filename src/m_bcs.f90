module m_bcs! Bcs of conserved variables

use m_parameters

implicit none

private; public :: s_bcs

contains

subroutine s_bcs()

do j=1,3
!Left
qp(j,-2:0)=qp(j,1)
qpp(j,-2:0)=qpp(j,1)
qpm(j,-2:0)=qpm(j,1)
!Right
qp(j,n+1:n+3)=qp(j,n)
qpp(j,n+1:n+3)=qpp(j,n)
qpm(j,n+1:n+3)=qpm(j,n)
end do 

end subroutine 
end module 

