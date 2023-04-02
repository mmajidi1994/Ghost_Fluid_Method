module m_qpqs! conserved variables

use m_parameters
use m_weno

implicit none

private; public :: s_qpqs

contains

subroutine s_qpqs()

   call s_weno5
  
   qsp(1,-2:n+3)=qpp(1,-2:n+3)
   qsp(2,-2:n+3)=qpp(1,-2:n+3)*qpp(2,-2:n+3) 
   qsp(3,-2:n+3)=qpp(3,-2:n+3)/(gama-1.0d0) + &
   0.5d0*(qpp(1,-2:n+3))*(qpp(2,-2:n+3)**2.0) 
  
   qsm(1,-2:n+3)=qpm(1,-2:n+3)
   qsm(2,-2:n+3)=qpm(1,-2:n+3)*qpm(2,-2:n+3)
   qsm(3,-2:n+3)=qpm(3,-2:n+3)/(gama-1.0d0) + &
   0.5d0*(qpm(1,-2:n+3))*(qpm(2,-2:n+3)**2.0)
 
end subroutine

end module

