module m_startup

use m_parameters

implicit none

private; public :: s_initial_condition

contains 

subroutine s_initial_condition()
   
   rho(-2:dp-1) = 1.0d0;
   rho(dp :n+3)= 1.0d0;
   p(-2:dp-1)= 0.4d0;
   p(dp :n+3)= 0.4d0;
   u(-2:dp-1)= -2.0d0;
   u(dp :n +3)= 2.0d0;
   e(-2:n+3)= p(-2:n+3)/(gama-1.0d0)+0.5d0*rho(-2:n+3)*(u(-2:n+3)**2.0)
   
   !conserved variables
   qs(1,-2:n+3)=rho(-2:n+3)
   qs(2,-2:n+3)=rho(-2:n+3)*u(-2:n+3)
   qs(3,-2:n+3)=e(-2:n+3)
   !primitive variable 
   qp(1,-2:n+3)=rho(-2:n+3)
   qp(2,-2:n+3)=u(-2:n+3)
   qp(3,-2:n+3)=p(-2:n+3)

end subroutine

end module

