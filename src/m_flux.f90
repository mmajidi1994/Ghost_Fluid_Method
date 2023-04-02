module m_flux

use m_parameters
use m_qpqs 

implicit none

private; public :: s_flux

contains
subroutine s_flux()

  call s_qpqs
 
  do i=0,n+1
    
    fr(1,i)=qsp(2,i+1)
    fr(2,i)=
(1.0-0.5*(gama-1))*(qsp(2,i+1)**2.0d0)/(qsp(1,i+1)+eps)+(gama-1.0)*qsp(3,i+1) 
    fr(3,i)=qpp(2,i+1)*(qsp(3,i+1)+qpp(3,i+1))
   
    fl(1,i)=qsm(2,i)
    fl(2,i)=
(1.0-0.5*(gama-1))*(qsm(2,i)**2.0d0)/(qsm(1,i)+eps)+(gama-1.0)*qsm(3,i) 
    fl(3,i)=qpm(2,i)*(qsm(3,i)+qpm(3,i))
 
  end do 

end subroutine 
end module 

