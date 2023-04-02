module m_rhs
use m_parameters
use m_rimenann_solvers

private;public::s_rhs
contains 
subroutine s_rhs 

 call s_hllc
 
   do i=1,n
    rhs(1,i) =
-(flx(1,i)*((rr(i)+0.5*dr)**gp)-flx(1,i-1)*((rr(i)-0.5*dr)**gp))/dr  
    rhs(2,i) =
-(flx(2,i)*((rr(i)+0.5*dr)**gp)-flx(2,i-1)*((rr(i)-0.5*dr)**gp))/dr & 
    + gp*qp(3,i)*(rr(i)**(gp-1))
    rhs(3,i) =
-(flx(3,i)*((rr(i)+0.5*dr)**gp)-flx(3,i-1)*((rr(i)-0.5*dr)**gp))/dr 
  end do

end subroutine
end module 


