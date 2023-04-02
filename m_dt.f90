module m_dt
use m_parameters

private;public::s_dt
contains 
subroutine s_dt()
 
  do i=1,n
  sn(i)=abs(qp(2,i))+sqrt(gama*qp(3,i)/qp(1,i))
  end do
 
  snmax=maxval(sn(1:n))
  dt=cn*dr/snmax
  alpha=dr/dt 

end subroutine
end module 


