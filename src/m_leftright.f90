module m_leftright

use m_parameters

implicit none

private; public :: s_leftright

contains 

subroutine s_leftright()

  do i=0,n+1 
  
   rol(i)=qpm(1,i)
   ror(i)=qpp(1,i+1)
   ul(i)=qpm(2,i)
   ur(i)=qpp(2,i+1)
   al(i)=sqrt(gama*qpm(3,i)/(qsm(1,i)+eps))
   ar(i)=sqrt(gama*qpp(3,i+1)/(qpp(1,i+1)+eps))
   pl(i)=qpm(3,i)
   pr(i)=qpp(3,i+1)
   sl(i)=min(ul(i)-al(i),ur(i)-ar(i))
   sr(i)=max(ul(i)+al(i),ur(i)+ar(i))
   ss(i)=(pr(i)-pl(i)+rol(i)*ul(i)*(sl(i)-ul(i))&
   -ror(i)*ur(i)*(sr(i)-ur(i)))/(rol(i)*(sl(i)-ul(i))-ror(i)*(sr(i)-ur(i)))
   ds(1,i)=0
   ds(2,i)=1.0d0
   ds(3,i)=ss(i)
   plr(i)=0.5d0*(pl(i)+pr(i)+rol(i)*(sl(i)-ul(i))*(ss(i)-ul(i))+&
   ror(i)*(sr(i)-ur(i))*(ss(i)-ur(i)))
 end do

end subroutine 
end module


