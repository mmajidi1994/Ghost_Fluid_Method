MODULE m_riemann_solvers

use m_parameters
use m_leftright
use m_flux

IMPLICIT NONE

PRIVATE; PUBLIC :: s_lax, s_hll, s_hllc

contains 
!****************************
subroutine s_lax

  call s_flux 

  do i=0,n
    do j=1,3
     flx(j,i) = 0.5d0*(fp(j,i)+fm(j,i)-alpha*(qsp(j,i+1)-qsm(j,i)))
    end do 
  end do       

end subroutine 
!****************************
subroutine s_hll

call s_qpqs 
call s_flux
call s_leftright

  do i=0,n+1
    
   flx(1,i) =
(sr(i)*fl(1,i)-sl(i)*fr(1,i)+sr(i)*sl(i)*(qsp(1,i+1)-qsm(1,i)))/(sr(i)-sl(i)+eps)
   flx(2,i) =
(sr(i)*fl(2,i)-sl(i)*fr(2,i)+sr(i)*sl(i)*(qsp(2,i+1)-qsm(2,i)))/(sr(i)-sl(i)+eps)
   flx(3,i) =
(sr(i)*fl(3,i)-sl(i)*fr(3,i)+sr(i)*sl(i)*(qsp(3,i+1)-qsm(3,i)))/(sr(i)-sl(i)+eps)
   
  end do 

  do i=0,n+1
     if (sl(i).ge.0) then
     flx(1:3,i)=fl(1:3,i)
     else if (sr(i).le.0) then 
     flx(1:3,i)=fr(1:3,i) 
    end if 
  end do

end subroutine
!****************************
subroutine s_hllc

call s_qpqs 
call s_flux
call s_leftright

do j=1,3
  do i=0,n+1
    
   fls(j,i)=(ss(i)*(sl(i)*qsm(j,i)-fl(j,i))+sl(i)*plr(i)*ds(j,i))/(sl(i)-ss(i)) 
   frs(j,i)=(ss(i)*(sr(i)*qsp(j,i+1)-fr(j,i))+sr(i)*plr(i)*ds(j,i))/(sr(i)-ss(i))
   
  end do 
end do 

  do i=0,n+1
     if (sl(i)>=0) then
     flx(1:3,i)=fl(1:3,i)
     else if (sl(i)<=0 .and. ss(i)>=0) then 
     flx(1:3,i)=fls(1:3,i) 
     else if (ss(i)<=0 .and. sr(i)>=0) then 
     flx(1:3,i)=frs(1:3,i) 
     else if (sr(i)<=0) then 
     flx(1:3,i)=fr(1:3,i)
      
    end if 
  end do

end subroutine 
end module  

