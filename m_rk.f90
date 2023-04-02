module m_rk
use m_parameters
use m_qpqs
use m_rhs 

private; public::s_tvdrk3, s_rk4 

contains 

subroutine s_tvdrk3 ! TVDRK3 to march in time

call s_rhs
do i=1,n
qsold(:,i)=qs(:,i)*(rr(i)**gp)
end do
!1st Step

do i=1,n
qsnew(:,i)=(dt)*rhs(:,i)+ qsold(:,i)
  
   qs(:,i)=qsnew(:,i)/(rr(i)**gp)
   qp(1,i)=qs(1,i)
   qp(2,i)=(qs(2,i)/(qs(1,i)+eps))
   qp(3,i)=(gama-1.0d0)*(qs(3,i)-0.5d0*(qs(2,i)**2.0)/(qs(1,i)+eps))
end do  

call s_rhs

!2nd Step
do i=1,n
qsnew(:,i)=
(0.25d0*dt)*rhs(:,i)+0.75d0*qsold(:,i)+0.25d0*qs(:,i)*(rr(i)**gp)
  
   qs(:,i)=qsnew(:,i)/(rr(i)**gp)
   qp(1,i)=qs(1,i)
   qp(2,i)=(qs(2,i)/(qs(1,i)+eps))
   qp(3,i)=(gama-1.0d0)*(qs(3,i)-0.5d0*(qs(2,i)**2.0)/(qs(1,i)+eps))
end do 
call s_rhs 

!3rd Step
do i=1,n
qsnew(:,i)= (2.0d0*dt/3.0d0)*rhs(:,i)+(1.0d0/3.0d0)*qsold(:,i)&
+(2.0d0/3.0d0)*qs(:,i)*(rr(i)**gp)
   
   qs(:,i)=qsnew(:,i)/(rr(i)**gp)
   qp(1,i)=qs(1,i)
   qp(2,i)=(qs(2,i)/(qs(1,i)+eps))
   qp(3,i)=(gama-1.0d0)*(qs(3,i)-0.5d0*(qs(2,i)**2.0)/(qs(1,i)+eps))
end do 

call s_qpqs

end subroutine
!*****************************************
subroutine s_rk4
call s_rhs

qsold(:,:)=qs(:,:)

!1st Step

qsnew1(:,:)= (dt)*rhs(:,:)
qs(:,:)=0.5d0*qsnew1(:,:)+qsold(:,:)

   qp(1,-2:n+3)=qs(1,-2:n+3)
   qp(2,-2:n+3)=(qs(2,-2:n+3)/(qs(1,-2:n+3)+eps))*(rr(-2:n+3)**gp)
   qp(3,-2:n+3)=(gama-1.0d0)*(qs(3,-2:n+3)-0.5d0*(qs(2,-2:n+3)**2.0)/(qs(1,-2:n+3)+eps))

call s_rhs

!2nd Step

qsnew2(:,:)= (dt)*rhs(:,:)
qs(:,:)=0.5d0*qsnew2(:,:)+qsold(:,:)

   qp(1,-2:n+3)=qs(1,-2:n+3)
   qp(2,-2:n+3)=(qs(2,-2:n+3)/(qs(1,-2:n+3)+eps))*(rr(-2:n+3)**gp)
   qp(3,-2:n+3)=(gama-1.0d0)*(qs(3,-2:n+3)-0.5d0*(qs(2,-2:n+3)**2.0)/(qs(1,-2:n+3)+eps))

call s_rhs

!3rd Step

qsnew3(:,:)= (dt)*rhs(:,:)
qs(:,:)=qsnew3(:,:)+qsold(:,:)

   qp(1,-2:n+3)=qs(1,-2:n+3)
   qp(2,-2:n+3)=(qs(2,-2:n+3)/(qs(1,-2:n+3)+eps))*(rr(-2:n+3)**gp)
   qp(3,-2:n+3)=(gama-1.0d0)*(qs(3,-2:n+3)-0.5d0*(qs(2,-2:n+3)**2.0)/(qs(1,-2:n+3)+eps))

call s_rhs

!4th Step

qsnew4(:,:)= (dt)*rhs(:,:)
qs(:,:)=qsold(:,:)+(1.0/6.0)*(qsnew1(:,:)+2.0d0*qsnew2(:,:)+2.0d0*qsnew3(:,:)+qsnew4(:,:))                              

   qp(1,-2:n+3)=qs(1,-2:n+3)
   qp(2,-2:n+3)=(qs(2,-2:n+3)/(qs(1,-2:n+3)+eps))*(rr(-2:n+3)**gp)
   qp(3,-2:n+3)=(gama-1.0d0)*(qs(3,-2:n+3)-0.5d0*(qs(2,-2:n+3)**2.0)/(qs(1,-2:n+3)+eps))

do i=1,n

  rho(i)=qp(1,i)/(rr(i)**gp)
  u(i)=qp(2,i)/(rr(i)**gp)
  p(i)=qp(3,i)/(rr(i)**gp)

end do

   rho(0) = rho(1)
   rho(n+1)=rho(n)
   u(0)=u(1)  
   u(n+1)=u(n)
   p(0)=p(1)
   p(n+1)=p(n)
end subroutine

end module 

