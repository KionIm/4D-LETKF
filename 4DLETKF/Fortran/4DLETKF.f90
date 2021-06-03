module mt19937
  implicit none
  integer,private :: N, N1, M, MATA, UMASK, LMASK, TMASKB, TMASKC
  parameter( &
             & N = 624, &
             & N1 = 625, &
             & M = 397, &
             & MATA = -1727483681, &
             & UMASK = -2147483647, &
             & LMASK = 2147483647, &
             & TMASKB = -1658038656, &
             & TMASKC = -272236544 &
             & )
  integer,private :: mti = N1, mt(0:N-1), mag01(0:1) = (/0, MATA/)

  contains

  subroutine sgrnd(seed)
    integer,intent(in) :: seed
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
    mt(0) = iand(seed, -1)
    do mti = 1, N - 1
      mt(mti) = iand(69069 * mt(mti - 1), -1)
    end do
  end subroutine sgrnd

  real(8) function grnd()
    integer :: y, kk

    if(mti >= N) then
!                   generate N words at one time
      if(mti == N + 1) then
!                        if sgrnd() has not been called,
        call sgrnd(4357)
!                          a default initial seed is used
      endif

      do kk = 0, N - M - 1
        y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
        mt(kk) = ieor(ieor(mt(kk + M), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      do kk = N - M, N - 2
        y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
        mt(kk) = ieor(ieor(mt(kk + (M - N)), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      y = ior(iand(mt(N - 1), UMASK), iand(mt(0), LMASK))
      mt(N - 1) = ieor(ieor(mt(M - 1), ishft(y, -1)), mag01(iand(y, 1)))
      mti = 0
    endif

    y = mt(mti)
    mti = mti + 1
    y = ieor(y, ishft(y, -11))
    y = ieor(y, iand(ishft(y, 7), TMASKB))
    y = ieor(y, iand(ishft(y, 15), TMASKC))
    y = ieor(y, ishft(y, -18))

    if(y < 0) then
      grnd = (dble(y) + 2.0d0 ** 32) / (2.0d0 ** 32)
    else
      grnd = dble(y) / (2.0d0 ** 32)
    endif
  end function grnd

end module mt19937

module definition
  implicit NONE
  integer(8),parameter::JM = 40
  real(8),parameter::F = 8d0
  real(8),parameter::dt = 0.01d0 !6h
  real(8),dimension(JM)::Var,dVar
  Integer(8),parameter::NumSt = 3
  Integer(8),parameter::ObsNum = 30
  Integer(8),parameter::Member = 24
  Integer(8),parameter::AssimilationTime = 4 !year
  Integer(8),parameter::DAWindow = 1 !day
  Integer(8),parameter::Nota = 60 !print day on Console
  Real(8),parameter::Inflation = 1.8

  Real(8),dimension(Numst,Numst)::Rini,H,H6,R !np.eye(3)/
  Real(8),dimension(Member)::Sum12,dyRObs,wbar
  Real(8),dimension(Member,Member)::Sum13,dyRdy,Eye,Pa_ny,EyePa_ny,EyeD,UT,W_a
  Real(8),dimension(JM,Member)::Var_am,dX,Var_am2
  Real(8),dimension(JM)::Var_t,meanVar_am
  Real(8),dimension(JM,4*DAWindow)::Var_o
  Integer(8),dimension(JM-ObsNum,4*DAWindow)::b_list
  Real(8),dimension(4*DAWindow,JM)::Xbar_l
  Real(8),dimension(4*DAWindow,JM,Member)::dX_f_l
  Integer(8),dimension(NumSt)::Location
  Real(8)::Var_am_sub,trPa
  Real(8),dimension(JM,JM)::Pa
  Real(8),dimension(NumSt,NumSt) :: EyE3

  Real(8),dimension(NumSt)::Var_o_loc,Obs,xbar_lloc,Var_am_sub2,Obs2
  Real(8),dimension(NumSt,Member)::delX_fm_loc
  Real(8),allocatable :: delY_m_loc(:,:),delY_m_locT(:,:),delY_m_locTRInv(:,:)

  integer(8)::LDA=Member
  integer(8),parameter::LWORK=3*Member
  double precision,dimension( Member )::W
  double precision,dimension(LWORK)::WORK
  integer(8) INFO

   
end module definition

PROGRAM model
  use mt19937
  use definition
  IMPLICIT NONE

  Call DataMaker()

  CAll modelIni()

  CALL Assimilation()

END PROGRAM model

!-------------------------------------------------------------------
subroutine Assimilation()
  use definition
  implicit NONE
  integer(8)::i,i2,i3,i4
  integer(8)::mode

  do i = 1,int(AssimilationTime*365/real(DAWindow))
     mode = mod(i,Nota)
     If (mode==0) then 
       print*,i*DAWindow, "day,RMSE=",sqrt(sum((meanVar_am-Var_t)**2)/real(JM))
       !print*,i*DAWindow, "day,RMSE=",sqrt(sum((Var_o(:,4)-Var_t)**2)/real(JM))
     end if 
     Sum12(:) = 0
     Sum13(:,:) = 0
     do i2 = 1,4*DAWindow
        do i3 = 1,5
           do i4 = 1,Member
              Call RK4(Var_am(:,i4),JM,dt,F)
           end do
        end do
        READ(100_8) Var_t
        READ(200_8) Var_o(:,i2)
        READ(300_8) b_list(:,i2)
        do i3 = 1,JM
          Xbar_l(i2,i3)= sum(Var_am(i3,:))/real(Member)
        end do
        do i3 = 1,Member
          dX_f_l(i2,:,i3) = Var_am(:,i3) - Xbar_l(i2,:)
        end do 
     end do
     Location = (/1,2,3/)
     CALL LETKF4Main(Location)
     do i2 = 2,JM-2
       Location(:) = Location(:) + 1
       CALL LETKF4Main(Location)
     end do
     Location = (/JM-1,JM,1_8/)
     CALL LETKF4Main(Location)
     Location = (/JM,1_8,2_8/)
     CALL LETKF4Main(Location)
     Var_am = Var_am2
     do i2 = 1,JM
        meanVar_am(i2) = Sum(Var_am(i2,:))/real(Member)
     end do
     do i2 = 1,Member
        dX(:,i2) = Var_am(:,i2) - meanVar_am(:)
     end do
     Pa = (1/real((Member-1)))*matmul(dX,TRANSPOSE(dX))
     Write(500_8) sqrt(sum((meanVar_am-Var_t)**2)/real(JM))
     trPa = 0
     do i2 = 1,JM
        trPa = trPa + Pa(i2,i2)
     end do
     trPa = sqrt(trPa/real(JM))
     Write(600_8) trPa

  end do

end subroutine Assimilation

subroutine LETKF4Main(Location_)
use definition 
implicit NONE
integer(8)::i,i2,Count
Integer(8),dimension(JM-ObsNum)::b
Integer(8),dimension(NumSt)::Location_



do i = 1,4*DAWindow
   do i2 = 1,NumSt
      Var_o_loc(i2) = Var_o(Location_(i2),i)
      delX_fm_loc(i2,:) = dX_f_l(i,Location_(i2),:)
      xbar_lloc(i2) = xbar_l(i,Location_(i2))
   end do
   R = Rini
   H6 = EyE3
   b = b_list(:,i)
   call LocHVar_oR(H6,Var_o_loc,R,Count,Location_,b)

   !Count= 2
   allocate (delY_m_loc(Count,Member))
   allocate (delY_m_locT(Member,Count))
   allocate (delY_m_locTRInv(Member,Count))
   delY_m_loc = matmul(H6(1:Count,:),delX_fm_loc)  !(2,1) = (2,3)(3,M)
   !call matinv(R(1:Count,1:Count),Count) 
   delY_m_locT = TRANSPOSE(delY_m_loc) !(M,2) = (2,M).T
   delY_m_locTRInv = matmul(delY_m_locT,R(1:Count,1:Count)) !(M,2) = (M,2)(2,2)
   dyRdy = matmul(delY_m_locTRInv,delY_m_loc) !(M,M) = (M,2)(2,M)
   Obs2(:) = 0d0
   Obs2(1:Count) = matmul(H6(1:Count,:),xbar_lloc)
   Obs(1:Count) = Var_o_loc(1:Count) - Obs2(1:Count)!matmul(H6(1:Count,:),xbar_lloc) !(2,1) = (2,3)(3,1)
   dyRObs = matmul(delY_m_locTRInv,Obs(1:Count)) !(M,1)
   Sum13 = Sum13 + dyRdy
   Sum12 = Sum12 + dyRObs

   deallocate (delY_m_loc)
   deallocate (delY_m_locT)
   deallocate (delY_m_locTRInv)
end do

Pa_ny = (Member-1)*Eye/Inflation + Sum13 !(M,M)
call matinv(Pa_ny,Member)
wbar = matmul(Pa_ny,Sum12) !(M)

Pa_ny = (Member-1)*Pa_ny

CALL DSYEV ('V','U',Member,Pa_ny,Member,W,WORK,LWORK,INFO)

EyePa_ny(:,:) = 0d0
do i = 1,Member
   EyePa_ny(i,i) = sqrt(W(i))
end do

EyePa_ny = matmul(Pa_ny,EyePa_ny)
UT = TRANSPOSE(Pa_ny)
W_a = matmul(EyePa_ny,UT)

do i = 1,Member
  Var_am_sub2 = matmul(delX_fm_loc,wbar + W_a(:,i))
  Var_am_sub = Xbar_l(4*DAWindow,Location_(2)) + Var_am_sub2(2)
  Var_am2(Location_(2),i) = Var_am_sub
end do

end subroutine LETKF4Main


subroutine LocHVar_oR(H6_,Var_o_loc_,R_,Count,Location_,b_list_)
  use definition
  implicit none
  Integer(8)::Count,Count2
  Real(8),dimension(NumSt,NumSt)::H6_,R_
  Real(8),dimension(NumSt)::Var_o_loc_
  Integer(8),dimension(NumSt)::Location_
  Integer(8),dimension(JM-ObsNum)::b_list_
  Integer(8)::i,i2

  Count = 3
  do i = 1,NumSt
    Count2 = 1
    do i2 = 1,JM-ObsNum
      if(Location_(i)==b_list_(i2)) then
        Count = Count - 1
        exit
      end if
      if(Count2 == JM-ObsNum) then
        Var_o_loc_(i+Count-3) = Var_o_loc_(i)
        H6_(i+Count-3,:) = H6_(i,:)
        R_(i+Count-3,i+Count-3) = R_(i,i)
      end if
      Count2 = Count2 + 1
    end do
  end do

end subroutine 



!---------------------------------------------------------------------------
subroutine modelIni()
  use mt19937
  use definition
  IMPLICIT NONE
  Integer(8)::i

  call R_locIni()

  OPEN(100_8,FILE='Var_t.dat',FORM='UNFORMATTED')
  OPEN(200_8,FILE='Var_o.dat',FORM='UNFORMATTED')
  OPEN(300_8,FILE='b_list.dat',FORM='UNFORMATTED')
  OPEN(400_8,FILE='Var_m.dat',FORM='UNFORMATTED')
  OPEN(500_8,FILE='Var_aRMSE.dat',FORM='UNFORMATTED')
  OPEN(600_8,FILE='trPa.dat',FORM='UNFORMATTED')

  Eye(:,:) = 0d0
  do i = 1,Member
    Read(400_8) Var_am(:,i) 
    Eye(i,i) = 1d0
  end do

  EyE3(:,:) = 0d0
  do i = 1,NumSt
    EyE3(i,i) = 1d0
  end do

end subroutine modelIni

subroutine GaspariCohn(r,Gas)
  implicit NONE
  real(8)::r,Gas

  if (r<0) then
    r = -r
  end if
  if (abs(r)<=1) then
    Gas = 1 - (1/4d0)*r**5 + (1/2d0)*r**4 + (5/8d0)*r**3 - (5/3d0)*r**2
  elseif (abs(r)<=2) then
    Gas = (1/12d0)*r**5 - (1/2d0)*r**4 + (5/8d0)*r**3 + (5/3d0)*r**2 - 5*r + 4 - (2/3d0)/r
  else 
    Gas = 0
  end if
end subroutine GaspariCohn

subroutine R_locIni()
use definition
implicit NONE
integer(8)::i
real(8)::r_,Gas_

Rini(:,:) = 0d0
do i = 1,NumSt
  r_ = 4d0*i/real(NumSt+1) - 2d0
  call GaspariCohn(r_,Gas_)
  Rini(i,i) = 1/Gas_
end do

end subroutine R_locIni

!----------------------------------------------------------------------------

subroutine lorenz96(Var,dVar,F,JM)
  IMPLICIT NONE
  integer(8),intent(in)::JM
  real(8),intent(in)::F
  real(8),dimension(JM),intent(in)::Var
  real(8),dimension(JM),intent(out)::dVar
  integer(8)::j

  dVar(1) = (Var(2)-Var(JM-1))*Var(JM) - Var(1) + F
  dVar(2) = (Var(3)-Var(JM))*Var(1) - Var(2) + F
  dVar(JM) = (Var(1)-Var(JM-2))*Var(JM-1) - Var(JM) + F

  Do j = 3,JM-1
    dVar(j) = (Var(j+1)-Var(j-2))*Var(j-1) - Var(j)  + F
  end do

end subroutine lorenz96

subroutine RK4(Var,JM,dt,F)
  implicit NONE
  integer(8),intent(in)::JM
  real(8),intent(in)::dt,F
  real(8),dimension(JM),intent(inout)::Var
  real(8),dimension(JM)::k1,k2,k3,k4,k
  
  call lorenz96(Var,k1,F,JM)
  call lorenz96(Var + 0.5d0*k1*dt,k2,F,JM)
  call lorenz96(Var + 0.5d0*k2*dt,k3,F,JM)
  call lorenz96(Var + k3*dt,k4,F,JM)

  k = (k1+2d0*k2+2d0*k3+k4)/6d0
  Var = Var + k*dt
end subroutine RK4

subroutine DataMaker()
  use mt19937
  use definition
  implicit none
  real(8),dimension(JM)::VarIni,Var_o_
  real(8)::pi=4d0*atan(1d0)
  real(8)::rnd
  integer(8)::i,i2

  integer(8)::accu = 2  !accustum(year)
  integer(8)::DataTime = AssimilationTime !DataAmount(year)
  
  integer(8),dimension(JM-ObsNum)::b_list_

  VarIni(:) = F 
  VarIni(1) = VarIni(1) + 0.008

  do i = 1,4*365*accu
    call RK4(VarIni,JM,dt,F)
  end do

  OPEN(100_8,FILE='Var_t.dat',FORM='UNFORMATTED')
  OPEN(200_8,FILE='Var_o.dat',FORM='UNFORMATTED')
  OPEN(300_8,FILE='b_list.dat',FORM='UNFORMATTED')
  OPEN(400_8,FILE='Var_m.dat',FORM='UNFORMATTED')

  do i = 1,4*365*DataTime
     do i2 = 1,5
       call RK4(VarIni,JM,dt,F)
     end do 
     Write(100_8) VarIni
     do i2 = 1,JM
       rnd = sqrt(-2*log(grnd()))*cos(2*pi*grnd())
       Var_o_(i2) = VarIni(i2) + rnd
     end do
     Write(200_8) Var_o_
     call ObsSelect(b_list_,ObsNum,JM)
     Write(300_8) b_list_
  end do

  do i = 1,Member
      do i2 = 1,4*365
        call RK4(VarIni,JM,dt,F)
      end do
      Write(400_8) VarIni
  end do

  close(100_8)
  close(200_8)
  close(300_8)
  close(400_8)
end subroutine DataMaker

subroutine ObsSelect(list,ObsNum,JM)
implicit NONE
integer(8),intent(in)::ObsNum,JM
real :: x
real(8),dimension(JM)::A
Integer(8),dimension(JM-ObsNum)::B
Integer(8),dimension(JM-ObsNum),intent(out)::list
integer(8)::i,j,k

do i = 1,JM
  call random_number(x)
  A(i) = x
end do

do i = 1,JM-ObsNum
  B(i) = maxloc(A,1)
  A(B(i)) = 0d0
end do

do i=1,JM-ObsNum-1
  do j=i,JM-ObsNum
    if( B(i) > B(j)) then
      k = B(i)
      B(i) = B(j)
      B(j) = k
    endif
  enddo
enddo

list = B
end subroutine ObsSelect

subroutine matinv(a,n)
  Implicit none
  real(8) :: a(n,n),c,dum
  integer(8) :: n,i,j,k
  
  do k=1,n
  c=a(k,k)
  a(k,k)=1
  do j=1,n
  a(k,j)=a(k,j)/c
  enddo
  do i=1,n
  if( i/=k ) then
  dum=a(i,k)
  a(i,k)=0
  do j=1,n
  a(i,j)=a(i,j)-dum*a(k,j)
  enddo
  endif
  enddo
  enddo
  endsubroutine matinv
  
  
