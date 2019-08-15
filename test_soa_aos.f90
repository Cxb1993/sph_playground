Module OM
   Implicit none
   Real(8), Parameter :: zero = 0.0d0, one=1.0d0
   Type t_Position
      Real(8) :: X=zero, Y=zero, Z=zero
   End Type t_Position
   Type t_Velocity
      Real(8) :: u=zero, v=zero, w=zero
   End Type t_Velocity
   Type t_Acceleration
      Real(8) :: a=zero, b=zero, c=zero
   End Type t_Acceleration
   Type t_Force
      Real(8) :: F=zero, G=zero, H=zero
   End Type t_Force

   Type t_SOA
      Real(8), Allocatable :: mass(:)
      Real(8), Allocatable :: X(:), Y(:), Z(:)
      Real(8), Allocatable :: u(:), v(:), w(:)
      Real(8), Allocatable :: a(:), b(:), c(:)
      Real(8), Allocatable :: F(:), G(:), H(:)
   Contains
      Procedure :: Set => Set_SOA
      Procedure :: Release => Release_SOA
   End Type t_SOA

   Type t_AOS
      Real(8) :: mass
      Type(t_Position) :: Position
      Type(t_Velocity) :: Velocity
      Type(t_Acceleration) :: Acceleration
      Type(t_Force) :: Force
   End Type t_AOS

Contains

   Subroutine Set_SOA(obj, N, iError)
      Class(t_SOA), Intent(OUT) :: obj
      Integer, Intent(IN) :: N
      Integer, Intent(OUT) :: iError
      iError = 0
      if (Allocated(obj%mass)) then
         DeAllocate(obj%mass,             &
                    obj%X, obj%Y, obj%Z,  &
                    obj%u, obj%v, obj%w,  &
                    obj%a, obj%b, obj%c,  &
                    obj%F, obj%G, obj%H)
      endif
      Allocate(obj%mass(N),                     &
               obj%X(N), obj%Y(N), obj%Z(N),    &
               obj%u(N), obj%v(N), obj%w(N),    &
               obj%a(N), obj%b(N), obj%c(N),    &
               obj%F(N), obj%G(N), obj%H(N),    &
               stat=iError)

      obj%mass = zero
      obj%X = zero
      obj%Y = zero
      obj%Z = zero
      obj%u = zero
      obj%v = zero
      obj%w = zero
      obj%a = zero
      obj%b = zero
      obj%c = zero
      obj%F = zero
      obj%G = zero
      obj%H = zero

   End Subroutine Set_SOA

   Subroutine Release_SOA(obj)
      Class(t_SOA), Intent(OUT) :: obj
      if (Allocated(obj%mass) ) then
         DeAllocate(obj%mass,             &
                    obj%X, obj%Y, obj%Z,  &
                    obj%u, obj%v, obj%w,  &
                    obj%a, obj%b, obj%c,  &
                    obj%F, obj%G, obj%H)
      endif
   End Subroutine Release_SOA

   Subroutine Set_AOS(obj, N, iError)
      Type(t_AOS), Intent(OUT), Allocatable :: obj(:)
      Integer, Intent(IN)  :: N
      Integer, Intent(OUT) :: iError
      iError = 0
      Allocate(obj(N), stat=iError)
   End Subroutine Set_AOS

   Subroutine Release_AOS(obj)
      Type(t_AOS), Intent(OUT), Allocatable :: obj(:)
   End Subroutine Release_AOS

   Subroutine Set_MultiArray(POS,VEL,ACC,FRS,MSS,N)
     Real(8), Allocatable, Intent(INOUT) :: POS(:,:),VEL(:,:),ACC(:,:),FRS(:,:), MSS(:)
     Integer, Intent(IN) :: N

     Allocate(POS(3,N))
     Allocate(VEL(3,N))
     Allocate(ACC(3,N))
     Allocate(FRS(3,N))
     Allocate(MSS(N))

     Call Random_Number(POS)
     Call Random_Number(VEL)
     Call Random_Number(ACC)
     Call Random_Number(FRS)
     Call Random_Number(MSS)
     MSS = one + MSS
     FRS = one + FRS
   End Subroutine Set_MultiArray

   Subroutine Release_MultiArray(POS,VEL,ACC,FRS,MSS)
     Real(8), Allocatable, Intent(INOUT) :: POS(:,:),VEL(:,:),ACC(:,:),FRS(:,:), MSS(:)

     DeAllocate(POS)
     DeAllocate(VEL)
     DeAllocate(ACC)
     DeAllocate(FRS)
     DeAllocate(MSS)
   End Subroutine Release_MultiArray

   Subroutine Set_SuperArray(VAR,N)
     Real(8), Allocatable, Intent(INOUT) :: VAR(:,:)
     Integer, Intent(IN) :: N

     Allocate(VAR(3*4+1,N))
     Call Random_Number(VAR)
     VAR(10:12,:) = one + VAR(10:12,:)
     VAR(13,:)    = one + VAR(13,:)
   End Subroutine Set_SuperArray

   Subroutine Release_SuperArray(VAR)
     Real(8), Allocatable, Intent(INOUT) :: VAR(:,:)
     DeAllocate(VAR)
   End Subroutine

   Subroutine Get_Indexes(obj,N)
     Integer, Intent(INOUT), Allocatable :: obj(:)
     Integer, Intent(IN) :: N
     Real(8), Allocatable :: t(:)
     Integer :: i

     Allocate(t(N))
     call random_number(t)
     do i = 1,N
       obj(i) = 1 + FLOOR((N-1)*t(i))
     end do
   End Subroutine Get_Indexes
End Module OM

Program Test
   Use OM
   Use omp_lib
   Implicit none
   Type(t_SOA) :: SOA
   Type(t_AOS), Allocatable :: AOS(:)
   Real(8), Allocatable :: POS(:,:),VEL(:,:),ACC(:,:),FRS(:,:), MSS(:), SAR(:,:)

   Integer :: N
   Integer :: i, j, k
   Integer :: iError
   Integer, Allocatable :: Indexes(:)

   Real(8) :: ts, te, dt

   dt = 1.0d-02
   N = 1e4
   print*,"==================================================="
   do i=1,8

      Call Random_Seed()
      N = N*10
      Write(*,"(A,e7.1)") "Array size N=", real(N)
      print*,"==================================================="

      Call SOA%Set(N, iError)
      Call Random_Number(SOA%mass)
      Call Random_Number(SOA%F)
      Call Random_Number(SOA%G)
      Call Random_Number(SOA%H)
      SOA%mass = one + SOA%mass
      SOA%F = one + SOA%F
      SOA%G = one + SOA%G
      SOA%H = one + SOA%H
      ts = omp_get_wtime()
         do j = 1, N
            SOA%a(j) = SOA%F(j) / SOA%mass(j)
            SOA%b(j) = SOA%G(j) / SOA%mass(j)
            SOA%c(j) = SOA%H(j) / SOA%mass(j)
            SOA%u(j) = SOA%u(j) + (SOA%a(j) * dt)
            SOA%v(j) = SOA%v(j) + (SOA%b(j) * dt)
            SOA%w(j) = SOA%w(j) + (SOA%c(j) * dt)
            SOA%X(j) = SOA%X(j) + (SOA%u(j) * dt)
            SOA%Y(j) = SOA%Y(j) + (SOA%v(j) * dt)
            SOA%Z(j) = SOA%Z(j) + (SOA%w(j) * dt)
         end do
      te = omp_get_wtime()
      Write(*,*) "       Sequential  SOA 1loop:", te - ts
      Call SOA%Release()

      Call SOA%Set(N, iError)
      Call Random_Number(SOA%mass)
      Call Random_Number(SOA%F)
      Call Random_Number(SOA%G)
      Call Random_Number(SOA%H)
      SOA%mass = one + SOA%mass
      SOA%F = one + SOA%F
      SOA%G = one + SOA%G
      SOA%H = one + SOA%H
      ts = omp_get_wtime()
         do j = 1, N
            SOA%a(j) = SOA%F(j) / SOA%mass(j)
            SOA%b(j) = SOA%G(j) / SOA%mass(j)
            SOA%c(j) = SOA%H(j) / SOA%mass(j)
         end do

         do j = 1, N
            SOA%u(j) = SOA%u(j) + (SOA%a(j) * dt)
            SOA%v(j) = SOA%v(j) + (SOA%b(j) * dt)
            SOA%w(j) = SOA%w(j) + (SOA%c(j) * dt)
         end do

         do j = 1, N
            SOA%X(j) = SOA%X(j)  +  (SOA%u(j) * dt)
            SOA%Y(j) = SOA%Y(j)  +  (SOA%v(j) * dt)
            SOA%Z(j) = SOA%Z(j)  +  (SOA%w(j) * dt)
         end do
      te = omp_get_wtime()
      Write(*,*) "       Sequential SOA 3loops:", te - ts
      Call SOA%Release()

      Call Set_AOS(AOS, N, iError)
      Call Random_Number(AOS%mass)
      Call Random_Number(AOS%Force%F)
      Call Random_Number(AOS%Force%G)
      Call Random_Number(AOS%Force%H)
      AOS%mass = one + AOS%mass
      AOS%Force%F = one + AOS%Force%F
      AOS%Force%G = one + AOS%Force%G
      AOS%Force%H = one + AOS%Force%H
      ts = omp_get_wtime()
      do j = 1, N
         AOS(j)%Acceleration%a = AOS(j)%Force%F / AOS(j)%mass
         AOS(j)%Acceleration%b = AOS(j)%Force%G / AOS(j)%mass
         AOS(j)%Acceleration%c = AOS(j)%Force%H / AOS(j)%mass

         AOS(j)%Velocity%u = AOS(j)%Velocity%u + (AOS(j)%Acceleration%a * dt)
         AOS(j)%Velocity%v = AOS(j)%Velocity%v + (AOS(j)%Acceleration%b * dt)
         AOS(j)%Velocity%w = AOS(j)%Velocity%w + (AOS(j)%Acceleration%c * dt)

         AOS(j)%Position%X = AOS(j)%Position%X  +  (AOS(j)%Velocity%u * dt)
         AOS(j)%Position%Y = AOS(j)%Position%Y  +  (AOS(j)%Velocity%v * dt)
         AOS(j)%Position%Z = AOS(j)%Position%Z  +  (AOS(j)%Velocity%w * dt)
      end do
      te = omp_get_wtime()
      Write(*,*) "       Sequential        AOS:", te - ts
      Call Release_AOS(AOS)

      Call Set_MultiArray(POS,VEL,ACC,FRS,MSS,N)
      ts = omp_get_wtime()
      do j = 1, N
          ACC(:,j) = FRS(:,j)/MSS(j)
          VEL(:,j) = VEL(:,j) + ACC(:,j)*dt
          POS(:,j) = POS(:,j) + VEL(:,j)*dt
      end do
      te = omp_get_wtime()
      Write(*,*) "       Sequential        NDA:", te - ts
      Call Release_MultiArray(POS,VEL,ACC,FRS,MSS)

      Call Set_SuperArray(SAR,N)
      ts = omp_get_wtime()
      do j = 1, N
          SAR(7:9,j) = SAR(10:12,j)/SAR(13,j)
          SAR(4:6,j) = SAR(4:6,j) + SAR(7:9,j)*dt
          SAR(1:3,j) = SAR(1:3,j) + SAR(4:6,j)*dt
      end do
      te = omp_get_wtime()
      Write(*,*) "       Sequential        SAR:", te - ts
      Call Release_SuperArray(SAR)

      print*,"---------------------------------------------------"

      Call SOA%Set(N, iError)
      Call Random_Number(SOA%mass)
      Call Random_Number(SOA%F)
      Call Random_Number(SOA%G)
      Call Random_Number(SOA%H)
      SOA%mass = one + SOA%mass
      SOA%F = one + SOA%F
      SOA%G = one + SOA%G
      SOA%H = one + SOA%H
      ts = omp_get_wtime()
         !$omp parallel do default(none)&
         !$omp shared(SOA,dt)
         do j = 1, N
            SOA%a(j) = SOA%F(j) / SOA%mass(j)
            SOA%b(j) = SOA%G(j) / SOA%mass(j)
            SOA%c(j) = SOA%H(j) / SOA%mass(j)
            SOA%u(j) = SOA%u(j) + (SOA%a(j) * dt)
            SOA%v(j) = SOA%v(j) + (SOA%b(j) * dt)
            SOA%w(j) = SOA%w(j) + (SOA%c(j) * dt)
            SOA%X(j) = SOA%X(j) + (SOA%u(j) * dt)
            SOA%Y(j) = SOA%Y(j) + (SOA%v(j) * dt)
            SOA%Z(j) = SOA%Z(j) + (SOA%w(j) * dt)
         end do
         !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP Sequential  SOA 1loop:", te - ts
      Call SOA%Release()

      Call SOA%Set(N, iError)
      Call Random_Number(SOA%mass)
      Call Random_Number(SOA%F)
      Call Random_Number(SOA%G)
      Call Random_Number(SOA%H)
      SOA%mass = one + SOA%mass
      SOA%F = one + SOA%F
      SOA%G = one + SOA%G
      SOA%H = one + SOA%H
      ts = omp_get_wtime()
        !$omp parallel do default(none)&
        !$omp shared(SOA,dt)
         do j = 1, N
            SOA%a(j) = SOA%F(j) / SOA%mass(j)
            SOA%b(j) = SOA%G(j) / SOA%mass(j)
            SOA%c(j) = SOA%H(j) / SOA%mass(j)
         end do
         !$omp end parallel do

         !$omp parallel do default(none)&
         !$omp shared(SOA,dt)
         do j = 1, N
            SOA%u(j) = SOA%u(j) + (SOA%a(j) * dt)
            SOA%v(j) = SOA%v(j) + (SOA%b(j) * dt)
            SOA%w(j) = SOA%w(j) + (SOA%c(j) * dt)
         end do
         !$omp end parallel do

         !$omp parallel do default(none)&
         !$omp shared(SOA,dt)
         do j = 1, N
            SOA%X(j) = SOA%X(j)  +  (SOA%u(j) * dt)
            SOA%Y(j) = SOA%Y(j)  +  (SOA%v(j) * dt)
            SOA%Z(j) = SOA%Z(j)  +  (SOA%w(j) * dt)
         end do
         !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP Sequential SOA 3loops:", te - ts
      Call SOA%Release()

      Call Set_AOS(AOS, N, iError)
      Call Random_Number(AOS%mass)
      Call Random_Number(AOS%Force%F)
      Call Random_Number(AOS%Force%G)
      Call Random_Number(AOS%Force%H)
      AOS%mass = one + AOS%mass
      AOS%Force%F = one + AOS%Force%F
      AOS%Force%G = one + AOS%Force%G
      AOS%Force%H = one + AOS%Force%H
      ts = omp_get_wtime()
      !$omp parallel do default(none)&
      !$omp shared(AOS,dt)
      do j = 1, N

         AOS(j)%Acceleration%a = AOS(j)%Force%F / AOS(j)%mass
         AOS(j)%Acceleration%b = AOS(j)%Force%G / AOS(j)%mass
         AOS(j)%Acceleration%c = AOS(j)%Force%H / AOS(j)%mass

         AOS(j)%Velocity%u = AOS(j)%Velocity%u + (AOS(j)%Acceleration%a * dt)
         AOS(j)%Velocity%v = AOS(j)%Velocity%v + (AOS(j)%Acceleration%b * dt)
         AOS(j)%Velocity%w = AOS(j)%Velocity%w + (AOS(j)%Acceleration%c * dt)

         AOS(j)%Position%X = AOS(j)%Position%X  +  (AOS(j)%Velocity%u * dt)
         AOS(j)%Position%Y = AOS(j)%Position%Y  +  (AOS(j)%Velocity%v * dt)
         AOS(j)%Position%Z = AOS(j)%Position%Z  +  (AOS(j)%Velocity%w * dt)

      end do
      !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP Sequential        AOS:", te - ts
      Call Release_AOS(AOS)

      Call Set_MultiArray(POS,VEL,ACC,FRS,MSS,N)
      ts = omp_get_wtime()
      !$omp parallel do default(none)&
      !$omp shared(ACC,FRS,MSS,VEL,POS,dt)
      do j = 1, N
          ACC(:,j) = FRS(:,j)/MSS(j)
          VEL(:,j) = VEL(:,j) + ACC(:,j)*dt
          POS(:,j) = POS(:,j) + VEL(:,j)*dt
      end do
      !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP Sequential        NDA:", te - ts
      Call Release_MultiArray(POS,VEL,ACC,FRS,MSS)

      Call Set_SuperArray(SAR,N)
      ts = omp_get_wtime()
      !$omp parallel do default(none)&
      !$omp shared(SAR,dt)
      do j = 1, N
          SAR(7:9,j) = SAR(10:12,j)/SAR(13,j)
          SAR(4:6,j) = SAR(4:6,j) + SAR(7:9,j)*dt
          SAR(1:3,j) = SAR(1:3,j) + SAR(4:6,j)*dt
      end do
      !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP Sequential        SAR:", te - ts
      Call Release_SuperArray(SAR)

      print*,"---------------------------------------------------"

      Allocate(Indexes(N))
      call Get_Indexes(Indexes,N)

      Call SOA%Set(N, iError)
      Call Random_Number(SOA%mass)
      Call Random_Number(SOA%F)
      Call Random_Number(SOA%G)
      Call Random_Number(SOA%H)
      SOA%mass = one + SOA%mass
      SOA%F = one + SOA%F
      SOA%G = one + SOA%G
      SOA%H = one + SOA%H
      ts = omp_get_wtime()
       do k = 1, N
          j = Indexes(k)
          SOA%a(j) = SOA%F(j) / SOA%mass(j)
          SOA%b(j) = SOA%G(j) / SOA%mass(j)
          SOA%c(j) = SOA%H(j) / SOA%mass(j)
          SOA%u(j) = SOA%u(j) + (SOA%a(j) * dt)
          SOA%v(j) = SOA%v(j) + (SOA%b(j) * dt)
          SOA%w(j) = SOA%w(j) + (SOA%c(j) * dt)
          SOA%X(j) = SOA%X(j) + (SOA%u(j) * dt)
          SOA%Y(j) = SOA%Y(j) + (SOA%v(j) * dt)
          SOA%Z(j) = SOA%Z(j) + (SOA%w(j) * dt)
       end do
      te = omp_get_wtime()
      Write(*,*) "           Random  SOA 1loop:", te - ts
      Call SOA%Release()

      call Get_Indexes(Indexes,N)
      Call SOA%Set(N, iError)
      Call Random_Number(SOA%mass)
      Call Random_Number(SOA%F)
      Call Random_Number(SOA%G)
      Call Random_Number(SOA%H)
      SOA%mass = one + SOA%mass
      SOA%F = one + SOA%F
      SOA%G = one + SOA%G
      SOA%H = one + SOA%H
      ts = omp_get_wtime()
      do k = 1, N
         j = Indexes(k)
          SOA%a(j) = SOA%F(j) / SOA%mass(j)
          SOA%b(j) = SOA%G(j) / SOA%mass(j)
          SOA%c(j) = SOA%H(j) / SOA%mass(j)
       end do
       do k = 1, N
          j = Indexes(k)
          SOA%u(j) = SOA%u(j) + (SOA%a(j) * dt)
          SOA%v(j) = SOA%v(j) + (SOA%b(j) * dt)
          SOA%w(j) = SOA%w(j) + (SOA%c(j) * dt)
       end do
       do k = 1, N
          j = Indexes(k)
          SOA%X(j) = SOA%X(j)  +  (SOA%u(j) * dt)
          SOA%Y(j) = SOA%Y(j)  +  (SOA%v(j) * dt)
          SOA%Z(j) = SOA%Z(j)  +  (SOA%w(j) * dt)
       end do
      te = omp_get_wtime()
      Write(*,*) "           Random SOA 3loops:", te - ts
      Call SOA%Release()

      call Get_Indexes(Indexes,N)
      Call Set_AOS(AOS, N, iError)
      Call Random_Number(AOS%mass)
      Call Random_Number(AOS%Force%F)
      Call Random_Number(AOS%Force%G)
      Call Random_Number(AOS%Force%H)
      AOS%mass = one + AOS%mass
      AOS%Force%F = one + AOS%Force%F
      AOS%Force%G = one + AOS%Force%G
      AOS%Force%H = one + AOS%Force%H
      ts = omp_get_wtime()

      do k = 1, N
         j = Indexes(k)

         AOS(j)%Acceleration%a = AOS(j)%Force%F / AOS(j)%mass
         AOS(j)%Acceleration%b = AOS(j)%Force%G / AOS(j)%mass
         AOS(j)%Acceleration%c = AOS(j)%Force%H / AOS(j)%mass

         AOS(j)%Velocity%u = AOS(j)%Velocity%u + (AOS(j)%Acceleration%a * dt)
         AOS(j)%Velocity%v = AOS(j)%Velocity%v + (AOS(j)%Acceleration%b * dt)
         AOS(j)%Velocity%w = AOS(j)%Velocity%w + (AOS(j)%Acceleration%c * dt)

         AOS(j)%Position%X = AOS(j)%Position%X  +  (AOS(j)%Velocity%u * dt)
         AOS(j)%Position%Y = AOS(j)%Position%Y  +  (AOS(j)%Velocity%v * dt)
         AOS(j)%Position%Z = AOS(j)%Position%Z  +  (AOS(j)%Velocity%w * dt)

      end do
      te = omp_get_wtime()
      Write(*,*) "           Random        AOS:", te - ts
      Call Release_AOS(AOS)

      Call Set_MultiArray(POS,VEL,ACC,FRS,MSS,N)
      call Get_Indexes(Indexes,N)
      ts = omp_get_wtime()
      do k = 1, N
         j = Indexes(k)
          ACC(:,j) = FRS(:,j)/MSS(j)
          VEL(:,j) = VEL(:,j) + ACC(:,j)*dt
          POS(:,j) = POS(:,j) + VEL(:,j)*dt
      end do
      te = omp_get_wtime()
      Write(*,*) "           Random        NDA:", te - ts
      Call Release_MultiArray(POS,VEL,ACC,FRS,MSS)

      Call Set_SuperArray(SAR,N)
      Call Get_Indexes(Indexes,N)
      ts = omp_get_wtime()
      do k = 1, N
         j = Indexes(k)
          SAR(7:9,j) = SAR(10:12,j)/SAR(13,j)
          SAR(4:6,j) = SAR(4:6,j) + SAR(7:9,j)*dt
          SAR(1:3,j) = SAR(1:3,j) + SAR(4:6,j)*dt
      end do
      te = omp_get_wtime()
      Write(*,*) "           Random        SAR:", te - ts
      Call Release_SuperArray(SAR)

      print*,"---------------------------------------------------"

      call Get_Indexes(Indexes,N)
      Call SOA%Set(N, iError)
      Call Random_Number(SOA%mass)
      Call Random_Number(SOA%F)
      Call Random_Number(SOA%G)
      Call Random_Number(SOA%H)
      SOA%mass = one + SOA%mass
      SOA%F = one + SOA%F
      SOA%G = one + SOA%G
      SOA%H = one + SOA%H
      ts = omp_get_wtime()
        !$omp parallel do default(none)&
        !$omp shared(SOA,Indexes,dt)&
        !$omp private(j,k)
         do k = 1, N
            j = Indexes(k)
            SOA%a(j) = SOA%F(j) / SOA%mass(j)
            SOA%b(j) = SOA%G(j) / SOA%mass(j)
            SOA%c(j) = SOA%H(j) / SOA%mass(j)
            SOA%u(j) = SOA%u(j) + (SOA%a(j) * dt)
            SOA%v(j) = SOA%v(j) + (SOA%b(j) * dt)
            SOA%w(j) = SOA%w(j) + (SOA%c(j) * dt)
            SOA%X(j) = SOA%X(j) + (SOA%u(j) * dt)
            SOA%Y(j) = SOA%Y(j) + (SOA%v(j) * dt)
            SOA%Z(j) = SOA%Z(j) + (SOA%w(j) * dt)
         end do
         !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP     Random  SOA 1loop:", te - ts
      Call SOA%Release()

      call Get_Indexes(Indexes,N)
      Call SOA%Set(N, iError)
      Call Random_Number(SOA%mass)
      Call Random_Number(SOA%F)
      Call Random_Number(SOA%G)
      Call Random_Number(SOA%H)
      SOA%mass = one + SOA%mass
      SOA%F = one + SOA%F
      SOA%G = one + SOA%G
      SOA%H = one + SOA%H
      ts = omp_get_wtime()
        !$omp parallel do default(none)&
        !$omp shared(SOA,Indexes)&
        !$omp private(j,k)
        do k = 1, N
           j = Indexes(k)
            SOA%a(j) = SOA%F(j) / SOA%mass(j)
            SOA%b(j) = SOA%G(j) / SOA%mass(j)
            SOA%c(j) = SOA%H(j) / SOA%mass(j)
         end do
         !$omp end parallel do

         !$omp parallel do default(none)&
         !$omp shared(SOA,Indexes,dt)&
         !$omp private(j,k)
         do k = 1, N
            j = Indexes(k)
            SOA%u(j) = SOA%u(j) + (SOA%a(j) * dt)
            SOA%v(j) = SOA%v(j) + (SOA%b(j) * dt)
            SOA%w(j) = SOA%w(j) + (SOA%c(j) * dt)
         end do
         !$omp end parallel do

         !$omp parallel do default(none)&
         !$omp shared(SOA,Indexes,dt)&
         !$omp private(j,k)
         do k = 1, N
            j = Indexes(k)
            SOA%X(j) = SOA%X(j)  +  (SOA%u(j) * dt)
            SOA%Y(j) = SOA%Y(j)  +  (SOA%v(j) * dt)
            SOA%Z(j) = SOA%Z(j)  +  (SOA%w(j) * dt)
         end do
         !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP     Random SOA 3loops:", te - ts
      Call SOA%Release()

      call Get_Indexes(Indexes,N)
      Call Set_AOS(AOS, N, iError)
      Call Random_Number(AOS%mass)
      Call Random_Number(AOS%Force%F)
      Call Random_Number(AOS%Force%G)
      Call Random_Number(AOS%Force%H)
      AOS%mass = one + AOS%mass
      AOS%Force%F = one + AOS%Force%F
      AOS%Force%G = one + AOS%Force%G
      AOS%Force%H = one + AOS%Force%H
      ts = omp_get_wtime()
      !$omp parallel do default(none)&
      !$omp shared(AOS,Indexes,dt)&
      !$omp private(j,k)
      do k = 1, N
         j = Indexes(k)

         AOS(j)%Acceleration%a = AOS(j)%Force%F / AOS(j)%mass
         AOS(j)%Acceleration%b = AOS(j)%Force%G / AOS(j)%mass
         AOS(j)%Acceleration%c = AOS(j)%Force%H / AOS(j)%mass

         AOS(j)%Velocity%u = AOS(j)%Velocity%u + (AOS(j)%Acceleration%a * dt)
         AOS(j)%Velocity%v = AOS(j)%Velocity%v + (AOS(j)%Acceleration%b * dt)
         AOS(j)%Velocity%w = AOS(j)%Velocity%w + (AOS(j)%Acceleration%c * dt)

         AOS(j)%Position%X = AOS(j)%Position%X  +  (AOS(j)%Velocity%u * dt)
         AOS(j)%Position%Y = AOS(j)%Position%Y  +  (AOS(j)%Velocity%v * dt)
         AOS(j)%Position%Z = AOS(j)%Position%Z  +  (AOS(j)%Velocity%w * dt)

      end do
      !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP     Random        AOS:", te - ts
      Call Release_AOS(AOS)

      Call Set_MultiArray(POS,VEL,ACC,FRS,MSS,N)
      call Get_Indexes(Indexes,N)
      ts = omp_get_wtime()
      !$omp parallel do default(none)&
      !$omp shared(Indexes,POS,VEL,ACC,FRS,MSS,N,dt)&
      !$omp private(j,k)
      do k = 1, N
         j = Indexes(k)
          ACC(:,j) = FRS(:,j)/MSS(j)
          VEL(:,j) = VEL(:,j) + ACC(:,j)*dt
          POS(:,j) = POS(:,j) + VEL(:,j)*dt
      end do
      te = omp_get_wtime()
      Write(*,*) "OpenMP     Random        NDA:", te - ts
      Call Release_MultiArray(POS,VEL,ACC,FRS,MSS)

      Call Set_SuperArray(SAR,N)
      Call Get_Indexes(Indexes,N)
      ts = omp_get_wtime()
      !$omp parallel do default(none)&
      !$omp shared(Indexes,SAR,N,dt)&
      !$omp private(j,k)
      do k = 1, N
         j = Indexes(k)
          SAR(7:9,j) = SAR(10:12,j)/SAR(13,j)
          SAR(4:6,j) = SAR(4:6,j) + SAR(7:9,j)*dt
          SAR(1:3,j) = SAR(1:3,j) + SAR(4:6,j)*dt
      end do
      !$omp end parallel do
      te = omp_get_wtime()
      Write(*,*) "OpenMP     Random        SAR:", te - ts
      Call Release_SuperArray(SAR)

      print*,"==================================================="

      DeAllocate(Indexes)
   end do

End Program Test
