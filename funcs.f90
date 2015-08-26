module globalvar
    implicit none
!**************************************************************************
!
! Function: define some global variables
!
!**************************************************************************
    real,allocatable    :: lambda(:),vecij(:,:),Uobs(:),matA(:,:),matb(:)
    integer             :: Nlocations
    real                :: numsample,numreqloc
    real                :: miuU,sill,nugget,range

end module globalvar


subroutine readin
    use globalvar
    implicit none 
!**************************************************************************
!
! Function: read samples values in the file "Cooper[Samples].txt 
!           and allocate variables
!
!**************************************************************************
!
! Local variables
    integer     :: inunit
    character   :: inputfile*30
    logical     :: succeed

    inunit=1200

! Read control file and some variables
    inputfile = 'control.txt'
    inquire(file=inputfile,exist=succeed)
    if(.ne. succeed) then
! TODO  write error msg in diagnoistics file
        stop
    end if
    open(inunit,file=inputfile)
    read(inunit,*)numsample
    read(inunit,*)numreqloc
    read(inunit,*)sill
    read(inunit,*)nugget
    read(inunit,*)range
    close(inunit)

!   Read sample values in the file "Cooper[Samples].txt"
    inputfile = 'Cooper[Samples].txt'
    inquire(file=inputfile,exist=succeed)
    if(.ne. succeed) then
! TODO  write error msg in diagnoistics file
        stop
    end if
    open(inunit,file=inputfile)
    read(inunit,*)
    read(inunit,*) Nlocations
    if(numsample .gt. Nlocations) then 
    !TODO write error msg in the diagnostics file
        stop
    end if
    Nlocations = numsample

    allocate(vecij(3,Nlocations),Uobs(Nlocations),matA(Nlocations,Nlocations),
             matb(Nlocations),lambda(Nlocations))
    read(inunit,*)
    do i=1,Nlocations
        read(inunit,*)vecij(1,i),vecij(2,i),vecij(3,i),Uobs(i)
    end do
    close(inunit)

    return 
end subroutine readin


subroutine AverageMiuU
    use globalvar
    implicit none 
!**************************************************************************
!
! Function: Calculating the average Miu_u of a set of U values
!                   Stage 1 ---- 1.1
!
!**************************************************************************
!
! Local variables
    integer :: i
    real    :: sumUobs

    sumUobs = 0.0

    do i=1,Nlocations
        sumUobs = sumUobs + Uobs(i)
    end do
    miuU = sumUobs / real(Nlocations)

    return 
end subroutine AverageMiuU
        
        
subroutine distancefunc(veci,vecj,distancesl)
    implicit none
!**************************************************************************
!
! Function: Calculating the distance between a pair of points i and j in 3D 
!           space. ( L[V_i, V_j] )
!                   Stage 1  ---- 1.2
!
!**************************************************************************

! Arguments
    real    :: veci(3),vecj(3),distancesl

! Local variables
    

    distancesl = (veci(1)-vecj(1))**2
    distancesl = distancesl+(veci(2)-vecj(2))**2
    distancesl = distancesl+(veci(3)-vecj(3))**2

    distancesl = sqrt(distancesl)

    return 
end subroutine distancefunc


subroutine Upred_comp(Upred)
    use globalvar
    implicit none
!**************************************************************************
!
! Function: Compute U^pred at location v_0 
!
!**************************************************************************
!
! Arguments
    real    :: Upred
! Local variables
    real    :: Upred_temp 
    integer :: i, j

! Compute eqn. (1) for all locations
    Upred_temp = 0.0
    do i=1,Nlocations
        Upred_temp=Upred_temp+(Uobs(i)-miuU)*lambda(i)
    end do
    Upred=miuU+Upred_temp

    return
end subroutine Upred_comp 


subroutine LUDecomposition
    use globalvar
    implicit none 
!**************************************************************************
!
! Function: Solve linear equation A * Lambda = b with LU decomposition
!           method ( eqn. (2) )
!
!**************************************************************************
!
! Local variables

! TODO compute Lambda


end subroutine LUDecomposition


subroutine coeffformatxA
    use globalvar
    implicit none 
!**************************************************************************
!
! Function: compute coefficients of matrix A 
!
!**************************************************************************
!
! Local variables
    integer :: i, j
    real    :: distances,gammal

!  Compute coefficients of Matrix A with eqn. (3)
    do i=1,Nlocations
        do j=1,Nlocations
            call CovarianceComp(vecij(:,j),vecij(:,i),nugget,sill,range,covaij)
            matA(j,i) = covaij
        end do
    end do

    return 
end subroutine coeffformatxA


subroutine coeffformatxb(vec_zero)
    use globalvar
    implicit none 
!**************************************************************************
!
! Function: compute coefficients of matrix b
!
!**************************************************************************
!
! Arguments
    real    :: vec_zero(3)

! Local variables
    integer :: i, j
    real    :: distances,gammal

!  Compute coefficients of Vector b with eqn. (4)
    do i=1,Nlocations
        call CovarianceComp(vec_zero,vecij(:,i),nugget,sill,range,covaij)
        matb(i) = covaij
    end do

    return 
end subroutine coeffformatxb


subroutine CovarianceComp(veci,vecj,nugget,sill,range,covaij)
    implicit none
!**************************************************************************
!
! Function: compute covariance between Vector I and J
!
!**************************************************************************
!
! Arguments
    real    :: veci(3),vecj(3),nugget,sill,range,covaij

! Local variables
    real     :: distancesij

! Compute distances between Veci and Vecj ( L[V_i, V_j] )
    call distancefunc(veci,vecj,distancesij)

! Compute function Gamma[L] with eqn. (6)
    call gammalfunc(distancesij,nugget,sill,range,covaij)

    return
end subroutine CovarianceComp



subroutine gammalfunc(l,nugget,sill,range,gammaL)
    implicit none 
!*************************************************************************
!
! Function: compute Gamma(L) (eqn. (6))
!
!*************************************************************************

! Arguments:
    real    :: l,nugget,sill,range,gammaL,

! Local variables
    
    if( l .gt. 0.0 ) then
        gammaL = sill * exp(-3.0*l/range)
    else
        gammaL = nugget+sill
    end if
    return 
end subroutine gammalfunc
