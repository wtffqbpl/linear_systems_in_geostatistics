module globalvar
    implicit none
!**************************************************************************
!
! Function: define some global variables
!
!**************************************************************************
! ����˵����
!   lambda         --- kriging coefficients
!	vecij			--- ��������(ά����[3,N])
!	Uobs			---	U^obs(N)
!	matA			--- eqn.(2) �о���A��ϵ��, [N,N]
!	matb			---	eqn.(2) ���Ҷ���b [N]
!	rhserr			--- RHS error [N]
!	diatancesA		---	Matrix of distances
!	Nlocations		---	�����ܵ���
!	numsample		---	��������Ҫ����������
!	numreqloc		---	��������Ҫ���Ƶĵ���
!	miuU			---	����U^obs��ƽ��ֵ
!	sill,nugget,range	---	eqn.(6)�еļ�������

	real,allocatable    :: lambda(:),vecij(:,:),Uobs(:),matA(:,:), &
						   matb(:),rhserr(:),distancesA(:,:)
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
!           and allocate memories for variables
!
!**************************************************************************
!
! Local variables
    integer     :: inunit,i
    character   :: inputfile*30
    logical     :: succeed

    inunit=1200

! Read control file and some variables
    inputfile = 'control.txt'
    inquire(file=inputfile,exist=succeed)
    if(.not. succeed) then
		write(*,*)"û�п����ļ�control.txt"
		call errmsg(2)
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
    inputfile = 'Copper[Samples].txt'
    inquire(file=trim(inputfile),exist=succeed)
    if(.not. succeed) then
		write(*,*)"û��Copper[Samples].txt�ļ�"
		call errmsg(3)
! TODO  write error msg in diagnoistics file
        stop
    end if
    open(inunit,file=inputfile)
    read(inunit,*)
    read(inunit,*) Nlocations
    if(numsample .gt. Nlocations) then 
		write(*,*)"��Ҫ������������Ŀǰ�����ṩ��������"
		write(*,*)"Ŀǰ���������У�",Nlocations
		write(*,*)"�����ļ���Ҫ����������:",numsample
		call errmsg(4)
        stop
    end if
    Nlocations = numsample

	! write Nlocations in diagnostics file
	call diagnostics(2,0)

    allocate(vecij(3,Nlocations),Uobs(Nlocations),matA(Nlocations,Nlocations), &
             matb(Nlocations),lambda(Nlocations),rhserr(Nlocations), &
			 distancesA(Nlocations,Nlocations))
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
! Function:  LU decomposition matrix A (A=LU)
!         matrix L and U are stored in A
!
!**************************************************************************
!
! Local variables
	integer	:: r,i,k
	

	! �ֽ�A=LU��L��U�����matA��
	do r=1,Nlocations
		!����U�ĵ�r��
		do i=r,Nlocations
			do k=1,r-1
				matA(r,i)=matA(r,i)-matA(r,k)*matA(k,i)
			end do
		end do

		!����L�ĵ�r��
		do i=r+1,Nlocations
			do k=1,r-1
				matA(i,r)=matA(i,r)-matA(i,k)*matA(k,r)
			end do
			matA(i,r)=matA(i,r)/matA(r,r)
		end do
	end do

	return
end subroutine LUDecomposition












subroutine solvelineareq
    use globalvar
    implicit none 
!**************************************************************************
!
! Function: Solve linear equation LU * Lambda = b
!           method ( eqn. (2) )
!
!**************************************************************************
!
! Local variables

! TODO compute Lambda
	integer	:: i,k
	
	!���Ly=b��y����x��
	do i=1,Nlocations
		lambda(i)=matb(i)
		do k=1,i-1
			lambda(i)=lambda(i)-matA(i,k)*lambda(k)
		end do
	end do
	!���Ux=y
	do i=Nlocations,1,-1
		do k=Nlocations,i+1,-1
			lambda(i)=lambda(i)-matA(i,k)*lambda(k)
		end do
		lambda(i)=lambda(i)/matA(i,i)
	end do


	return
end subroutine solvelineareq


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
    real    :: distances,gammal,covaij

!  Compute coefficients of Matrix A with eqn. (3)
    do i=1,Nlocations
        do j=1,Nlocations
			call distancefunc(vecij(:,j),vecij(:,i),distancesA(j,i))
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
    real    :: distances,gammal,covaij

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
    real    :: l,nugget,sill,range,gammaL

! Local variables
    
    if( l .gt. 0.0 ) then
        gammaL = sill * exp(-3.0*l/range)
    else
        gammaL = nugget+sill
    end if
    return 
end subroutine gammalfunc


subroutine rhserrcomp
	use globalvar
	implicit none 
!*************************************************************************
!
! Function: compute RHS error (RHSERR = b - A * Lambda)
!
!*************************************************************************
!
! Local variables
	integer :: i,j

	rhserr = 0.0

	do i=1,Nlocations
		do j=1,Nlocations
			rhserr(i) = rhserr(i) + matA(i,j)*lambda(j)
		end do
		rhserr(i) = matb(i)-rhserr(i)
	end do

	call diagnostics(9,0)

	return 
end subroutine rhserrcomp



subroutine errmsg(err_type)
	implicit none 
!*************************************************************************
!
! Function: print error msg
!
!*************************************************************************
!
! Arguments:
	integer :: err_type
! Local variables
	character :: errmsgfile*30
	integer	  :: errunit

	errmsgfile='Cooper[errors].txt'
	errunit = 2200

	open(errunit,file=errmsgfile)
	select case(err_type)
		case(1)
			write(errunit,*)"������Cooper[reqloc].txt�ļ�"
		case(2)
			write(errunit,*)"û�п����ļ�control.txt"
		case(3)
			write(errunit,*)"û��Copper[Samples].txt�ļ�"
		case(4)
			write(errunit,*)"��Ҫ������������Ŀǰ�����ṩ��������"
	end select

	return 
end subroutine errmsg


subroutine diagnostics(msg_iflag,ireq)
	use globalvar
	implicit none 
!*************************************************************************
!
! Function: write intermiediate calculations to a "diagnostics file"
!
!*************************************************************************
!
! Arguments:
	integer :: msg_iflag,ireq
! Local variables
	integer	:: diagunit,i,j
	character :: diagfile*30
	
	diagunit=1000
	diagfile='Copper[Diagnostics].txt'

	select case(msg_iflag)
		case(0)
			open(diagunit,file=diagfile)
		case(1)
			!open(diagunit,file=diagfile,
			write(diagunit,*)'Diagnostic file'
		case(2)
			write(diagunit,*)Nlocations
		case(3)
			write(diagunit,*)miuU
		case(4)
			write(diagunit,*)
			write(diagunit,*)"Matrix of distances"
		case(5)
			write(diagunit,*)
			write(diagunit,*)"Covariance matrix A"
			do i=1,Nlocations
				write(diagunit,'(100(E14.4))')(matA(i,j),j=1,Nlocations)
			end do
		case(10)
			write(diagunit,*)
			write(diagunit,*)"Matrix of distances"
			do i=1,Nlocations
				write(diagunit,'(100(E14.4))')(distancesA(i,j),j=1,Nlocations)
			end do
		case(6)
			write(diagunit,*)numreqloc
		case(7)
			write(diagunit,*)
			write(diagunit,*)ireq
			write(diagunit,*)'Covariance vector b'
			write(diagunit,'(100(E14.4))')(matb(i),i=1,Nlocations)
		case(8)
			write(diagunit,*)'Kriging coefficients Lambda'
			write(diagunit,'(100(E14.4))')(lambda(i),i=1,Nlocations)
		case(9)
			write(diagunit,*)'RHS error'
			write(diagunit,'(100(E14.4))')(rhserr(i),i=1,Nlocations)
	end select

	return 
end subroutine diagnostics

