program main
    use globalvar
    implicit none 
!**************************************************************************
!
! Global kriging algorithm
!
!**************************************************************************
!
! Local variables
    integer     :: requnit,resunit
    character   :: reqlocfile*30,resultsfile*30
    real        :: vecreqloc(3)
    integer :: nreqLocs
    logical     :: succeed

    requnit = 1300
    resunit = 1400
    reqlocfile  = 'Cooper[ReqLoc].txt'
    resultsfile = 'Cooper[results].txt'

!   Read dataset of measurements {U1_obs, u2_obs, ..., Un_obs} and 
!   corresponding locations {v1,v2,...,vn}
    call readin
    
!   Calculating the average Miu_u of a set of U values
    call AverageMiuU

!   Calculating the coefficients of matrix A
    call coeffformatxA

!   open results file
    open(resunit,file=resultsfile)

!   check "Cooper[ReqLoc].txt"
    inquire(file=reqlocfile,exist=succeed)
    if(.ne. succeed) then
    ! TODO print error msg in diagnostics file
    end if
!   Cooper[ReqLoc].txt exists, then open and read it
    open(requnit,file=reqlocfile)
    read(requnit,*)
    read(requnit,*) nreqLocs
    if(numreqloc .gt. nreqLocs) then
    !TODO print error msg in diagnostics file
        stop
    end if
    nreqLocs = numreqloc

!   write nReq in "Cooper[results].txt"
    write(resunit,*) nreqLocs
    read(requnit,*)
    do i=1,nreqLocs
        read(requnit,*) vecreqloc(1),vecreqloc(2),vecreqloc(3)
    ! compute coefficients of matrix b 
        call coeffformatxb(vecreqloc)
    ! LU decomposition method for solving A * Lambda = b 
        call LUDecomposition
    ! Compute U^pred at location VecReqLoc
        call Upred_comp(UpredReqLoc)
    ! Write ReqLocs and kriging estimator U^pred in Cooper[results].txt
        write(resunit,1401) i,vecreqloc(1),vecreqloc(2),vecreqloc(3),UpredReqLoc
    end do

    close(resunit)
    close(requnit)

!**************************************************************************
1401 format(3x,I3,4(2X,F18.10))

end program main
