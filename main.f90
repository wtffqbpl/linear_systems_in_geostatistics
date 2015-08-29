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
    integer     :: requnit,resunit,exunit
    character   :: reqlocfile*30,resultsfile*30,excelfile*30
    real        :: vecreqloc(3),UpredReqLoc
    integer     :: nreqLocs,i
    logical     :: succeed

    requnit = 1300
    resunit = 1400
	exunit  = 1500
    reqlocfile  = 'Copper[ReqLoc].txt'
    resultsfile = 'Copper[results].txt'
	excelfile	= 'results_excel.csv'

! open diagnostic file and write simple title
	call diagnostics(0,0)
	call diagnostics(1,0)

!   Read dataset of measurements {U1_obs, u2_obs, ..., Un_obs} and 
!   corresponding locations {v1,v2,...,vn}
    call readin
    
!   Calculating the average Miu_u of a set of U values
    call AverageMiuU
!   Write anverage Miu_u in diagnostics file
	call diagnostics(3,0)

!   Calculating the coefficients of matrix A
    call coeffformatxA
!	write distances and Covariance matrix A in diagnostics file
	call diagnostics(10,0)
	call diagnostics(5,0)

!	LU decomposition matrix A and store matrix L and U in A
	call  LUDecomposition

!   open results file and results for excel
    open(resunit,file=resultsfile)
	open(exunit,file=excelfile)

!   check "Copper[ReqLoc].txt"
    inquire(file=reqlocfile,exist=succeed)
    if(.not. succeed) then
		write(*,*)"不存在Cooper[reqloc].txt文件"
		call errmsg(1)
    ! TODO print error msg in diagnostics file
		stop
    end if
!   Copper[ReqLoc].txt exists, then open and read it
    open(requnit,file=reqlocfile)
    read(requnit,*)
    read(requnit,*) nreqLocs
    if(numreqloc .gt. nreqLocs) then
    !TODO print error msg in diagnostics file
        stop
    end if
    nreqLocs = numreqloc

!   write nReq in "Copper[results].txt" and excel output file
	write(exunit,*)"Num,X,Y,Z,U^pred"
    write(resunit,*) nreqLocs
!	write nReq in diagnostics file
	call diagnostics(6,0)

    read(requnit,*)
    do i=1,nreqLocs
        read(requnit,*) vecreqloc(1),vecreqloc(2),vecreqloc(3)

!      ex(1.9) 添加Local kriging算法的地方
!		在这里添加子程序，在样本中寻找与所需坐标最相近的坐标点K，
!		利用K点来进行如下的操作即可完成Local kriging算法
    ! compute coefficients of matrix b 
        call coeffformatxb(vecreqloc)
	! write covariance vector b in diagnostics file
		call diagnostics(7,i)
    ! LU decomposition method for solving A * Lambda = b 
        call solvelineareq
	! write Kriging coefficients Lambda in diagnostics file
		call diagnostics(8,0)
    ! Compute U^pred at location VecReqLoc
        call Upred_comp(UpredReqLoc)
    ! Write ReqLocs and kriging estimator U^pred in Copper[results].txt
        write(resunit,1401) i,vecreqloc(1),vecreqloc(2),vecreqloc(3),UpredReqLoc
		write(exunit, 1501) i,vecreqloc(1),vecreqloc(2),vecreqloc(3),UpredReqLoc

	! Compute RHS error and ouput in Copper[diagnostics].txt
		call rhserrcomp
    end do

    close(resunit)
    close(requnit)
	close(exunit)

!**************************************************************************
1401 format(3x,I3,4(2X,F18.10))
1501 format(I3,',',4(F18.10,','))

end program main
