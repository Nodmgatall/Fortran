module math_stuff
    contains
function gauss(constraints, GCM) result(y)
    INTEGER, INTENT(IN) :: constraints
    REAL, DIMENSION(:,:) :: GCM
    REAL, DIMENSION(constraints) :: y
    INTEGER :: counter
    REAL :: lol
    print*,
do i=1 , constraints
        print*, GCM(i,:)
    end do
    print*,

    !! GAUSS STUFF
    do i=1,constraints -1
        do k=i+1, constraints
            counter = i + 1
            do while(GCM(i,i) == 0.0)
                GCM([counter,i],:) = GCM([i,counter],:)
                counter = counter + 1
                print*, "switched"
                if (counter > constraints) then
                    print*, "Called EXIT"
                    call EXIT(0)
                end if
            end do
            print*, "mu;t"
            xmult  = GCM(k,i) / GCM(i,i)
            print*, xmult
            GCM(k,i) = 0
            do j=i+1,constraints + 1
                print*,k, j, i, GCM(k,j)," - ", xmult," * " ,GCM(i,j)
                GCM(k,j) = GCM(k,j) - (xmult * GCM(i,j))
                print*,"= " ,GCM(k,j)
            end do
        end do
    end do 
        
    print*,
    do l = 1, constraints
            print*, GCM(l,:)
    end do
    
    do i=constraints, 1 ,-1
            lol = GCM(i, constraints + 1)
            do j=i+1, constraints
            lol = lol - GCM(i,j) * y(j)
            end do
            print*, "lol" , lol
            y(i) = lol/GCM(i,i)
        print*,
    end do
end function    
end module math_stuff
PROGRAM SIMPLEX_REV
    use math_stuff 

    INTEGER :: variables, constraints
    INTEGER :: counter
    INTEGER :: status
    INTEGER :: current_entry_variable_index
    INTEGER :: current_exit_variable_index
    !! EINGANGSDATEN
    REAL,dimension(:), allocatable :: function_array
    REAL,dimension(:,:), allocatable :: A
    REAL,dimension(:), allocatable :: b
    REAL,dimension(:), allocatable :: c_T
    INTEGER,dimension(:), allocatable :: variable_indices
    !!OTHER STUFF
    REAL, dimension(:,:), allocatable :: yT_x_B
    REAL, dimension(:,:), allocatable :: B_x_d
    REAL, dimension(:), allocatable :: y
    REAL, dimension(:), allocatable :: v_yT_x_ai
    REAL, dimension(:), allocatable :: d
    REAL, dimension(:), allocatable :: t_buffer
    REAL, dimension(:), allocatable :: b_x_d_buffer
    INTEGER :: index
    INTEGER :: iteration
    !!GET INPUT
    print*, "enter number of variables and constraints"
    !read(*,*) variables, constraints
    print *, variables, constraints
    variables = 3
    constraints = 3
    !!SETUP
    !!ALLOCATE EINGANGSDATEN
    allocate(function_array(1:variables), stat = status)
    allocate(A(1:constraints,1:variables + constraints), stat = status)
    allocate(b(constraints), stat = status)
    allocate(c_T(variables + constraints), stat = status)
    allocate(variable_indices(constraints + variables), stat = status)

    !!BUFFER AND ITERATION VARIABLES
    allocate(y(constraints),stat = status)
    allocate(yT_x_B(constraints, constraints + 1), stat = status)
    allocate(B_x_d(constraints, constraints + 1), stat = status)
    allocate(v_yT_x_ai(variables), stat = status)
    allocate(d(constraints),stat = status)
        ! BUFFER FOR FINDING BIGGEST T
    allocate(t_buffer(constraints), stat = status)
    allocate(b_x_d_buffer(constraints), stat = status)


    !START ACTUAL PROGRAM
    print*, "please enter function to maximize"
    !read(*,*) function_array 
    print*, "please enter constraints"
    !read(*,*) A(:,1:variables + 1)

    !!TEST 3, RESULT SHOULD BE:  x1 = 0, x2 = 2.2, x3 = 1.6 z = =0,6

    function_array = [1, -1, 1]
    A(1,:) = [2, -1, 2, 4]
    A(2,:) = [2, -3, 1, -5]
    A(3,:) = [-1, 1, -2, -1]
    
    
    !!TEST 2, RESULT SHOULD BE: x3 = 10000, z = 10000 / WORKS
    !function_array = [100, 10, 1]
    !A(1,:) = [1, 0 ,0, 1]
    !A(2,:) = [20, 1, 0, 100]
    !A(3,:) = [200, 20, 1 , 10000]

    !!TEST 1 RESULT SHOULD BE x1 = 2 x3 = 1 z = 13 / WORKS
    !function_array = [5,4,3]
    !A(:,1) = [2, 4, 3]
    !A(:,2) = [3, 1, 4]
    !A(:,3) = [1, 2, 2]
    !A(:,4) = [5, 11, 8]
    !!SETUP A_B, b and c_T
    b = A(:,variables + 1)

    A(:,variables + 1:) = 0
    do i = 1, constraints
        A(i, i + variables) = 1
    end do

    c_T = 0.0
    c_T(1:variables) = function_array
    do i = 1 , constraints + variables
        variable_indices(i) = i
    end do
    
    print*, "Setup done"
    print*, "A:"
    do i = 1, variables
    print*, A(i,:)
    end do
    print*, "b:",b
    print*, "c_T",c_T
    iteration = 0
    do
    if (iteration == 5) CALL EXIT()
    iteration = iteration + 1
    print*, "------Starting iteration--------" , iteration
    yT_x_B =0
    yT_x_B(:,1:constraints) = transpose(A(:,variables + 1:))
    yT_x_B(:,constraints + 1) = c_T(variables + 1:)
    print*, "yT_x_B"
    do i = 1, constraints
    print*, yT_x_B(i,:)
    end do
    y = gauss(constraints, yT_x_B)
    print*, "y"
    print*, y
    print*, "A culs"
    do i = 1, variables
    print*,  A(:,i)
    v_yT_x_ai(i) = DOT_PRODUCT(y, A(:,i))
    end do
    print*, "v_yT_x_ai"
    print*, v_yT_x_ai
    v_yT_x_ai = c_T(1:variables) - v_yT_x_ai
    print*, v_yT_x_ai
    if(maxval(v_yT_x_ai) < 0) then
        print*,"RESULT: "
        print*, b
        print*, c_T
        print*, variable_indices
        max_value = 0
        do i = 1 , constraints
        if(variable_indices(i + variables) <= variables) then
        max_value = max_value + (function_array(variable_indices(i + variables)) * b(i))
        end if
        end do
        print*, "max = ", max_value 
        CALL EXIT()
    end if
    current_entry_variable_index = maxloc(v_yT_x_ai, INTEGER)
    print*,current_entry_variable_index

    B_x_d(:,1:constraints) = A(:,variables+1:)
    B_x_d(:,constraints + 1) = A(:,current_entry_variable_index) 

    do i = 1 , constraints
    print*, B_x_d(i,:)
    end do

    d = gauss(constraints, B_x_d)
    print*, "d"
    print*, d

    do i = 1 , constraints
        if(d(i) /= 0.0) then
        t_buffer(i) = b(i)/d(i)
        else
            t_buffer(i) = 0
        end if
        end do
    print*, "t_buffer"
    print*, t_buffer
    print*, 
   
    print*, b

    do i = 1 , constraints
        current_exit_variable_index = maxloc(t_buffer, INTEGER)
        b_x_d_buffer = b - t_buffer(current_exit_variable_index) * d
        print*,b_x_d_buffer
        if(minval(b_x_d_buffer) == 0.0) then
            exit
        end if
        t_buffer(current_exit_variable_index) = 0
        end do

    b_x_d_buffer(current_exit_variable_index) = t_buffer(current_exit_variable_index)
    print*, b_x_d_buffer

    print*, current_exit_variable_index
    print*, current_entry_variable_index
    A(:,[current_entry_variable_index, variables + current_exit_variable_index]) &
        = A(:,[variables + current_exit_variable_index, current_entry_variable_index])
   do i = 1, variables
    print*, A(i,:)
    end do
    test = c_T(current_entry_variable_index)
    print*, "test", test
    c_T(current_entry_variable_index) = c_T(current_exit_variable_index + variables)
    c_T(current_exit_variable_index + variables) = test
    test = variable_indices(current_entry_variable_index)
    variable_indices(current_entry_variable_index) = variable_indices(current_exit_variable_index + variables)
    variable_indices(current_exit_variable_index + variables) = test
    print*,
    print*, c_T
    b = b_x_d_buffer
        end do
    print*, "Results:"
    print*,
    print*, variable_indices
    print*, c_T
    !constraints_matrix(1,1) = -1
    !constraints_matrix(1,2) = 1
    !constraints_matrix(1,3) = 2
    !constraints_matrix(1,4) = 3
    !constraints_matrix(2,1) = 1
    !constraints_matrix(2,2) = -2
    !constraints_matrix(2,3) = -3
    !constraints_matrix(2,4) = 0
    !constraints_matrix(3,1) = 0
    !constraints_matrix(3,2) = 2
    !constraints_matrix(3,3) = 1
    !constraints_matrix(3,4) = -8

    
   ! constraints_matrix(1,1) = 2
   ! constraints_matrix(1,2) = 0
   ! constraints_matrix(1,3) = 0
   ! constraints_matrix(1,4) = 5
   ! constraints_matrix(2,1) = 4
   ! constraints_matrix(2,2) = 1
   ! constraints_matrix(2,3) = 0
   ! constraints_matrix(2,4) = 0
   ! constraints_matrix(3,1) = 3
   ! constraints_matrix(3,2) = 0
   ! constraints_matrix(3,3) = 1
   ! constraints_matrix(3,4) = 0
   ! y = gauss(constraints,constraints_matrix)

    print*,"y"
    !print*,y
!! GAUSS STUFF END
!! start deciding which row will be a
    !for all row do y*row
    !do i = 1, variables
    !a_v(1) = y(1)
    !print*, i
    !print*, constraints_matrix(:,i)
    !print*, DOT_PRODUCT(y , constraints_matrix(:,i))
    !end do
    !save result
    !lol = c_a - result
    !take biggest numbe from lol that is bigger 0
    !if there is none we are done

    deallocate(function_array, stat = status)
    deallocate(A, stat = status)
    deallocate(b, stat = status)
    deallocate(c_T, stat = status)
    deallocate(y, stat = status)
    deallocate(yT_x_B, stat = status)
    deallocate(v_yT_x_ai, stat = status)
    deallocate(B_x_d, stat = status)
    deallocate(d, stat = status)
    
    
    
    
    
    
    
    
!!
END PROGRAM SIMPLEX_REV


