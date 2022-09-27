program rmd2
    use, intrinsic :: iso_fortran_env, dp=>real64
    implicit none
    real(dp), parameter   :: PI = 4.0_dp * atan(1.0_dp)
    integer               :: n_atoms, n_dims, n_grid_points, n_iters, n_neighs,&
                             n_samples, i, j, k
    real(dp)              :: grid_extent, cutoff, dt, sig, cutoff_width, u, ke, t, uO
    real(dp), allocatable :: cell(:,:), r(:,:), neigh_vecs(:,:),&
                             plane_norms(:,:), inv_proj(:),&
                             x(:), q(:,:), dudr(:,:), p(:), q_mean(:), at(:,:),&
                             vt(:,:), vth(:,:), p_hist(:,:), xO(:), pO(:), duOdr(:,:),&
                             q_meanO(:)

        
    call read_input_file(n_atoms, n_dims, n_grid_points,&
    n_iters, n_samples, grid_extent, cutoff, dt, sig, cutoff_width, cell, r, q)

    call calc_x(x, n_grid_points, grid_extent)    
    call calc_x(xO, n_grid_points, PI)    

    call calc_neigh_vecs(neigh_vecs, n_neighs, cell, cutoff, n_dims) 
    
    call calc_cell_attribs(plane_norms, inv_proj, cell, n_dims)
    
    call wrap(r, cell, plane_norms, inv_proj, n_atoms, n_dims)

    allocate(p(n_grid_points))
    allocate(dudr(n_dims, n_atoms))
    allocate(pO(n_grid_points))
    allocate(duOdr(n_dims, n_atoms))
    allocate(at(n_dims, n_atoms))
    allocate(vt(n_dims, n_atoms))
    allocate(vth(n_dims, n_atoms))
    allocate(q_mean(n_grid_points))
    allocate(q_meanO(n_grid_points))
    allocate(p_hist(n_grid_points,n_iters/10))

    q_mean = mean(q)

    call evaluate_config(u, dudr, p, r, x, q_mean, n_atoms,&
    n_dims, n_grid_points, n_neighs, sig, cutoff_width, neigh_vecs)
    at(:,:) = -dudr(:,:)
    vt = 0.0_dp
    vth = 0.0_dp

    do i=1, n_iters
        vth = vt + 0.5_dp*at*dt
        do j=1, n_atoms
            r(:,j) = r(:,j) + vth(:,j)*dt
            do k=1, n_dims
                t = dot_product(plane_norms(:,k), r(:,j))*inv_proj(k)
                if (t .gt. 1.0_dp) then
                    r(:,j) = r(:,j) - cell(:,k)
                else if (t .lt. 0.0_dp) then
                    r(:,j) = r(:,j) + cell(:,k)
                end if
            end do
        end do
        
        call evaluate_config(u, dudr, p, r, x, q_mean, n_atoms, n_dims,&
        n_grid_points, n_neighs, sig, cutoff_width, neigh_vecs)
        call calc_ke(ke, vth, n_atoms)

        if (modulo(i,10) .eq. 0) then
            open(unit=1, file="p_hist.dat", status="old", access="sequential",&
            form="formatted", position="append", action="write")
            write(1, *) p
            close(1)
        end if

        at = -dudr
        vt = vth + 0.5_dp*at*dt
        write(*,*) i, u+ke, u, ke
    end do


    deallocate(cell)
    deallocate(r)
    deallocate(neigh_vecs)
    deallocate(plane_norms)
    deallocate(inv_proj)
    deallocate(p)
    deallocate(dudr)
    deallocate(q_mean)
    deallocate(p_hist)

contains

    subroutine calc_ke(ke, v, n_atoms)
        real(dp), allocatable, intent(in)   :: v(:,:)
        real(dp), intent(inout)             :: ke
        integer, intent(in)                 :: n_atoms
        integer                             :: i

        ke = 0.0_dp

        do i=1, n_atoms
            ke = ke + 0.5_dp* dot_product(v(:,i),v(:,i))
        end do
    end subroutine calc_ke

    subroutine evaluate_config(u, dudr, p, r, x, q, n_atoms, n_dims,&
    n_grid_points, n_neighs, sig, cutoff_width, neigh_vecs, uO, duOdr, pO, xO, qO)
        real(dp), allocatable, intent(in)   :: r(:,:), x(:), q(:), neigh_vecs(:,:)
        real(dp), allocatable, optional, intent(in)   :: xO(:), qO(:)
        real(dp), allocatable, intent(inout):: dudr(:,:), p(:)
        real(dp), allocatable, optional, intent(inout):: duOdr(:,:), pO(:)
        real(dp), intent(inout)             :: u
        real(dp), optional, intent(inout)   :: uO
        real(dp), intent(in)                :: sig, cutoff_width
        integer, intent(in)                 :: n_atoms, n_dims, n_grid_points, n_neighs

        real(dp), allocatable               :: dpdr(:,:,:), gauss_basis(:),&
                                               dr(:), dndr(:,:), pp(:),&
                                               gauss_basisO(:), ppO(:), dr2(:),&
                                               dujaldr(:), dOdr(:),&
                                               dcombidr(:), dnOdr(:,:),&
                                               dpOdr(:,:,:)
        integer                             :: l, j, i, k, a, b, counter
        real(dp)                            :: djl, n, log_ratio, Ojal, dal, dotprod,&
                                               scalarprod, dja, Ujal, nO, cutoff_func,&
                                               cutoff_func_2, combi_func

        allocate(dpdr(n_grid_points, n_dims, n_atoms))
        allocate(dpOdr(n_grid_points, n_dims, n_atoms))
        allocate(dndr(n_dims, n_atoms))
        allocate(dnOdr(n_dims, n_atoms))
        allocate(gauss_basis(n_grid_points))
        allocate(gauss_basisO(n_grid_points))
        allocate(dcombidr(n_dims))
        allocate(dr(n_dims))
        allocate(dr2(n_dims))
        allocate(dujaldr(n_dims))
        allocate(dOdr(n_dims))
        allocate(pp(n_grid_points))
        allocate(ppO(n_grid_points))
        
        p = 0.0_dp
        pp = 1e-12_dp
        ppO = 1e-12_dp
        n = 0.0_dp
        nO = 0.0_dp
        u = 0.0_dp
        dudr = 0.0_dp
        dpdr = 0.0_dp
        dndr = 0.0_dp
        dpOdr = 0.0_dp
        dnOdr = 0.0_dp
        if (present(uO)) then
            uO = 0.0_dp
        end if
        if (present(duOdr)) then
            duOdr = 0.0_dp
        end if
        if (present(pO)) then
            pO = 0.0_dp
        end if
        counter = 0
        
        do l=1, n_atoms
            do j=1, n_atoms
                do k=1, n_neighs
                    if (.not.((j .eq. l) .and. (k .eq. 1))) then
                        dr = r(:,l)-(r(:,j)+neigh_vecs(:,k))
                        djl = norm2(dr)
                        if (djl .lt. cutoff) then
                            gauss_basis(:) = exp(-0.5_dp*((x(:)-djl)/sig)**2)
                            pp(:) = pp(:) + gauss_basis(:)
                            do i=1, n_dims
                                dpdr(:,i,l) = dpdr(:,i,l) + 2.0_dp * dr(i) *&
                                gauss_basis(:) * (x(:)-djl)/(djl*sig*sig)
                            end do
                        else if (djl .lt. cutoff+cutoff_width) then
                            gauss_basis(:) = exp(-0.5_dp*((x(:)-djl)/sig)**2)
                            cutoff_func = 0.5_dp*(1.0_dp+cos(PI/cutoff_width*&
                                          (djl-cutoff)))
                            pp(:) = pp(:) + gauss_basis(:)*cutoff_func
                            do i=1, n_dims
                                dpdr(:,i,l) = dpdr(:,i,l) + dr(i)/djl * gauss_basis(:)*&
                                (2.0_dp*(x(:)-djl)/(sig*sig)*cutoff_func - &
                                PI/(cutoff_width)*sin(PI/cutoff_width*(djl-cutoff)))
                            end do
                        end if
                        if (present(uO)) then
                            do a=1, n_atoms
                                do b=1, n_neighs
                                    if (.not. ((a .eq. j) .and. (b .eq. k)) .or.&
                                        .not. ((a .eq. l) .and. (b .eq. 1))) then
                                        dr2 = r(:,l)-(r(:,a)+neigh_vecs(:,b))
                                        dal = norm2(dr2)
                                        if ((djl .lt. cutoff) .and. (dal .lt. cutoff)) then
                                            dotprod = dot_product(dr, dr2)
                                            scalarprod = (dal*djl)
                                            dja = norm2(r(:,j)+neigh_vecs(:,k)-&
                                                  r(:,a)-neigh_vecs(:,b))
                                            Ujal = dotprod/scalarprod
                                            Ojal = acos(Ujal)
                                            gauss_basisO(:) = exp(-0.5_dp*((xO(:)-Ojal)&
                                            /sig)**2)
                                            ppO(:) = ppO(:) + gauss_basisO(:)
                                            dujaldr(:) = scalarprod*(dr(:)+dr2(:))+dotprod*&
                                            2.0_dp*(((dja-dal)/djl)*dr(:)+&
                                            ((dja-djl)/dal)*dr2(:))
                                            dodr(:) = -1.0_dp/sqrt(1-ujal*ujal)*dujaldr(:)
                                            do i=1, n_dims
                                                dpOdr(:,i,l) = dpOdr(:,i,l) +&
                                                (x(:)-Ojal)/(sig*sig)*gauss_basis(:)*&
                                                dOdr(i)
                                            end do
                                            counter = counter + 1
                                        else if ((djl .gt. cutoff .and.&
                                                djl .lt. cutoff+cutoff_width) .and.&
                                                (dal .gt. cutoff .and.&
                                                dal .lt. cutoff+cutoff_width)) then
                                            cutoff_func_2 = 0.5_dp*(1.0_dp +&
                                                               cos(PI/cutoff_width*&
                                                               (dal-cutoff)))
                                            combi_func = sqrt(cutoff_func*cutoff_func_2)
                                            dcombidr = 0.0_dp
                                            dotprod = dot_product(dr, dr2)
                                            scalarprod = (dal*djl)
                                            dja = norm2(r(:,j)+neigh_vecs(:,k)-&
                                                  r(:,a)-neigh_vecs(:,b))
                                            Ujal = dotprod/scalarprod
                                            Ojal = acos(Ujal)
                                            gauss_basisO(:) = exp(-0.5_dp*((xO(:)-Ojal)&
                                            /sig)**2)
                                            ppO(:) = ppO(:) + gauss_basisO(:)*combi_func
                                            dujaldr(:) = scalarprod*(dr(:)+dr2(:))+dotprod*&
                                            2.0_dp*(((dja-dal)/djl)*dr(:)+&
                                            ((dja-djl)/dal)*dr2(:))
                                            dodr(:) = -1.0_dp/sqrt(1-ujal*ujal)*dujaldr(:)
                                            do i=1, n_dims
                                                dcombidr(i) = dcombidr(i)+1.0_dp/&
                                                combi_func*(&
                                                cutoff_func*(-PI/cutoff_width*sin(&
                                                PI/cutoff_width*(djl-cutoff))*dr(i)/djl)+&
                                                cutoff_func_2*(-PI/cutoff_width*sin(&
                                                PI/cutoff_width*(dal-cutoff))*dr2(i)/dal))
                                                dpOdr(:,i,l) = dpOdr(:,i,l) +&
                                                (x(:)-Ojal)/(sig*sig)*gauss_basis(:)*&
                                                dOdr(i)*combi_func+dcombidr(i)*&
                                                gauss_basis(:)
                                            end do
                                        end if
                                    end if
                                end do
                            end do
                        end if
                    end if
                end do
            end do
            do j=1, n_dims
                do i=1, n_grid_points-1
                    dndr(j,l) = dndr(j,l) - (x(i+1)-x(i))*&
                    (dpdr(i+1,j,l)+dpdr(i,j,l))
                    if (present(uO)) then
                        dnOdr(j,l) = dnOdr(j,l) - (xO(i+1)-xO(i))*&
                        (dpOdr(i+1,j,l)+dpOdr(i,j,l))
                    end if
                end do
            end do
        end do

        do i=1, n_grid_points-1
            n = n + (x(i+1)-x(i))*(pp(i+1)+pp(i))
            if (present(uO)) then
                nO = nO + (xO(i+1)-xO(i))*(ppO(i+1)+ppO(i))
            end if
        end do
        
        n = 2.0_dp/n
        if (present(uO)) then
            nO = 2.0_dp/nO
        end if

        dndr = 0.5_dp * dndr * n * n
        if (present(uO)) then
            dnOdr = 0.5_dp * dnOdr * nO * nO
        end if

        p = pp * n        
        if (present(uO)) then
            pO = ppO * nO        
        end if

        do l=1, n_atoms
            do j=1, n_dims
                do i=1, n_grid_points
                    dpdr(i,j,l) = n * dpdr(i,j,l) + pp(i)*dndr(j,l)
                    dudr(j,l) = dudr(j,l) + (1.0_dp + log(p(i)/q(i)))*dpdr(i,j,l)
                    if (present(uO)) then
                        dpOdr(i,j,l) = nO * dpOdr(i,j,l) + ppO(i)*dnOdr(j,l)
                        duOdr(j,l) = duOdr(j,l) + (1.0_dp + log(pO(i)/qO(i)))*dpOdr(i,j,l)
                    end if
                end do
            end do
        end do

        do i=1, n_grid_points
            u = u + p(i)*log(p(i)/q(i))
            if (present(uO)) then
                uO = uO + pO(i)*log(pO(i)/qO(i))
            end if
        end do
        
        !write(*,*) counter

        !open(unit=1, file="ppO.dat", status="replace", access="sequential",&
        !form="formatted", action="write")
        !write(1, *) ppO(:)
        !close(1)

        deallocate(dpdr)
        deallocate(dpOdr)
        deallocate(dndr)
        deallocate(dnOdr)
        deallocate(gauss_basis)
        deallocate(gauss_basisO)
        deallocate(dcombidr)
        deallocate(dr)
        deallocate(dr2)
        deallocate(dujaldr)
        deallocate(dOdr)
        deallocate(pp)
        deallocate(ppO)

    end subroutine evaluate_config

    function mean(array)
        real(dp), intent(in)  :: array(:,:)
        real(dp), allocatable :: mean(:)
        integer               :: i, array_shape(2)

        array_shape = shape(array)

        allocate(mean(array_shape(1)))

        do i=1, array_shape(2)
            mean(:) = mean(:) + array(:,i)
        end do  
        
        mean = mean/real(array_shape(2), dp)

    end function mean

    function cross(a, b)
        real(dp),  intent(in) :: a(:), b(:)
        real(dp) :: cross(3)

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

    end function cross

    subroutine calc_x(x, n_grid_points, grid_extent)
        real(dp), allocatable, intent(inout):: x(:)
        real(dp), intent(in)                :: grid_extent
        integer, intent(in)                 :: n_grid_points

        real(dp)                            :: dx
        integer                             :: i        

        allocate(x(n_grid_points))

        dx = grid_extent/real(n_grid_points - 1, dp)    

        do i=1, n_grid_points
            x(i) = 0.0_dp + real(i-1, dp) * dx
        end do

    end subroutine calc_x

    subroutine wrap(r, cell, plane_norms, inv_proj, n_atoms, n_dims)
        real(dp), allocatable, intent(inout):: r(:,:)
        real(dp), allocatable, intent(in)   :: cell(:,:), inv_proj(:), plane_norms(:,:)
        integer, intent(in)                 :: n_atoms, n_dims

        real(dp)                            :: t
        integer                             :: i, j

        do i=1, n_atoms
            do j=1, n_dims
                t = dot_product(plane_norms(:,j), r(:,i))*inv_proj(j)
                if (t .gt. 1.0_dp) then
                    r(:,i) = r(:,i) - cell(:,j)
                else if (t .lt. 0.0_dp) then
                    r(:,i) = r(:,i) + cell(:,j)
                end if
            end do
        end do

    end subroutine wrap

    subroutine calc_cell_attribs(plane_norms, inv_proj, cell, n_dims)
        real(dp), allocatable, intent(in)   :: cell(:,:)
        real(dp), allocatable, intent(inout):: plane_norms(:,:), inv_proj(:)
        integer, intent(in)                 :: n_dims

        integer                             :: i

        allocate(plane_norms(n_dims, n_dims))
        allocate(inv_proj(n_dims))

        do i=1, n_dims
            plane_norms(:,i) = cross(cell(:,modulo(i, n_dims)+1),&
            cell(:,modulo(i+1, n_dims)+1))
            inv_proj(i) = 1.0_dp/dot_product(plane_norms(:,i), cell(:,i))
        end do

    end subroutine calc_cell_attribs

    subroutine calc_neigh_vecs(neigh_vecs, n_neighs, cell, cutoff, n_dims)
        real(dp), allocatable, intent(in)   :: cell(:,:)
        real(dp), intent(in)                :: cutoff
        integer, intent(in)                 :: n_dims
        real(dp), allocatable, intent(inout):: neigh_vecs(:,:)
        integer, intent(inout)              :: n_neighs

        integer, allocatable                :: neigh_degrees(:)
        real(dp), allocatable               :: ext_cell(:,:)
        integer                             :: i, j, k, c, n_ext_cell,&
                                               iend, jend, kend, skipped

        allocate(neigh_degrees(n_dims))

        n_neighs = 1
        n_ext_cell = 0
        do i=1, n_dims
            neigh_degrees(i) = ceiling(norm2(cell(:,i))/cutoff)
            n_neighs = n_neighs * (2 * neigh_degrees(i) + 1)
            n_ext_cell = n_ext_cell + (2 * neigh_degrees(i) + 1)
        end do
        
        allocate(neigh_vecs(n_dims, n_neighs))
        allocate(ext_cell(n_dims, n_ext_cell))        

        neigh_vecs = 0.0_dp
        
        c = 1
        do i=1, n_dims
            do j=1, 2*neigh_degrees(i)+1
                ext_cell(:,c) = cell(:,i) * real(-neigh_degrees(i) + (j-1), dp)
                c = c + 1
            end do
        end do

        skipped = 0
        c = 2
        iend = 2*neigh_degrees(1)+1
        jend = iend + 2*neigh_degrees(2)+1
        kend = jend + 2*neigh_degrees(3)+1
        do i=1, iend 
            do j=iend+1, jend
                do k=jend+1, kend
                    if ((c .ne. (n_neighs+1)/2+1).or.(skipped .eq. 1)) then
                        neigh_vecs(:,c) = ext_cell(:,i) + ext_cell(:,j) + ext_cell(:,k)
                    else
                        skipped = 1
                        c = c - 1
                    end if
                    c = c + 1
                end do
            end do
        end do

        deallocate(neigh_degrees)
    end subroutine calc_neigh_vecs

    subroutine read_input_file(n_atoms, n_dims, n_grid_points,&
    n_iters, n_samples, grid_extent, cutoff, dt, sig, cutoff_width, cell, r, q)
        integer, intent(inout)              :: n_atoms, n_dims, n_grid_points, n_iters,&
                                               n_samples
        real(dp), intent(inout)             :: grid_extent, cutoff, dt, sig,&
                                               cutoff_width
        real(dp), allocatable, intent(inout):: cell(:,:), r(:,:), q(:,:)    
        
        character(len=50)                   :: in_filename, q_filename
        integer                             :: i    

        call get_command_argument(1, in_filename)
        open(unit=1, file=in_filename, status="old", access="sequential",&
        form="formatted", action="read")

        read(1, *) n_atoms
        read(1, *) n_dims
        read(1, *) n_grid_points
        read(1, *) n_iters
        read(1, *) n_samples
        read(1, *) grid_extent
        read(1, *) cutoff
        read(1, *) dt
        read(1, *) sig
        read(1, *) cutoff_width

        allocate(cell(n_dims, n_dims))
        allocate(r(n_dims, n_atoms))

        do i=1, n_dims
            read(1, *) cell(:,i)
        end do

        do i=1, n_atoms
            read(1, *) r(:,i)
        end do

        close(1)

        call get_command_argument(2, q_filename)
        open(unit=1, file=q_filename, status="old", access="sequential",&
        form="formatted", action="read")
        
        allocate(q(n_grid_points, n_samples))

        do i=1, n_samples
            read(1, *) q(:,i)
        end do

        close(1)
        
    end subroutine read_input_file
end program rmd2
