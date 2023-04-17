program save_gridnc
    use datastruct, only: &
        grid_structure, &
        scalar_field, &
        vectorinterpol_methods, &
        vector_field_cart

    use lmeshpack
    use smeshpack
    use netcdf


    implicit none
    type(grid_structure)    :: mesh

    integer ::ncid,status
    integer ::dimid_nt,dimid_nv,dimid_ne, &
        dimid_nnc,dimid_nnv,dimid_nne, &
        dimid_nno,dimid_two_grf,dimid_max_chdom, &
        dimid_cell_grf,dimid_edge_grf,dimid_vert_grf, &
        dimid_max_stored_decompositions

    integer :: varid_clon,varid_clat,varid_clon_vertices, &
        varid_clat_vertices,varid_vlon,varid_vlat, &
        varid_vlon_vertices,varid_vlat_vertices
    character(len = 100) :: fileg,filen
    character(len = 2) :: glevel

    !--------------------------------------------------
    glevel = 'g2';
    write(fileg,'(2A)') '../../grid/gridSCVT/',glevel
    write(filen,'(4A)') trim(fileg),'/gridSCVT',glevel,'.nc' 


    write(*,'(2A)') 'Reading Grid: ',fileg
    call loadpts(mesh, fileg);
    call calcTopoVar(mesh);

    write(*,'(2A)') 'Saving: ',filen
    !----------------------------------------------------
    ! Opening file
    status = nf90_create(filen, NF90_NETCDF4, ncid)
    


    !----------------------------------------------------
    ! Define Dimensions
    status = nf90_def_dim(ncid, 'cell', mesh%nt, dimid_nt)
    status = nf90_def_dim(ncid, 'vertex', mesh%nv, dimid_nv)
    status = nf90_def_dim(ncid, 'edge', mesh%ne, dimid_ne)
    status = nf90_def_dim(ncid, 'nc', 2, dimid_nnc)
    status = nf90_def_dim(ncid, 'nv', 3, dimid_nnv)
    status = nf90_def_dim(ncid, 'ne', 6, dimid_nne)
    status = nf90_def_dim(ncid, 'no', 4, dimid_nno)
    status = nf90_def_dim(ncid, 'two_grf', 2, dimid_two_grf)
    status = nf90_def_dim(ncid, 'max_chdom', 1, dimid_max_chdom)
    status = nf90_def_dim(ncid, 'cell_grf', 14, dimid_cell_grf)
    status = nf90_def_dim(ncid, 'edge_grf', 24, dimid_edge_grf)
    status = nf90_def_dim(ncid, 'vert_grf', 24, dimid_vert_grf)
    status = nf90_def_dim(ncid, 'max_stored_decompositions', 4, &
        dimid_max_stored_decompositions)


    !----------------------------------------------------
    ! Define Variable
    status = nf90_def_var(ncid, 'clon', NF90_FLOAT, [dimid_nt], varid_clon)
    status = nf90_def_var(ncid, 'clat', NF90_FLOAT, [dimid_nt], varid_clat)
    status = nf90_def_var(ncid, 'clon_vertices', NF90_FLOAT, [dimid_nt,dimid_nnv], &
        varid_clon_vertices)
    status = nf90_def_var(ncid, 'clat_vertices', NF90_FLOAT, [dimid_nt,dimid_nnv], &
        varid_clat_vertices)
    status = nf90_def_var(ncid, 'vlon', NF90_FLOAT, [dimid_nv],varid_vlon)
    status = nf90_def_var(ncid, 'vlat', NF90_FLOAT, [dimid_nv],varid_vlat)
    status = nf90_def_var(ncid, 'vlon_vertices', NF90_FLOAT, [dimid_nv,dimid_nne], &
        varid_vlon_vertices)
    status = nf90_def_var(ncid, 'vlat_vertices', NF90_FLOAT, [dimid_nv,dimid_nne], &
        varid_vlat_vertices)
    status = nf90_def_var(ncid, 'elon', NF90_FLOAT, [dimid_ne],varid_elon)
    status = nf90_def_var(ncid, 'elat', NF90_FLOAT, [dimid_ne],varid_elat)
    status = nf90_def_var(ncid, 'elon_vertices', NF90_FLOAT, [dimid_ne,dimid_nno], &
        varid_elon_vertices)
    status = nf90_def_var(ncid, 'elat_vertices', NF90_FLOAT, [dimid_ne,dimid_nno], &
        varid_elat_vertices)
    status = nf90_def_var(ncid, 'cell_area', NF90_FLOAT, [dimid_nt],varid_cell_area)
    status = nf90_def_var(ncid, 'dual_area', NF90_FLOAT, [dimid_nv],varid_dual_area)
    status = nf90_def_var(ncid, 'phys_cell_id', NF90_FLOAT, [dimid_nt], &
        varid_phys_cell_id)
    status = nf90_def_var(ncid, 'phys_edge_id', NF90_FLOAT, [dimid_ne], &
        varid_phys_edge_id)
    status = nf90_def_var(ncid, 'lon_cell_centre', NF90_FLOAT, [dimid_nt], &
        varid_lon_cell_centre)
    status = nf90_def_var(ncid, 'lat_cell_centre', NF90_FLOAT, [dimid_nt], &
        varid_lat_cell_centre)
    status = nf90_def_var(ncid, 'longitude_vertices', NF90_FLOAT, [dimid_nv], &
        varid_longitude_vertices)
    status = nf90_def_var(ncid, 'latitude_vertices', NF90_FLOAT, [dimid_nv], &
        varid_latitude_vertices)
    status = nf90_def_var(ncid, 'lon_edge_centre', NF90_FLOAT, [dimid_ne], &
        varid_lon_edge_centre)
    status = nf90_def_var(ncid, 'lat_edge_centre', NF90_FLOAT, [dimid_ne], &
        varid_lat_edge_centre)
    status = nf90_def_var(ncid, 'edge_of_cell', NF90_FLOAT, [dimid_nne,dimid_nt], &
        varid_edge_of_cell) !TODO
    status = nf90_def_var(ncid, 'vertex_of_cell', NF90_FLOAT, [dimid_nnv,dimid_nt], &
        varid_edge_of_cell)
    status = nf90_def_var(ncid, 'adjacent_cell_of_edge', NF90_FLOAT, &
        [dimid_nnc,dimid_ne],varid_adjacent_cell_of_edge)
    status = nf90_def_var(ncid, 'edge_vertices', NF90_FLOAT, &
        [dimid_nnc,dimid_ne],varid_adjacent_cell_of_edge)
    status = nf90_def_var(ncid, 'cells_of_vertex', NF90_FLOAT, &
        [dimid_nne,dimid_nv],varid_cells_of_vertex)
    status = nf90_def_var(ncid, 'edges_of_vertex', NF90_FLOAT, &
        [dimid_nne,dimid_nv],varid_edges_of_vertex)
    status = nf90_def_var(ncid, 'vertices_of_vertex', NF90_FLOAT, &
        [dimid_nne,dimid_nv],varid_vertices_of_vertex)

double cell_area_p(cell) ;
double cell_elevation(cell) ;
int cell_sea_land_mask(cell) ;
double dual_area_p(vertex) ;
double edge_length(edge) ;
double edge_cell_distance(nc, edge) ;
double dual_edge_length(edge) ;
double edgequad_area(edge) ;
double edge_elevation(edge) ;
int edge_sea_land_mask(edge) ;
double edge_vert_distance(nc, edge) ;
double zonal_normal_primal_edge(edge) ;
double meridional_normal_primal_edge(edge) ;
double zonal_normal_dual_edge(edge) ;
double meridional_normal_dual_edge(edge) ;
int orientation_of_normal(nv, cell) ;
int cell_index(cell) ;
int parent_cell_index(cell) ;
int parent_cell_type(cell) ;
int neighbor_cell_index(nv, cell) ;
int child_cell_index(no, cell) ;
int child_cell_id(cell) ;
int edge_index(edge) ;
int edge_parent_type(edge) ;
int vertex_index(vertex) ;
int edge_orientation(ne, vertex) ;
int edge_system_orientation(edge) ;
int refin_c_ctrl(cell) ;
int index_c_list(two_grf, cell_grf) ;
int start_idx_c(max_chdom, cell_grf) ;
int end_idx_c(max_chdom, cell_grf) ;
int refin_e_ctrl(edge) ;
int index_e_list(two_grf, edge_grf) ;
int start_idx_e(max_chdom, edge_grf) ;
int end_idx_e(max_chdom, edge_grf) ;
int refin_v_ctrl(vertex) ;
int index_v_list(two_grf, vert_grf) ;
int start_idx_v(max_chdom, vert_grf) ;
int end_idx_v(max_chdom, vert_grf) ;
int parent_edge_index(edge) ;
int child_edge_index(no, edge) ;
int child_edge_id(edge) ;
int parent_vertex_index(vertex) ;
double cartesian_x_vertices(vertex) ;
double cartesian_y_vertices(vertex) ;
double cartesian_z_vertices(vertex) ;
double edge_middle_cartesian_x(edge) ;
double edge_middle_cartesian_y(edge) ;
double edge_middle_cartesian_z(edge) ;
double edge_dual_middle_cartesian_x(edge) ;
double edge_dual_middle_cartesian_y(edge) ;
double edge_dual_middle_cartesian_z(edge) ;
double edge_primal_normal_cartesian_x(edge) ;
double edge_primal_normal_cartesian_y(edge) ;
double edge_primal_normal_cartesian_z(edge) ;
double edge_dual_normal_cartesian_x(edge) ;
double edge_dual_normal_cartesian_y(edge) ;
double edge_dual_normal_cartesian_z(edge) ;
double cell_circumcenter_cartesian_x(cell) ;
double cell_circumcenter_cartesian_y(cell) ;
double cell_circumcenter_cartesian_z(cell) ;



    !----------------------------------------------------
    ! closing file
    status = nf90_close(ncid)
!    call check(status, 'close')

end program save_gridnc
