void update_pressure(simulation_data_t *simdata, process_s *process) {
  const double dtdx = simdata->params.dt / simdata->params.dx;

  const int numnodesx = NUMNODESX(simdata->pold);
  const int numnodesy = NUMNODESY(simdata->pold);
  const int numnodesz = NUMNODESZ(simdata->pold);
@@ -1119,7 +1158,10 @@ void update_pressure(simulation_data_t *simdata, process_s *process) {
  MPI_Request requestx;
  MPI_Request requesty;
  MPI_Request requestz;


  int size_direction[3*3];
  size_process(process->coords, process->world, size_direction);

  int m = numnodesx - 1;
  for (int p = 0; p < numnodesz; p++) {
    for (int n = 0; n < numnodesy; n++) {
@@ -1135,12 +1177,14 @@ void update_pressure(simulation_data_t *simdata, process_s *process) {
        dvz -= p > 0 ? GETVALUE(simdata->vzold, m, n, p - 1) : 0.0;

        process->px_bdy[0][p*numnodesy+n] = GETVALUE(simdata->pold, m, n, p);

        printf("DEBUG iteration nb = %d",p*numnodesy+n);
        printf("px_bdy size = %d", size_direction[3]*size_direction[6]);
        SETVALUE(simdata->pnew, m, n, p,
                 GETVALUE(simdata->pold, m, n, p) - rhoc2dtdx * (dvx + dvy + dvz));
    }
  }
  MPI_Isend(process->px_bdy[0], numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[LEFT], 0, process->world->cart_comm, &requestx);

  int n = numnodesy - 1;
  for (int p = 0; p < numnodesz; p++) {
    for (int m = 0; m < numnodesx; m++) {
@@ -1162,6 +1206,7 @@ void update_pressure(simulation_data_t *simdata, process_s *process) {
    }
  }
  MPI_Isend(process->py_bdy[0], numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[DOWN], 1, process->world->cart_comm, &requesty);

  int p = numnodesz - 1;
  for (int n = 0; n < numnodesy; n++) {
    for (int m = 0; m < numnodesx; m++) {
@@ -1207,9 +1252,9 @@ void update_pressure(simulation_data_t *simdata, process_s *process) {
      }
    }
  }
  MPI_Recv(process->vx_bdy, numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[LEFT], 3, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->vy_bdy, numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[DOWN], 4, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->vz_bdy, numnodesx*numnodesy, MPI_DOUBLE, process->neighbors[BACKWARD], 5, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->vx_bdy[1], numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[LEFT], 3, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->vy_bdy[1], numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[DOWN], 4, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->vz_bdy[1], numnodesx*numnodesy, MPI_DOUBLE, process->neighbors[BACKWARD], 5, process->world->cart_comm, MPI_STATUS_IGNORE);


  for (int p = 0; p < numnodesz; p++) {
@@ -1222,9 +1267,9 @@ void update_pressure(simulation_data_t *simdata, process_s *process) {
          double dvy = GETVALUE(simdata->vyold, m, 0, p);
          double dvz = GETVALUE(simdata->vzold, m, n, 0);

          dvx -= process->vx_bdy[p*numnodesy+n];
          dvy -= process->vy_bdy[p*numnodesx+m];
          dvz -= process->vz_bdy[n*numnodesx+m];
          dvx -= process->vx_bdy[1][p*numnodesy+n];
          dvy -= process->vy_bdy[1][p*numnodesx+m];
          dvz -= process->vz_bdy[1][n*numnodesx+m];

          if(m == 0){
            SETVALUE(simdata->pnew, 0, n, p,
@@ -1241,7 +1286,7 @@ void update_pressure(simulation_data_t *simdata, process_s *process) {
      }
    }
  }
  

}

void update_velocities(simulation_data_t *simdata, process_s *process) {
@@ -1291,14 +1336,14 @@ void update_velocities(simulation_data_t *simdata, process_s *process) {
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        process->vx_bdy[p*numnodesy+n] = GETVALUE(simdata->vxold, m, n, p);
        process->vx_bdy[0][p*numnodesy+n] = GETVALUE(simdata->vxold, m, n, p);

        SETVALUE(simdata->vxnew, m, n, p, prev_vx - dtdxrho * dpx);
        SETVALUE(simdata->vynew, m, n, p, prev_vy - dtdxrho * dpy);
        SETVALUE(simdata->vznew, m, n, p, prev_vz - dtdxrho * dpz);
    }
  }
  MPI_Isend(process->vx_bdy, numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[RIGHT], 3, process->world->cart_comm, &requestx);
  MPI_Isend(process->vx_bdy[0], numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[RIGHT], 3, process->world->cart_comm, &requestx);

  int n = numnodesy - 1;
  for (int p = 0; p < numnodesz; p++) {
@@ -1319,14 +1364,14 @@ void update_velocities(simulation_data_t *simdata, process_s *process) {
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        process->vy_bdy[p*numnodesx+m] = GETVALUE(simdata->vxold, m, n, p);
        process->vy_bdy[0][p*numnodesx+m] = GETVALUE(simdata->vxold, m, n, p);

        SETVALUE(simdata->vxnew, m, n, p, prev_vx - dtdxrho * dpx);
        SETVALUE(simdata->vynew, m, n, p, prev_vy - dtdxrho * dpy);
        SETVALUE(simdata->vznew, m, n, p, prev_vz - dtdxrho * dpz);
    }
  }
  MPI_Isend(process->vy_bdy, numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[UP], 4, process->world->cart_comm, &requesty);
  MPI_Isend(process->vy_bdy[0], numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[UP], 4, process->world->cart_comm, &requesty);

  int p = numnodesz - 1;
  for (int n = 0; n < numnodesy; n++) {
@@ -1347,14 +1392,14 @@ void update_velocities(simulation_data_t *simdata, process_s *process) {
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        process->vz_bdy[n*numnodesx+m] = GETVALUE(simdata->vxold, m, n, p);
        process->vz_bdy[0][n*numnodesx+m] = GETVALUE(simdata->vxold, m, n, p);

        SETVALUE(simdata->vxnew, m, n, p, prev_vx - dtdxrho * dpx);
        SETVALUE(simdata->vynew, m, n, p, prev_vy - dtdxrho * dpy);
        SETVALUE(simdata->vznew, m, n, p, prev_vz - dtdxrho * dpz);
    }
  }
  MPI_Isend(process->vz_bdy, numnodesy*numnodesx, MPI_DOUBLE, process->neighbors[FORWARD], 5, process->world->cart_comm, &requestz);
  MPI_Isend(process->vz_bdy[0], numnodesy*numnodesx, MPI_DOUBLE, process->neighbors[FORWARD], 5, process->world->cart_comm, &requestz);

  for (int p = 0; p < numnodesz - 1; p++) {
    for (int n = 0; n < numnodesy - 1; n++) {
@@ -1385,6 +1430,7 @@ void update_velocities(simulation_data_t *simdata, process_s *process) {
}






//OLD VERSION
MPI_Request requestx_p;
  MPI_Request requesty_p;
  MPI_Request requestz_p;
  
  MPI_Request requestx_v;
  MPI_Request requesty_v;
  MPI_Request requestz_v;

  MPI_Isend(process->vx_bdy[0], numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[RIGHT], 3, process->world->cart_comm, &requestx_v);
  MPI_Isend(process->vy_bdy[0], numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[UP], 4, process->world->cart_comm, &requesty_v);
  MPI_Isend(process->vz_bdy[0], numnodesy*numnodesx, MPI_DOUBLE, process->neighbors[FORWARD], 5, process->world->cart_comm, &requestz_v);
    
  MPI_Recv(process->vx_bdy[1], numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[LEFT], 3, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->vy_bdy[1], numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[DOWN], 4, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->vz_bdy[1], numnodesx*numnodesy, MPI_DOUBLE, process->neighbors[BACKWARD], 5, process->world->cart_comm, MPI_STATUS_IGNORE);

  int size_direction[3*3];
  size_process(process->coords, process->world, size_direction);

  int m = 0;
  for (int p = 0; p < numnodesz; p++) {
    for (int n = 0; n < numnodesy; n++) {
        double rhoc2dtdx = GETVALUE(simdata->rho, m, n, p) *
                           GETVALUE(simdata->c, m, n, p) *
                           GETVALUE(simdata->c, m, n, p) * dtdx;

        double dvx = GETVALUE(simdata->vxold, m, n, p);
        double dvy = GETVALUE(simdata->vyold, m, n, p);
        double dvz = GETVALUE(simdata->vzold, m, n, p);

        dvx -= process->vx_bdy[1][p*numnodesy+n];
        dvy -= n > 0 ? GETVALUE(simdata->vyold, m, n - 1, p) : process->vy_bdy[1][p*numnodesx + m];
        dvz -= p > 0 ? GETVALUE(simdata->vzold, m, n, p - 1) : process->vz_bdy[1][n*numnodesx + m];

        double value = GETVALUE(simdata->pold, m, n, p) - rhoc2dtdx * (dvx + dvy + dvz);

        process->px_bdy[0][p*numnodesy+n] = value;
        SETVALUE(simdata->pnew, m, n, p, value);
    }
  }
  MPI_Isend(process->px_bdy[0], numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[LEFT], 0, process->world->cart_comm, &requestx_p);

  int n = 0;
  for (int p = 0; p < numnodesz; p++) {
    for (int m = 0; m < numnodesx; m++) {
        double rhoc2dtdx = GETVALUE(simdata->rho, m, n, p) *
                           GETVALUE(simdata->c, m, n, p) *
                           GETVALUE(simdata->c, m, n, p) * dtdx;

        double dvx = GETVALUE(simdata->vxold, m, n, p);
        double dvy = GETVALUE(simdata->vyold, m, n, p);
        double dvz = GETVALUE(simdata->vzold, m, n, p);
        
        dvx -= m > 0 ? GETVALUE(simdata->vxold, m - 1, n, p) : process->vx_bdy[1][p*numnodesy+n];
        dvy -= process->vy_bdy[1][p*numnodesx + m];
        dvz -= p > 0 ? GETVALUE(simdata->vzold, m, n, p - 1) : process->vz_bdy[1][n*numnodesx + m];

        double value = GETVALUE(simdata->pold, m, n, p) - rhoc2dtdx * (dvx + dvy + dvz);

        process->py_bdy[0][p*numnodesx+m] = value;
        SETVALUE(simdata->pnew, m, n, p, value);
    }
  }
  MPI_Isend(process->py_bdy[0], numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[DOWN], 1, process->world->cart_comm, &requesty_p);
  
  int p = 0;
  for (int n = 0; n < numnodesy; n++) {
    for (int m = 0; m < numnodesx; m++) {
        double rhoc2dtdx = GETVALUE(simdata->rho, m, n, p) *
                           GETVALUE(simdata->c, m, n, p) *
                           GETVALUE(simdata->c, m, n, p) * dtdx;

        double dvx = GETVALUE(simdata->vxold, m, n, p);
        double dvy = GETVALUE(simdata->vyold, m, n, p);
        double dvz = GETVALUE(simdata->vzold, m, n, p);

        dvx -= m > 0 ? GETVALUE(simdata->vxold, m - 1, n, p) : process->vx_bdy[1][p*numnodesy+n];
        dvy -= n > 0 ? GETVALUE(simdata->vyold, m, n - 1, p) : process->vy_bdy[1][p*numnodesx + m];
        dvz -= process->vz_bdy[1][n*numnodesx + m];

        double value = GETVALUE(simdata->pold, m, n, p) - rhoc2dtdx * (dvx + dvy + dvz);

        process->pz_bdy[0][n*numnodesx+m] = value;
        SETVALUE(simdata->pnew, m, n, p, value);
    }
  }
  MPI_Isend(process->pz_bdy[0], numnodesy*numnodesx, MPI_DOUBLE, process->neighbors[BACKWARD], 2, process->world->cart_comm, &requestz_p);

  for (int p = 1; p < numnodesz; p++) {
    for (int n = 1; n < numnodesy; n++) {
      for (int m = 1; m < numnodesx; m++) {
        double c = GETVALUE(simdata->c, m, n, p);
        double rho = GETVALUE(simdata->rho, m, n, p);

        double rhoc2dtdx = rho * c * c * dtdx;

        double dvx = GETVALUE(simdata->vxold, m, n, p);
        double dvy = GETVALUE(simdata->vyold, m, n, p);
        double dvz = GETVALUE(simdata->vzold, m, n, p);

        dvx -= GETVALUE(simdata->vxold, m - 1, n, p);
        dvy -= GETVALUE(simdata->vyold, m, n - 1, p);
        dvz -= GETVALUE(simdata->vzold, m, n, p - 1);

        double prev_p = GETVALUE(simdata->pold, m, n, p);

        SETVALUE(simdata->pnew, m, n, p, prev_p - rhoc2dtdx * (dvx + dvy + dvz));
      }
    }
  }
}

void update_velocities(simulation_data_t *simdata, process_s *process) {
  const double dtdx = simdata->params.dt / simdata->params.dx;

  const int numnodesx = NUMNODESX(simdata->vxold);
  const int numnodesy = NUMNODESY(simdata->vxold);
  const int numnodesz = NUMNODESZ(simdata->vxold);

  /*
  FROM LECTURE #1 about the memory hierarchy we know that it's better to avoid fetch data that will not be used if we have to used it later on.
  So we preferer to use the fact that when fetching data we get the full cache line. Since the data structure is organize such that we store 'm' near each other (and then n an then p)
  it is better to do the 3 loop in the order p - n - m (we use the next m since we already got it from the cache line of m).

  We know this from the form of the 'INDEX 3D' macro : 
  #define INDEX3D(grid, m, n, p)                                                 
    ((size_t)grid.numnodesy * grid.numnodesx * (p) + grid.numnodesx * (n) + (m))

  We gain a factor of about 5 time faster !!
  */

  
  MPI_Recv(process->px_bdy[1], numnodesy*numnodesz, MPI_DOUBLE, process->neighbors[RIGHT], 0, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->py_bdy[1], numnodesx*numnodesz, MPI_DOUBLE, process->neighbors[UP], 1, process->world->cart_comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->pz_bdy[1], numnodesy*numnodesx, MPI_DOUBLE, process->neighbors[FORWARD], 2, process->world->cart_comm, MPI_STATUS_IGNORE);

  int p = numnodesz - 1;
  for (int n = 0; n < numnodesy - 1; n++) {
    for (int m = 0; m < numnodesx - 1; m++) {
        double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);

        double p_mnq = GETVALUE(simdata->pnew, m, n, p);
        double dpx, dpy, dpz;

        dpx = GETVALUE(simdata->pnew, m+1, n, p) - p_mnq;
        dpy = GETVALUE(simdata->pnew, m, n+1, p) - p_mnq;
        dpz = process->pz_bdy[1][n*numnodesx + m] - p_mnq;
        

        double prev_vx = GETVALUE(simdata->vxold, m, n, p);
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        double value_x = prev_vx - dtdxrho * dpx;
        double value_y = prev_vy - dtdxrho * dpy;
        double value_z = prev_vz - dtdxrho * dpz;
        
        SETVALUE(simdata->vxnew, m, n, p, value_x);
        SETVALUE(simdata->vynew, m, n, p, value_y);
        SETVALUE(simdata->vznew, m, n, p, value_z);

        process->vz_bdy[0][n*numnodesx+m] = value_z;
    }
  }
  int n = numnodesy - 1;
  for (int p = 0; p < numnodesz - 1; p++) {
    for (int m = 0; m < numnodesx - 1; m++) {
        double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);

        double p_mnq = GETVALUE(simdata->pnew, m, n, p);
        double dpx, dpy, dpz;

        dpx = GETVALUE(simdata->pnew, m+1, n, p) - p_mnq;
        dpy = process->py_bdy[1][p*numnodesx + m] - p_mnq;
        dpz = GETVALUE(simdata->pnew, m, n, p+1) - p_mnq;
        

        double prev_vx = GETVALUE(simdata->vxold, m, n, p);
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        double value_x = prev_vx - dtdxrho * dpx;
        double value_y = prev_vy - dtdxrho * dpy;
        double value_z = prev_vz - dtdxrho * dpz;
        
        SETVALUE(simdata->vxnew, m, n, p, value_x);
        SETVALUE(simdata->vynew, m, n, p, value_y);
        SETVALUE(simdata->vznew, m, n, p, value_z);

        process->vy_bdy[0][p*numnodesx+m] = value_y;
    }
  }
  int m = numnodesx - 1;
  for (int p = 0; p < numnodesz - 1; p++) {
    for (int n = 0; n < numnodesy - 1; n++) {
        double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);

        double p_mnq = GETVALUE(simdata->pnew, m, n, p);
        double dpx, dpy, dpz;

        dpx = process->px_bdy[1][p*numnodesy + n] - p_mnq;
        dpy = GETVALUE(simdata->pnew, m, n+1, p) - p_mnq;
        dpz = GETVALUE(simdata->pnew, m, n, p+1) - p_mnq;
        

        double prev_vx = GETVALUE(simdata->vxold, m, n, p);
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        double value_x = prev_vx - dtdxrho * dpx;
        double value_y = prev_vy - dtdxrho * dpy;
        double value_z = prev_vz - dtdxrho * dpz;
        
        SETVALUE(simdata->vxnew, m, n, p, value_x);
        SETVALUE(simdata->vynew, m, n, p, value_y);
        SETVALUE(simdata->vznew, m, n, p, value_z);

        process->vx_bdy[0][p*numnodesy+n] = value_x;
    }
  }
  p = numnodesz - 1;
  n = numnodesy - 1;
  for (int m = 0; m < numnodesx; m++) {
        double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);

        double p_mnq = GETVALUE(simdata->pnew, m, n, p);
        double dpx, dpy, dpz;

        if(m == numnodesx - 1)
            dpx = process->px_bdy[1][p*numnodesy + n] - p_mnq;
        else
            dpx = GETVALUE(simdata->pnew, m+1, n, p) - p_mnq;
        dpy = process->py_bdy[1][p*numnodesx + m] - p_mnq;
        dpz = process->pz_bdy[1][n*numnodesx + m] - p_mnq;
        

        double prev_vx = GETVALUE(simdata->vxold, m, n, p);
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        double value_x = prev_vx - dtdxrho * dpx;
        double value_y = prev_vy - dtdxrho * dpy;
        double value_z = prev_vz - dtdxrho * dpz;
        
        SETVALUE(simdata->vxnew, m, n, p, value_x);
        SETVALUE(simdata->vynew, m, n, p, value_y);
        SETVALUE(simdata->vznew, m, n, p, value_z);

        if(m == numnodesx - 1)
            process->vx_bdy[0][p*numnodesy+n] = value_x;
        process->vy_bdy[0][p*numnodesx+m] = value_y;
        process->vz_bdy[0][n*numnodesx+m] = value_z;
  }
  for (int p = 0; p < numnodesz - 1; p++) {
    for (int n = 0; n < numnodesy - 1; n++) {
      for (int m = 0; m < numnodesx - 1; m++) {
        double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);

        double p_mnq = GETVALUE(simdata->pnew, m, n, p);
        double dpx, dpy, dpz;

        dpx = GETVALUE(simdata->pnew, m+1, n, p) - p_mnq;
        dpy = GETVALUE(simdata->pnew, m, n+1, p) - p_mnq;
        dpz = GETVALUE(simdata->pnew, m, n, p+1) - p_mnq;
        

        double prev_vx = GETVALUE(simdata->vxold, m, n, p);
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        double value_x = prev_vx - dtdxrho * dpx;
        double value_y = prev_vy - dtdxrho * dpy;
        double value_z = prev_vz - dtdxrho * dpz;
        
        SETVALUE(simdata->vxnew, m, n, p, value_x);
        SETVALUE(simdata->vynew, m, n, p, value_y);
        SETVALUE(simdata->vznew, m, n, p, value_z);
      }
    }
  }
}