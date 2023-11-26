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
