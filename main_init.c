#include "gmsh.h"
#include "dg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mkl.h>

int main(int argc, char **argv){
  
  // ======================================================================
  // 1] Read Setup File
  // ======================================================================
  
  // Parameters in setupfile
  int iOutputGmsh;        // Number of time steps between gmsh output
  int N;                  // Degree of polynomial bases
  double FinalTime;       // Final time of the simulation
  double *_ParamDB;       // Database for parameters
  int PARAM;              // Number of lines in database
  
  // Read the parameters
  char fileName[] = "setup";
  FILE* setupfile = NULL;
  setupfile = fopen(fileName, "r");
  if (setupfile == NULL){
    printf("ERROR: Setup file 'setup' not opened\n");
    exit(1);
  }
  
  char line[128] = "";
  while(fgets(line, 128, setupfile) != NULL){
    if(strstr(line, "$OutputStep") != NULL)
      fscanf(setupfile, "%i\n", &iOutputGmsh);
    if(strstr(line, "$PolynomialDegree") != NULL)
      fscanf(setupfile, "%i\n", &N);
    if(strstr(line, "$Duration") != NULL)
      fscanf(setupfile, "%lf\n", &FinalTime);
    if(strstr(line, "$Param") != NULL){
      fscanf(setupfile, "%i\n", &PARAM);
      _ParamDB = malloc(3*PARAM*sizeof(double));
      for(int i=0; i<PARAM; i++)
        fscanf(setupfile, "%lf %lf %lf", &_ParamDB[3*i+0], &_ParamDB[3*i+1], &_ParamDB[3*i+2]);
    }
  }
  
  
  // ======================================================================
  // 2] Load Mesh + Compute Connectivity Element-to-Element
  // ======================================================================
  
  int K = 0;           // Number of elements in total mesh
  double *_VX = NULL;  // Coord 'x' of all vertices (original gmsh ordering) [allVertices]
  double *_VY = NULL;  // Coord 'y' of all vertices (original gmsh ordering) [allVertices]
  double *_VZ = NULL;  // Coord 'z' of all vertices (original gmsh ordering) [allVertices]
  int *_ETag = NULL;   // List of element tag (gmsh: physical tag) [K]
  int *_EMsh = NULL;   // List of orginal gmsh numbers of elements [K]
  int *_EToV = NULL;   // Element-to-vertex connectivity array [K,NVertTet]
  int *_EToE = NULL;   // Element-to-element connectivity array [K,NFacesTet]
  int *_EToF = NULL;   // Element-to-face connectivity array [K,NFacesTet]
  
  // Load mesh (output: K, _VX, _VY, _VZ, _ETag, _EMsh, _EToV)
  char meshName[] = "mesh3.msh";
  loadMeshGmsh(meshName, &K, &_VX, &_VY, &_VZ, &_ETag, &_EMsh, &_EToV);
  
  // Build connectivity arrays (output: _EToE, _EToF)
  buildEToE(K, _EToV, &_EToE, &_EToF);
  
  
  // ======================================================================
  // 3] Build Discrete Representation + Compute Connectivity Node-to-Node
  // ======================================================================
  
  // Parameters
  int Np  = (N+1)*(N+2)*(N+3)/6;  // Number of nodes per element
  int Nfp = (N+1)*(N+2)/2;        // Number of nodes per face
  
#include "data.h"
  
  // Local coordinates for reference nodes
  double *_r;  // Coord 'r' for each reference node [Np]
  double *_s;  // Coord 's' for each reference node [Np]
  double *_t;  // Coord 't' for each reference node [Np]
  switch(N){
    case 1: _r = ref_r_01; _s = ref_s_01; _t = ref_t_01; break;
    case 2: _r = ref_r_02; _s = ref_s_02; _t = ref_t_02; break;
    case 3: _r = ref_r_03; _s = ref_s_03; _t = ref_t_03; break;
    case 4: _r = ref_r_04; _s = ref_s_04; _t = ref_t_04; break;
    case 5: _r = ref_r_05; _s = ref_s_05; _t = ref_t_05; break;
    case 6: _r = ref_r_06; _s = ref_s_06; _t = ref_t_06; break;
    case 7: _r = ref_r_07; _s = ref_s_07; _t = ref_t_07; break;
    case 8: _r = ref_r_08; _s = ref_s_08; _t = ref_t_08; break;
  }
  
  // Global coordinates for physical nodes
  double *_x = malloc(K*Np*sizeof(double));  // Coord 'x' for each physical node [K,Np]
  double *_y = malloc(K*Np*sizeof(double));  // Coord 'y' for each physical node [K,Np]
  double *_z = malloc(K*Np*sizeof(double));  // Coord 'z' for each physical node [K,Np]
  for(int k=0, cntv=0, cntf=0; k<K; ++k){
    int vGlo1 = _EToV[k*NVertTet + 0];
    int vGlo2 = _EToV[k*NVertTet + 1];
    int vGlo3 = _EToV[k*NVertTet + 2];
    int vGlo4 = _EToV[k*NVertTet + 3];
    double x1=_VX[vGlo1], y1=_VY[vGlo1], z1=_VZ[vGlo1];
    double x2=_VX[vGlo2], y2=_VY[vGlo2], z2=_VZ[vGlo2];
    double x3=_VX[vGlo3], y3=_VY[vGlo3], z3=_VZ[vGlo3];
    double x4=_VX[vGlo4], y4=_VY[vGlo4], z4=_VZ[vGlo4];
    for(int n=0; n<Np; ++n){
      double lambda1 = - (1+_r[n]+_s[n]+_t[n]);
      double lambda2 = 1+_r[n];
      double lambda3 = 1+_s[n];
      double lambda4 = 1+_t[n];
      _x[k*Np + n] = 0.5*(lambda1*x1 + lambda2*x2 + lambda3*x3 + lambda4*x4 );
      _y[k*Np + n] = 0.5*(lambda1*y1 + lambda2*y2 + lambda3*y3 + lambda4*y4 );
      _z[k*Np + n] = 0.5*(lambda1*z1 + lambda2*z2 + lambda3*z3 + lambda4*z4 );
      ++cntv;
    }
  }
  
  // Tolerances
  double NODETOL = 1e-4;
  double TOL = 1.e-6;
  
  // Build faceNode mask
  int *_Fmask = malloc(NFacesTet*Nfp*sizeof(int));  // Indices of faces nodes in the list of element nodes
  for(int n=0, cnt1=0, cnt2=0, cnt3=0, cnt4=0; n<Np; ++n){
    if(fabs(1+_t[n])<NODETOL)
      _Fmask[0*Nfp + cnt1++] = n;
    if(fabs(1+_s[n])<NODETOL)
      _Fmask[1*Nfp + cnt2++] = n;
    if(fabs(1+_r[n]+_s[n]+_t[n])<NODETOL)
      _Fmask[2*Nfp + cnt3++] = n;
    if(fabs(1+_r[n])<NODETOL)
      _Fmask[3*Nfp + cnt4++] = n;
  }
  
  // Build faceNode-to-faceNode connectivity (mapP)
  int *_mapP = NULL;
  buildNfToNf(K, Nfp, Np, _Fmask, _EToE, _EToF, _x, _y, _z, TOL, &_mapP);
  
  
  // ======================================================================
  // 4] Build Geometric Factors & Build Elemental Matrices
  // ======================================================================
  
  // Geometric factors
  double *_rstxyz = NULL;   // ... for volume terms  [K,9]
  double *_Fscale = NULL;   // ... for surface terms [K,NFacesTet]
  double *_nx = NULL;       // x-component of normal to the face [K,NFacesTet]
  double *_ny = NULL;       // y-component of normal to the face [K,NFacesTet]
  double *_nz = NULL;       // z-component of normal to the face [K,NFacesTet]
  
  buildGeomFactors(K, _EToV, _VX,_VY,_VZ, &_rstxyz, &_Fscale, &_nx, &_ny, &_nz);
  
  // Elemental matrices
  double *_Dr = NULL;    // r-derivative matrix
  double *_Ds = NULL;    // s-derivative matrix
  double *_Dt = NULL;    // t-derivative matrix
  double *_LIFT = NULL;  // lift matrix
  
  buildElemMatrices(N, Np, Nfp, _r, _s, _t, _Fmask, &_Dr, &_Ds, &_Dt, &_LIFT);
  
  // ======================================================================
  // 5] Load Physical Parameters & Compute Time Step
  // ======================================================================
  
  // Physical parameters (c, rho)
  double *_paramMap = malloc(2*K*sizeof(double));
  for(int k=0; k<K; k++){
    int i=0;
    while(i<PARAM){
      if((int)_ParamDB[3*i+0] == _ETag[k]){
        _paramMap[2*k+0] = _ParamDB[3*i+1];
        _paramMap[2*k+1] = _ParamDB[3*i+2];
        break;
      }
      i++;
    }
    if(i==PARAM){
      printf("ERROR: ETag %i not found in parameter database (in setup file)\n", _ETag[k]);
      exit(1);
    }
  }
  free(_ParamDB);
  exportParamGmsh(_paramMap, K, 2, _EMsh);
  
  // Time step
  double dt = 1e9;
  for(int k=0; k<K; k++){
    for(int f=0; f<NFacesTet; ++f){
      double val = 1./(_paramMap[2*k+0] * (N+1)*(N+1) * _Fscale[k*NFacesTet + f]);
      if(val < dt) dt = val;
    }
  }
  double CFL = 0.75;                  // CFL
  dt *= CFL;                          // Time step
  int Nsteps = ceil(FinalTime/dt);    // Number of global time steps
  
  printf("INFO: DT = %f and Nsteps = %i\n", dt, Nsteps);
  printf("INFO: Nfields.Np.K = %i\n", 4*Np*K);
  
  
  // ======================================================================
  // 6] Memory storage for fields, RHS and residual at nodes
  // ======================================================================
  
  // Parameter
  int Nfields = 4;
  
  // Allocation
  double *_valQ = malloc(K*Np*Nfields*sizeof(double));  // Storage for values of fields at each node of each element
  double *_rhsQ = malloc(K*Np*Nfields*sizeof(double));  // Storage for RHS at each node of each element
  double *_resQ = malloc(K*Np*Nfields*sizeof(double));  // Storage for residual at each node of each element
  
  // Initialization
  for(int k=0; k<K; ++k){
    for(int n=0; n<Np; ++n){
      for(int iField=0; iField<Nfields; ++iField){
        _valQ[k*Np*Nfields + n*Nfields + iField] = 0.;
        _rhsQ[k*Np*Nfields + n*Nfields + iField] = 0.;
        _resQ[k*Np*Nfields + n*Nfields + iField] = 0.;
      }
      double x = _x[k*Np + n];
      double y = _y[k*Np + n];
      double z = _z[k*Np + n];
      _valQ[k*Np*Nfields + n*Nfields + 0] = exp(-(x*x+y*y+z*z)/0.1f);
    }
  }
  
  // Export initial solution
  if(iOutputGmsh > 0)
    exportSolGmsh(_valQ, K, N, Np, Nfields, _r,_s,_t, _EMsh, 0, 0.);
  
  
  // ======================================================================
  // 7] RUN
  // ======================================================================
  
  // Global time iteration
  unsigned MKL_INT64 timeBegin;
  mkl_get_cpu_clocks(&timeBegin);
  for(int nGlo=0; nGlo<Nsteps; ++nGlo){
    double runTime = nGlo*dt;  // Time at the beginning of the step
    
    // Local time iteration (5 sub-iteration for low-storage RK44)
    for(int nLoc=0; nLoc<5; ++nLoc){
      double a = rk4a[nLoc];
      double b = rk4b[nLoc];
      
      // ======================== (1) UPDATE RHS
      for(int k=0; k<K; ++k){
        
        double c   = _paramMap[2*k+0];
        double rho = _paramMap[2*k+1];
        
        // ======================== (1.1) UPDATE PENALTY VECTOR
        
        double *s_p_flux = malloc(Nfp*NFacesTet*sizeof(double));
        double *s_u_flux = malloc(Nfp*NFacesTet*sizeof(double));
        double *s_v_flux = malloc(Nfp*NFacesTet*sizeof(double));
        double *s_w_flux = malloc(Nfp*NFacesTet*sizeof(double));
        
        for(int f=0; f<NFacesTet; f++){
          
          // Fetch normal
          double nx = _nx[k*NFacesTet + f];
          double ny = _ny[k*NFacesTet + f];
          double nz = _nz[k*NFacesTet + f];
          
          // Fetch medium parameters
          int k2 = _EToE[k*NFacesTet + f];
          double cP   = _paramMap[2*k2 + 0];
          double rhoP = _paramMap[2*k2 + 1];
          double param_flux=c/(rhoP*cP+rho*c); //Ce paramètre ne varie pas dans chaque boucle suivante, on factorise donc ce paramètre et on le remplace dans les appels des flux
          
          // Compute penalty terms
          for(int nf=f*Nfp; nf<(f+1)*Nfp; nf++){
            
            int n1 = _Fmask[nf];                // Index of node in current element
            int n2 = _mapP[nf + k*NFacesTet*Nfp];  // Index of node in neighbor element
            
            // Load values 'minus' corresponding to current element
            double pM = _valQ[k*Np*Nfields + n1*Nfields + 0];
            double uM = _valQ[k*Np*Nfields + n1*Nfields + 1];
            double vM = _valQ[k*Np*Nfields + n1*Nfields + 2];
            double wM = _valQ[k*Np*Nfields + n1*Nfields + 3];
            double nMdotuM = (nx*uM + ny*vM + nz*wM);
            
            if(n2 >= 0){ // ... if there is a neighbor element ...
              
              // Load values 'plus' corresponding to neighbor element
              double pP = _valQ[k2*Np*Nfields + n2*Nfields + 0];
              double uP = _valQ[k2*Np*Nfields + n2*Nfields + 1];
              double vP = _valQ[k2*Np*Nfields + n2*Nfields + 2];
              double wP = _valQ[k2*Np*Nfields + n2*Nfields + 3];
              double nMdotuP = (nx*uP + ny*vP + nz*wP);
              
              // Penalty terms for interface between two elements
              
              

              s_p_flux[nf] = c / (1./(rhoP*cP) + 1./(rho*c)) * ( (nMdotuP-nMdotuM) - 1./(rhoP*cP) * (pP-pM));
              s_u_flux[nf] = nx * c/(rhoP*cP+rho*c) * ( (pP-pM) - rhoP*cP * (nMdotuP-nMdotuM) );
              s_v_flux[nf] = ny * c/(rhoP*cP+rho*c) * ( (pP-pM) - rhoP*cP * (nMdotuP-nMdotuM) );
              s_w_flux[nf] = nz * c/(rhoP*cP+rho*c) * ( (pP-pM) - rhoP*cP * (nMdotuP-nMdotuM) );
            }
            else{
              // Penalty terms for boundary of the domain
              double tmp = -2./(rhoP*cP) * pM;            // Homogeneous Dirichlet on 'p'
              //double tmp = 2*nMdotuM;                   // Homogeneous Dirichlet on 'u'
              //double tmp = nMdotuM - 1./(rhoP*cP) * pM; // ABC
              
              s_p_flux[nf] = -c / (1./(rhoP*cP) + 1./(rho*c)) * tmp;
              s_u_flux[nf] = nx * c/(rhoP*cP+rho*c) * tmp;
              s_v_flux[nf] = ny * c/(rhoP*cP+rho*c) * tmp;
              s_w_flux[nf] = nz * c/(rhoP*cP+rho*c) * tmp;
            }
          }
        }
        
        // ======================== (1.2) COMPUTING VOLUME TERMS
        
        // Load geometric factors
        double rx = _rstxyz[k*9+0];
        double ry = _rstxyz[k*9+1];
        double rz = _rstxyz[k*9+2];
        double sx = _rstxyz[k*9+3];
        double sy = _rstxyz[k*9+4];
        double sz = _rstxyz[k*9+5];
        double tx = _rstxyz[k*9+6];
        double ty = _rstxyz[k*9+7];
        double tz = _rstxyz[k*9+8];
        
        // Load fields
        double *s_p = malloc(Np*sizeof(double));
        double *s_u = malloc(Np*sizeof(double));
        double *s_v = malloc(Np*sizeof(double));
        double *s_w = malloc(Np*sizeof(double));
        for(int n=0; n<Np; ++n){
          s_p[n] = _valQ[k*Np*Nfields + n*Nfields + 0];
          s_u[n] = _valQ[k*Np*Nfields + n*Nfields + 1];
          s_v[n] = _valQ[k*Np*Nfields + n*Nfields + 2];
          s_w[n] = _valQ[k*Np*Nfields + n*Nfields + 3];
        }
        
        // Compute mat-vec product for surface term
        for(int n=0; n<Np; ++n){
          double dpdr = 0, dpds = 0, dpdt = 0;
          double dudr = 0, duds = 0, dudt = 0;
          double dvdr = 0, dvds = 0, dvdt = 0;
          double dwdr = 0, dwds = 0, dwdt = 0;
          for(int m = 0; m < Np; ++m){
            double Dr = _Dr[n + m*Np];
            dpdr += Dr * s_p[m];
            dudr += Dr * s_u[m];
            dvdr += Dr * s_v[m];
            dwdr += Dr * s_w[m];
            double Ds = _Ds[n + m*Np];
            dpds += Ds * s_p[m];
            duds += Ds * s_u[m];
            dvds += Ds * s_v[m];
            dwds += Ds * s_w[m];
            double Dt = _Dt[n + m*Np];
            dpdt += Dt * s_p[m];
            dudt += Dt * s_u[m];
            dvdt += Dt * s_v[m];
            dwdt += Dt * s_w[m];
          }
          
          double dpdx = rx*dpdr + sx*dpds + tx*dpdt;
          double dpdy = ry*dpdr + sy*dpds + ty*dpdt;
          double dpdz = rz*dpdr + sz*dpds + tz*dpdt;
          double dudx = rx*dudr + sx*duds + tx*dudt;
          double dvdy = ry*dvdr + sy*dvds + ty*dvdt;
          double dwdz = rz*dwdr + sz*dwds + tz*dwdt;
          double divU = dudx + dvdy + dwdz;
          
          // Compute RHS (only part corresponding to volume terms)
          _rhsQ[k*Np*Nfields + n*Nfields + 0] = -c*c*rho * divU;
          _rhsQ[k*Np*Nfields + n*Nfields + 1] = -1.f/rho * dpdx;
          _rhsQ[k*Np*Nfields + n*Nfields + 2] = -1.f/rho * dpdy;
          _rhsQ[k*Np*Nfields + n*Nfields + 3] = -1.f/rho * dpdz;
        } 
        
        free(s_p);
        free(s_u);
        free(s_v);
        free(s_w);
        
        // ======================== (1.3) COMPUTING SURFACE TERMS
        
        for(int n=0; n<Np; ++n){
          
          for(int f=0; f<NFacesTet; f++){
            
            // Compute mat-vec product for surface term
            double p_lift = 0.;
            double u_lift = 0.;
            double v_lift = 0.;
            double w_lift = 0.;
            for(int m=f*Nfp; m<(f+1)*Nfp; m++){
              double tmp = _LIFT[n*NFacesTet*Nfp + m];
              p_lift += tmp*s_p_flux[m];
              u_lift += tmp*s_u_flux[m];
              v_lift += tmp*s_v_flux[m];
              w_lift += tmp*s_w_flux[m];
            }
            
            // Load geometric factor
            double Fscale = _Fscale[k*NFacesTet + f];
            
            // Update RHS (with part corresponding to surface terms)
            _rhsQ[k*Np*Nfields + n*Nfields + 0] -= p_lift * Fscale;
            _rhsQ[k*Np*Nfields + n*Nfields + 1] -= u_lift * Fscale;
            _rhsQ[k*Np*Nfields + n*Nfields + 2] -= v_lift * Fscale;
            _rhsQ[k*Np*Nfields + n*Nfields + 3] -= w_lift * Fscale;
          }
        }
        
        free(s_p_flux);
        free(s_u_flux);
        free(s_v_flux);
        free(s_w_flux);
      }
      
      // ======================== (2) UPDATE RESIDUAL + FIELDS
      
      for(int k=0; k<K; ++k){
        for(int n=0; n<Np; ++n){
          for(int iField=0; iField<Nfields; ++iField){
            int id = k*Np*Nfields + n*Nfields + iField;
            _resQ[id] = a*_resQ[id] + dt*_rhsQ[id];
            _valQ[id] = _valQ[id] + b*_resQ[id];
          }
        }
      }
    }
    
    // Export solution
    if((iOutputGmsh > 0) && ((nGlo+1) % iOutputGmsh == 0))
      exportSolGmsh(_valQ, K, N, Np, Nfields, _r,_s,_t, _EMsh, nGlo+1, runTime+dt);
  }
  
  // Export final solution
  if(iOutputGmsh > 0)
    exportSolGmsh(_valQ, K, N, Np, Nfields, _r,_s,_t, _EMsh, Nsteps, Nsteps*dt);
  
  unsigned MKL_INT64 timeEnd;
  mkl_get_cpu_clocks(&timeEnd);
  double timeTotal4 = (double)(timeEnd-timeBegin)/mkl_get_clocks_frequency()/1e9;
  printf("%f\n",timeTotal4);
  // ======================================================================
  // 8] Deallocation
  // ======================================================================
  
  free(_VX);
  free(_VY);
  free(_VZ);
  free(_EToV);
  free(_ETag);
  free(_EMsh);
  free(_EToE);
  free(_EToF);
  free(_x);
  free(_y);
  free(_z);
  free(_Fmask);
  free(_mapP);
  free(_rstxyz);
  free(_Fscale);
  free(_nx);
  free(_ny);
  free(_nz);
  free(_Dr);
  free(_Ds);
  free(_Dt);
  free(_LIFT);
  free(_paramMap);
  free(_valQ);
  free(_rhsQ);
  free(_resQ);
  
  return 0;
}
