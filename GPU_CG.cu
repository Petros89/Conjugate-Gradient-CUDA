#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cublas.h>
#include <cublas_v2.h>
#include <sys/time.h>
#include <sys/resource.h>


#define NDim 3
#define MNE 8

#define BLOCKSIZE 256

#define I3D(ni,nj,i,j,k) (((nj*ni)*(k))+((ni)*(j)) + i)

const char *FEM3Dname = "Finite Element Formulation - 3D Structured Hexahedral";
const char *sSDKname  = "Conjugate Gradient Multi-Device CG Linear Solver";

// get time structure
double get_time()
{
  struct timeval timeval_time;
  gettimeofday(&timeval_time,NULL);
  return (double)timeval_time.tv_sec + (double)timeval_time.tv_usec*1e-6;
}

// get RAM structure
struct rusage r_usage;

void KM_matrix(double delta_x, double delta_y,double delta_z, double *Ke){
//  double N[8];
  double DN[3][8];
  double gxyz[3][3];
  double wg[3][3],wxyz;
  int i,j,k,i1,j1;
  double a,b,c;
  double x,y,z,sommeK;

  gxyz[0][2]=-0.774596669241483;
  gxyz[1][2]=0;
  gxyz[2][2]=0.774596669241483;

  wg[0][0]=2;
  wg[0][1]=1;
  wg[0][2]=0.555555555555556;
  wg[1][0]=0;
  wg[1][1]=1;
  wg[1][2]=0.888888888888889;
  wg[2][0]=0;
  wg[2][1]=0;
  wg[2][2]=0.555555555555556;

a=0.5*delta_x;
b=0.5*delta_y;
c=0.5*delta_z;
        printf("boucle\n");


  for (i1=0;i1<8;i1++){
      for (j1=0;j1<8;j1++){
sommeK=0;
          for (i=0;i<3;i++){
            for (j=0;j<3;j++){
              for (k=0;k<3;k++){

                x=gxyz[i][2];
                y=gxyz[j][2];
                z=gxyz[k][2];
                wxyz=wg[i][2]*wg[j][2]*wg[k][2];

                DN[0][0]=-(1-y)*(1-z)*(1/a)*0.125;   DN[1][0]=-(1-x)*(1-z)*(1/b)*0.125;    DN[2][0]=-(1-x)*(1-y)*(1/c)*0.125;
                DN[0][1]=(1-y)*(1-z)*(1/a)*0.125;    DN[1][1]=-(1+x)*(1-z)*(1/b)*0.125;    DN[2][1]=-(1+x)*(1-y)*(1/c)*0.125;
                DN[0][2]=(1+y)*(1-z)*(1/a)*0.125;    DN[1][2]=(1+x)*(1-z)*(1/b)*0.125;     DN[2][2]=-(1+x)*(1+y)*(1/c)*0.125;
                DN[0][3]=-(1+y)*(1-z)*(1/a)*0.125;   DN[1][3]=(1-x)*(1-z)*(1/b)*0.125;     DN[2][3]=-(1-x)*(1+y)*(1/c)*0.125;
                DN[0][4]=-(1-y)*(1+z)*(1/a)*0.125;   DN[1][4]=-(1-x)*(1+z)*(1/b)*0.125;    DN[2][4]=(1-x)*(1-y)*(1/c)*0.125;
                DN[0][5]=(1-y)*(1+z)*(1/a)*0.125;    DN[1][5]=-(1+x)*(1+z)*(1/b)*0.125;    DN[2][5]=(1+x)*(1-y)*(1/c)*0.125;
                DN[0][6]=(1+y)*(1+z)*(1/a)*0.125;    DN[1][6]=(1+x)*(1+z)*(1/b)*0.125;     DN[2][6]=(1+x)*(1+y)*(1/c)*0.125;
                DN[0][7]=-(1+y)*(1+z)*(1/a)*0.125;   DN[1][7]=(1-x)*(1+z)*(1/b)*0.125;     DN[2][7]=(1-x)*(1+y)*(1/c)*0.125;

                sommeK=sommeK+(DN[0][i1]*DN[0][j1]+DN[1][i1]*DN[1][j1]+DN[2][i1]*DN[2][j1])*wxyz*a*b*c;

              }
            }
          }
          Ke[i1+8*j1]=sommeK;
          printf("%f ",Ke[i1+8*j1]);
    }
      printf("\n");
  }


}


void MatVec(double *Ke, double *U, double *r, int *face, size_t numElems, int numNodes, int *edofMat)
{

  // Set residual to zero
  for (int n=0; n<numNodes; n++){
    r[n] = 0.0;
  }

  // Loop over all elements
  for(size_t e=0; e<numElems; e++)
  {

    // For each node in element
    for(int id=0; id<MNE; id++)
    {
      size_t node = edofMat[id+MNE*e];

      double Au = 0.0;
      for(int j=0; j<MNE; j++)
      {
        size_t node2 = edofMat[j+MNE*e];
        Au = Au + Ke[id*MNE+j]*U[node2];
      }

      r[node] = r[node] + Au;

    }
  }

  //handle contraints [boundary faces]
  for(size_t e=0; e<numElems; e++)
  {

    if(face[e]==1){//j==0

      U[edofMat[0+MNE*e]]=0;
      U[edofMat[1+MNE*e]]=0;
      U[edofMat[4+MNE*e]]=0;
      U[edofMat[5+MNE*e]]=0;

      r[edofMat[0+MNE*e]]=0;
      r[edofMat[1+MNE*e]]=0;
      r[edofMat[4+MNE*e]]=0;
      r[edofMat[5+MNE*e]]=0;
    }

    if(face[e]==2){//i==0
      U[edofMat[0+MNE*e]]=0;
      U[edofMat[4+MNE*e]]=0;
      U[edofMat[3+MNE*e]]=0;
      U[edofMat[7+MNE*e]]=0;

      r[edofMat[0+MNE*e]]=0;
      r[edofMat[4+MNE*e]]=0;
      r[edofMat[3+MNE*e]]=0;
      r[edofMat[7+MNE*e]]=0;
    }

    if(face[e]==3){//j==(nj-1)
      U[edofMat[2+MNE*e]]=0;
      U[edofMat[3+MNE*e]]=0;
      U[edofMat[6+MNE*e]]=0;
      U[edofMat[7+MNE*e]]=0;

      r[edofMat[2+MNE*e]]=0;
      r[edofMat[3+MNE*e]]=0;
      r[edofMat[6+MNE*e]]=0;
      r[edofMat[7+MNE*e]]=0;
    }

    if(face[e]==4){//i==(ni-1)
      U[edofMat[1+MNE*e]]=0;
      U[edofMat[2+MNE*e]]=0;
      U[edofMat[5+MNE*e]]=0;
      U[edofMat[6+MNE*e]]=0;

      r[edofMat[1+MNE*e]]=0;
      r[edofMat[2+MNE*e]]=0;
      r[edofMat[5+MNE*e]]=0;
      r[edofMat[6+MNE*e]]=0;
    }

    if(face[e]==5){//k==0
      U[edofMat[0+MNE*e]]=0;
      U[edofMat[1+MNE*e]]=0;
      U[edofMat[2+MNE*e]]=0;
      U[edofMat[3+MNE*e]]=0;

      r[edofMat[0+MNE*e]]=0;
      r[edofMat[1+MNE*e]]=0;
      r[edofMat[2+MNE*e]]=0;
      r[edofMat[3+MNE*e]]=0;
    }

    if(face[e]==6){//k==nk-1
      U[edofMat[4+MNE*e]]=0;
      U[edofMat[5+MNE*e]]=0;
      U[edofMat[6+MNE*e]]=0;
      U[edofMat[7+MNE*e]]=0;

      r[edofMat[4+MNE*e]]=0;
      r[edofMat[5+MNE*e]]=0;
      r[edofMat[6+MNE*e]]=0;
      r[edofMat[7+MNE*e]]=0;
    }


  }
}

__global__ void ResetRes(double *r, int numNodes)
{
  size_t n = blockIdx.x * blockDim.x + threadIdx.x;

  if (n<numNodes){
     r[n] = 0.0;
  }
}



__global__ void GPUMatVec(int ni, int nj,  double *Ke, double *U, double *r, size_t size)

{

  int i,j,k,id1,edofMat2[MNE];



  __shared__ double Au[MNE],sum;


  size_t e = blockIdx.x * BLOCKSIZE + threadIdx.x;


  if (e<size)
  {

    k=e/(ni*nj);
    j=(e%(nj*ni))/ni;
    i=(e%(nj*ni))%ni;

    id1=(ni+1)*(nj+1)*k+(ni+1)*j+i;

    //local DOFs
    edofMat2[0]=id1;
    edofMat2[1]=id1+1;
    edofMat2[2]=id1+1+(ni+1);
    edofMat2[3]=id1+(ni+1);
    edofMat2[4]=id1 + (ni+1)*(nj+1);
    edofMat2[5]=id1+1 + (ni+1)*(nj+1);
    edofMat2[6]=id1+1+(ni+1) + (ni+1)*(nj+1);
    edofMat2[7]=id1+(ni+1) + (ni+1)*(nj+1);

	for (int i=0; i<MNE; i++){
           for (int j=0; j<MNE; j++){
              Au[j]= Ke[i*MNE+j]*U[edofMat2[j]];
           }

           sum=Au[0]+Au[1]+Au[2]+Au[3]+Au[4]+Au[5]+Au[6]+Au[7];

           atomicAdd(&r[edofMat2[i]], sum);
       }
  }
}


__global__ void DirichletBC(double *U, double *r, size_t numElems, int *edofMat, int *face)

{

  size_t e = blockIdx.x * blockDim.x + threadIdx.x;

  if (e<numElems)
  {

    if(face[e]==1){//j==0
      U[edofMat[0+MNE*e]]=0.0;
      U[edofMat[1+MNE*e]]=0.0;
      U[edofMat[4+MNE*e]]=0.0;
      U[edofMat[5+MNE*e]]=0.0;

      r[edofMat[0+MNE*e]]=0.0;
      r[edofMat[1+MNE*e]]=0.0;
      r[edofMat[4+MNE*e]]=0.0;
      r[edofMat[5+MNE*e]]=0.0;
    }

    if(face[e]==2){//i==0
      U[edofMat[0+MNE*e]]=0.0;
      U[edofMat[4+MNE*e]]=0.0;
      U[edofMat[3+MNE*e]]=0.0;
      U[edofMat[7+MNE*e]]=0.0;

      r[edofMat[0+MNE*e]]=0.0;
      r[edofMat[4+MNE*e]]=0.0;
      r[edofMat[3+MNE*e]]=0.0;
      r[edofMat[7+MNE*e]]=0.0;
    }

    if(face[e]==3){//j==(nj.0-1)
      U[edofMat[2+MNE*e]]=0.0;
      U[edofMat[3+MNE*e]]=0.0;
      U[edofMat[6+MNE*e]]=0.0;
      U[edofMat[7+MNE*e]]=0.0;

      r[edofMat[2+MNE*e]]=0.0;
      r[edofMat[3+MNE*e]]=0.0;
      r[edofMat[6+MNE*e]]=0.0;
      r[edofMat[7+MNE*e]]=0.0;
    }

    if(face[e]==4){//i==(ni.0-1)
      U[edofMat[1+MNE*e]]=0.0;
      U[edofMat[2+MNE*e]]=0.0;
      U[edofMat[5+MNE*e]]=0.0;
      U[edofMat[6+MNE*e]]=0.0;

      r[edofMat[1+MNE*e]]=0.0;
      r[edofMat[2+MNE*e]]=0.0;
      r[edofMat[5+MNE*e]]=0.0;
      r[edofMat[6+MNE*e]]=0.0;
    }

    if(face[e]==5){//k==0
      U[edofMat[0+MNE*e]]=0.0;
      U[edofMat[1+MNE*e]]=0.0;
      U[edofMat[2+MNE*e]]=0.0;
      U[edofMat[3+MNE*e]]=0.0;

      r[edofMat[0+MNE*e]]=0.0;
      r[edofMat[1+MNE*e]]=0.0;
      r[edofMat[2+MNE*e]]=0.0;
      r[edofMat[3+MNE*e]]=0.0;
    }

    if(face[e]==6){//k==nk-.01
      U[edofMat[4+MNE*e]]=0.0;
      U[edofMat[5+MNE*e]]=0.0;
      U[edofMat[6+MNE*e]]=0.0;
      U[edofMat[7+MNE*e]]=0.0;

      r[edofMat[4+MNE*e]]=0.0;
      r[edofMat[5+MNE*e]]=0.0;
      r[edofMat[6+MNE*e]]=0.0;
      r[edofMat[7+MNE*e]]=0.0;
    }
  }

}


__global__ void UpdateVec(double *p, double *r, double beta,  int numNodes)

{
  size_t n = blockIdx.x * blockDim.x + threadIdx.x;

  if (n<numNodes){

      p[n] = r[n] + beta * p[n];

    }
}




int main(int argc, char **argv)
{




  //////////////////////////////////////// End Read Mesh ////////////////////////////////////////
//GET GPU DEVICE
printf("Starting [%s]...\n", FEM3Dname);
int i,j,k,i0,idx,id1;

//size of the domain
const double domain_x=1;
const double domain_y=1;
const double domain_z=1;

//number of elements in each direction
const int ni=250;
const int nj=250;
const int nk=250;

printf("Mesh Element Size %e\n", float(ni*nj*nk));
//volume multiplier
const double Qvol=10;

//space steps
const double delta_x=(double)domain_x/ni;
const double delta_y=(double)domain_y/nj;
const double delta_z=(double)domain_z/nk;


//Total number of elements & nodes for the 3D structured uniform grid
const int numElems=ni*nj*nk;
const int numNodes=(ni+1)*(nj+1)*(nk+1);


printf("numNodes: %d\n", numNodes);
printf("numElems: %d\n", numElems);




//////////////////////////////////////// Solve Poisson ////////////////////////////////////////
  //define misc vars
  double rzold = 0.0;
  double alpha = 0.0;
  double beta  = 0.0;
  double malpha = -alpha;


  //Allocate CPU memory using MallocManaged
  double *Ke = (double *) malloc(MNE*MNE * sizeof(double));

  int *edofMat = (int *) malloc(MNE*numElems * sizeof(int));
  int *face = (int *) malloc(numElems * sizeof(int));

  double *U = (double *) malloc(numNodes * sizeof(double));
  double *r = (double *) malloc(numNodes * sizeof(double));
  double *b = (double *) malloc(numNodes * sizeof(double));
  double *p = (double *) malloc(numNodes * sizeof(double));
  double *Ap = (double *) malloc(numNodes * sizeof(double));


  KM_matrix(delta_x,delta_y,delta_z, Ke);


  // Assign right hand side
  for(size_t n=0; n<numNodes; n++)
  {
    b[n] = delta_x*delta_y*delta_z*Qvol;
  }


  /* Node/element connectivity */

  for(size_t e=0; e<numElems; e++)
  {
    //define loc2glob indexes
    k=e/(ni*nj);
    j=(e%(nj*ni))/ni;
    i=(e%(nj*ni))%ni;
    idx=ni*nj*k+ni*j+i;
    id1=(ni+1)*(nj+1)*k+(ni+1)*j+i;

    //local DOFs
    edofMat[0+MNE*e]=id1;
    edofMat[1+MNE*e]=id1+1;
    edofMat[2+MNE*e]=id1+1+(ni+1);
    edofMat[3+MNE*e]=id1+(ni+1);
    edofMat[4+MNE*e]=id1 + (ni+1)*(nj+1);
    edofMat[5+MNE*e]=id1+1 + (ni+1)*(nj+1);
    edofMat[6+MNE*e]=id1+1+(ni+1) + (ni+1)*(nj+1);
    edofMat[7+MNE*e]=id1+(ni+1) + (ni+1)*(nj+1);

  }

  printf("Last Node of Domain IS::  %d\n", edofMat[6+MNE*(numElems-1)]);
  /* Impose Boundary Conditions at 6 faces of the box */

  for(size_t e=0; e<numElems; e++)
  {
    k=e/(ni*nj);
    j=(e%(nj*ni))/ni;
    i=(e%(nj*ni))%ni;
    idx=ni*nj*k+ni*j+i;
    id1=(ni+1)*(nj+1)*k+(ni+1)*j+i;

    if(i==0 || j==0 || i==(ni-1) || j==(nj-1) || k==0 || k==(nk-1)){

    if(j==0) {
      face[idx]=1;
      b[edofMat[0+MNE*idx]]=0;  b[edofMat[1+MNE*idx]]=0;  b[edofMat[4+MNE*idx]]=0;   b[edofMat[5+MNE*idx]]=0;
    }

    if(i==0) {
      face[idx]=2;
      b[edofMat[0+MNE*idx]]=0; b[edofMat[4+MNE*idx]]=0; b[edofMat[3+MNE*idx]]=0;  b[edofMat[7+MNE*idx]]=0;
    }

    if(j==(nj-1)) {
      face[idx]=3;
      b[edofMat[2+MNE*idx]]=0; b[edofMat[3+MNE*idx]]=0; b[edofMat[6+MNE*idx]]=0;  b[edofMat[7+MNE*idx]]=0;
    }

    if(i==(ni-1)) {
      face[idx]=4;
        b[edofMat[1+MNE*idx]]=0; b[edofMat[2+MNE*idx]]=0; b[edofMat[5+MNE*idx]]=0;  b[edofMat[6+MNE*idx]]=0;
    }

    if(k==0) {
      face[idx]=5;
      b[edofMat[0+MNE*idx]]=0; b[edofMat[1+MNE*idx]]=0; b[edofMat[2+MNE*idx]]=0;  b[edofMat[3+MNE*idx]]=0;
    }

    if(k==(nk-1)) {
      face[idx]=6;
      b[edofMat[4+MNE*idx]]=0; b[edofMat[5+MNE*idx]]=0; b[edofMat[6+MNE*idx]]=0;  b[edofMat[7+MNE*idx]]=0;
    }

    }

  }


  //Initialize MatVec on the CPU
  MatVec(Ke, U, r, face, numElems, numNodes, edofMat);

  for(size_t n=0; n<numNodes; n++)
  {

      r[n] = b[n] - r[n];
      p[n] = r[n];
  }



  // control phi
  rzold = 0.0;

  for(size_t n=0; n<numNodes; n++)
  {
    rzold = rzold + r[n]*r[n];
  }
  printf("rzold=%lg\n",rzold);
  printf("***************** CPU PART IS OVER ***************\n");
  //End Initialization



/*------let's let the NVIDIA to take care allocation of matrices using the Unified Memory API----------*/


//GET GPU DEVICE
printf("\nStarting [%s]...\n", sSDKname);

  //Check out my device properties
	
  int nDevices;
  cudaGetDeviceCount(&nDevices);
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
  }


  int *face_d, *edofMat_d;
  double *U_d, *r_d, *p_d;
  double *Ap_d, *Ke_d;

  double *pAp_d;
  double *rznew_d;
  double *rNorm_d;

  // Allocate memory
  cudaMalloc((void **) &Ke_d, MNE*MNE * sizeof(double));
  cudaMalloc((void **) &edofMat_d, MNE*numElems * sizeof(int));
  cudaMalloc((void **) &face_d, numElems * sizeof(int));

  cudaMallocManaged((void **) &U_d, numNodes * sizeof(double));
  cudaMallocManaged((void **) &r_d, numNodes * sizeof(double));
  cudaMallocManaged((void **) &p_d, numNodes * sizeof(double));
  cudaMallocManaged((void **) &Ap_d, numNodes * sizeof(double));

  //additional variable for reduction operations
  cudaMallocManaged((void **) &pAp_d, numNodes * sizeof(double));
  cudaMallocManaged((void **) &rznew_d, numNodes * sizeof(double));
  cudaMallocManaged((void **) &rNorm_d, numNodes * sizeof(double));



  int  NBLOCKS   = (numElems+BLOCKSIZE-1) / BLOCKSIZE; //round up if n is not a multiple of blocksize
  printf("\nNumber of Blocks: %d\n",NBLOCKS);
  printf("numNodes: %d\n", numNodes);
  printf("numElems: %d\n", numElems);
  int  NBLOCKS2 = (numNodes+BLOCKSIZE-1) / BLOCKSIZE; //round up if n is not a multiple of blocksize
  const double TOL = 1.0e-6;




  //Send data to GPU [Ke_d,U_d,r_d,p_d,Ap_d,face_d,edofMat_d]
  cudaMemcpy(Ke_d, Ke, MNE*MNE * sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(U_d, U, numNodes * sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(r_d, r, numNodes * sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(p_d, p, numNodes * sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(Ap_d, Ap, numNodes * sizeof(double),cudaMemcpyHostToDevice);

  cudaMemcpy(face_d, face, numElems * sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(edofMat_d, edofMat, (numElems*MNE) * sizeof(int),cudaMemcpyHostToDevice);



  //define some values
  int iter = 0;
  const int Max_iters = 3000;


  r_usage.ru_maxrss = 0.0; //set memory to zero
  double t = get_time(); // get time stamp

    //create handle for Ddot on gpu
    cublasHandle_t handle;
    cublasCreate(&handle);

  //==============PCG Iterations Start Here=========================//
  for(iter=0; iter<Max_iters; iter++)
  {


    //reset result [Ap_d==0]
    cudaMemset(Ap_d, 0.0, numNodes*sizeof(double));



    //Apply Matrix-Free MatVec kernel [Loop over elements]

       //apply the matvec on gpuID=gpu[0,1,2,3]
    GPUMatVec<<<NBLOCKS,BLOCKSIZE>>>(ni, nj, Ke_d, p_d, Ap_d,numElems);
       //cudaDeviceSynchronize();


    //Apply Diriclhlet BC [U==r==0] on faces [Loop over elements]
    DirichletBC<<<NBLOCKS,BLOCKSIZE>>>(p_d, Ap_d, numElems, edofMat_d, face_d);
    //cudaDeviceSynchronize();



    //Reset pAp value
    pAp_d[0] = 0.0;
    //Apply Ddot
    cublasDdot(handle, numNodes, p_d, 1, Ap_d, 1, pAp_d);
    //cudaDeviceSynchronize();



    // control alpha
    alpha = rzold/pAp_d[0];
    //printf("rzold = %1.12f pAp = %1.12f   alpha = %1.12f\n",rzold,pAp_d[0], alpha);



    //Update solution vectors[U,r]
    //daxpy (+alpha)
    cublasDaxpy(handle, numNodes, &( alpha), p_d, 1, U_d, 1); //y = ax + y
    //daxpy (-alpha)
    malpha = -alpha;
    cublasDaxpy(handle, numNodes, &(malpha), Ap_d, 1, r_d, 1); //y = -ax + y
    //cudaDeviceSynchronize();



    //Reset Norm
    rNorm_d[0] = 0.0;
    //Apply Dnorm
    cublasDnrm2(handle, numNodes, r_d, 1, rNorm_d);
    //cudaDeviceSynchronize();


    //Check residual norm and my alpha
    //printf("iter = %i, rNorm = %1.16lf, alpha = %lf \n", iter, rNorm_d[0], alpha);



    // converge if reached Tolerance 1e-6
    if(rNorm_d[0] <= TOL) break;




    //Update rznew
    rznew_d[0] = 0.0;
    cublasDdot(handle, numNodes, r_d, 1, r_d, 1, rznew_d);
    //cudaDeviceSynchronize();



    // control beta
    beta = rznew_d[0]/rzold;



    //update vector [p] for next iters...
    UpdateVec<<<NBLOCKS2,BLOCKSIZE>>>(p_d, r_d, beta, numNodes);
    //cudaDeviceSynchronize();



    //update rzold for next iters...
    rzold = rznew_d[0];


  }
    //kill handle
    cublasDestroy(handle);



  // get memory usage
  getrusage(RUSAGE_SELF, &r_usage);

  // StdOut prints
  printf("\n");
  printf("*** PCG Solver Converged after %d total iterations ***\n", iter+1);
  printf("\n");
  printf("Residual Reduction: [%1.4e] <== Tolerance [%1.4e]\n",rNorm_d[0],TOL);
  printf("\n");
  printf("Solver Statistics [Wall Time, Memory Usage]\n");
  // Compute Elapsed Time for accumulative iters
  printf("-------------------------------------------\n");
  printf("|   Elapsed Time: [%g sec]                 \n",get_time()-t);
  printf("|   Elapsed Time: [%g sec/iter]            \n",(get_time()-t)/(iter+1));
  printf("-------------------------------------------\n");
  printf("|  Memory Usage: [%lf MB]                  \n", r_usage.ru_maxrss*0.001);
  printf("-------------------------------------------\n");



  // Post-Processing

  //Send solution [U] back to host
  cudaMemcpy(U, U_d, numNodes * sizeof(double), cudaMemcpyDeviceToHost);
  for (int n=0; n<numNodes; n++){
     if (U[n] < 0) U[n] = -U[n];
  }

  // Write to CSV file to visualize with paraview
  FILE *solution;
  solution = fopen("solutions.csv", "w");
  fprintf(solution, "i,j,k,u\n");

  int ni2=ni+1;
  int nj2=nj+1;
  double max=-1;
  for (k=0; k<(nk+1); k++) {
    for (j=0; j<(nj+1); j++) {
      for (i=0; i<(ni+1); i++) {
          i0=I3D(ni2,nj2,i,j,k);
                  fprintf(solution,"%d,%d,%d,%lg\n",i,j,k, U[i0]);
                  if(U[i0]>max) max=U[i0];
              }
            }
          }
    printf("Tmax=%lg\n",max);
  fclose(solution);


  //free gpu memory
  cudaFree(U_d);
  cudaFree(r_d);
  cudaFree(p_d);
  cudaFree(Ap_d);
  cudaFree(Ke_d);
  cudaFree(face_d);
  cudaFree(edofMat_d);
  cudaFree(pAp_d);
  cudaFree(rznew_d);
  cudaFree(rNorm_d);


  //free cpu memory
  free(U);
  free(b);
  free(r);
  free(Ke);
  free(p);
  free(Ap);
  //free(z);
  free(face);
  free(edofMat);
  //free(ndofMat);

  //Reset GPU Device
  cudaDeviceReset();

  return EXIT_SUCCESS;
}
