#include "vismgpu.h"

clStruct OCL;

clStruct::clStruct(void ){
	// cout<<"init clStruct"<<endl;
	platform = getPlatform();
	// device = getthisDevice(platform,devicecode);
	device = getDevice(platform, CL_DEVICE_TYPE_GPU);//get GPU

	parallel.device = 1; //alaways GPU
	parallel.context = getContext(device);

	parallel.command_queue = getCommandQueue(parallel.context, device);
	parallel.local_work_size = new size_t[3];
	parallel.global_work_size = new size_t[3];


	program = getProgram("/Users/zzirui/c25/src/cfangpu.cl", parallel.context, device, "");
	getLJcl = getKernel(program,"getLJ");
	getLJm2cl = getKernel(program,"getLJm2");
	getphicl = getKernel(program,"getphi");
	getphi2cl = getKernel(program,"getphi2");
	getLJsolcl = getKernel(program,"getLJsol");
	getoutsideGBkernel = getKernel(program,"getousideboxGB");
	getGBcomponent = getKernel(program,"getGBcomponent");
	getCFAcomp = getKernel(program,"getCFAcomp");
}

void setparallel(ParallelStruct &parallel, int *nx, int dim){
   int i, n, limit = static_cast<int>(round(log2(256.0))); //limit = 8

   parallel.dim = dim;
   if (parallel.device == 0)// for cpu
   {
      for (i = 0; i < parallel.dim; i++)
         parallel.local_work_size[i] = 1; //local_work_size = 1
   }
   else //for gpu
   {
      n = limit/parallel.dim; //n = 8/3 = 2
      for (i = 0; i < parallel.dim; i++)
         parallel.local_work_size[i] = static_cast<int>(round(exp(n*log(2.0)))); //local_work_size = 2^n = 4
   }
   for (i = 0; i < parallel.dim; i++)
      if ((nx[i]+1)%(parallel.local_work_size[i]) == 0) // if nx+1 is integer multiple of local_work_size
         parallel.global_work_size[i] = nx[i]+1; // then global_work_size is nx+1
      else
         parallel.global_work_size[i] = (static_cast<int>((nx[i]+1)/parallel.local_work_size[i])+1)*parallel.local_work_size[i];
         //global_work_size is interger multiple of local_work_size, larger than nx
}

void print_parallel(ParallelStruct &parallel){
	cout<<"local work size = ";
	for (int i = 0; i < parallel.dim; i++){
		cout<<parallel.local_work_size[i]<<" ";
 	}
 	cout<<endl;

 	cout<<"global work size = ";
	for (int i = 0; i < parallel.dim; i++){
		cout<<parallel.global_work_size[i]<<" ";
 	}
 	cout<<endl;
}

#if 0
// parallel initialize surface of single solute
void ocl_initsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid, bool tightfit, clStruct &clstruct){
	int i, r, tindex[grid.dim];
	double x[grid.dim];
	int thearg;
	cl_int err;

	clstruct.parallel.global_work_size[0] = grid.nx[0]+1;
	clstruct.parallel.global_work_size[1] = grid.nx[1]+1;
	clstruct.parallel.global_work_size[2] = grid.nx[2]+1;
	clstruct.parallel.local_work_size = NULL;
	// double start_time = wtime();

	cl_mem phi = clCreateBuffer(clstruct.parallel.context, CL_MEM_READ_WRITE,
	(grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)*sizeof(double),
	NULL,&err);
	errchk("ocl_initsurf:clCreateBuffer",err);
	if (tightfit) {	// tight fit
		for (i = 0; i < sol.dnum; i++){
			thearg = 0;
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(cl_mem),&phi);
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(int),&i);
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(sol.d[i][0]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(sol.d[i][1]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(sol.d[i][2]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(sol.sigma[i]));//average with water done in kernel
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(int),&(grid.nx[0]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(int),&(grid.nx[1]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(int),&(grid.nx[2]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(grid.dx[0]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(grid.dx[1]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(grid.dx[2]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(grid.a[0]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(grid.a[1]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clSetKernelArg(clstruct.getphicl,thearg,sizeof(double),&(grid.a[2]));
			thearg++;
			errchk("clSetKernelArg",err);
			err = clEnqueueNDRangeKernel(clstruct.parallel.command_queue, clstruct.getphicl, grid.dim,NULL,
										 clstruct.parallel.global_work_size,
										 clstruct.parallel.local_work_size,0,NULL,NULL);
			errchk("ocl_initsurf:clEnqueueNDRangeKernel",err);
			clFinish(clstruct.parallel.command_queue);
		}
	readbufferdouble(surf.phi,phi,grid.nx[0],grid.nx[1],grid.nx[2],
					  clstruct.parallel.command_queue);
	}
	else{ // loose fit
		for (i = 0; i < grid.dim; i++)
			tindex[i] = 0;
		while (tindex[0] <= grid.nx[0]){
			setvalarray(surf.phi,tindex,-1.0);
			(tindex[grid.dim-1])++;
			for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--){
				tindex[i] = 0;
				(tindex[i-1])++;
			}
		}
	}
	// double run_time = wtime() - start_time;
	// cout<<run_time<<" seconds."<<endl;
	clReleaseMemObject(phi);
}
#endif
// #if 0
// calculate LJ field of sol in parallel 
void parLJ(SoluteStruct &sol, double*** LJ, GridData &grid){
	// cout << "parallel FORMING LJ" << endl;
	// double start_time = wtime();
	setparallel(OCL.parallel, grid.nx, grid.dim);
	cl_int err;
	//opencl memory object
	cl_mem LJmem = clCreateBuffer(OCL.parallel.context,CL_MEM_READ_WRITE,
							  (grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)*sizeof(double),
							  NULL,&err);
	errchk("parLJ:clCreateBuffer", err);
	for (int i = 0; i < sol.dnum; i++)
	{
		int thearg = 0;
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(cl_mem),&LJmem);
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(int),&i);
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(sol.d[i][0]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(sol.d[i][1]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(sol.d[i][2]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(sol.epsilon[i]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(sol.sigma[i]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(int),&(grid.nx[0]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(int),&(grid.nx[1]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(int),&(grid.nx[2]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(grid.dx[0]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(grid.dx[1]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(grid.dx[2]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(grid.a[0]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(grid.a[1]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getLJcl,thearg,sizeof(double),&(grid.a[2]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clEnqueueNDRangeKernel(OCL.parallel.command_queue, OCL.getLJcl, grid.dim, NULL,
									OCL.parallel.global_work_size,
									OCL.parallel.local_work_size,
									0, NULL, NULL);
		errchk("parLJ:clEnqueueNDRangeKernel",err);
		clFinish(OCL.parallel.command_queue);
	}

	readbufferdouble(LJ, LJmem, grid.nx[0], grid.nx[1], grid.nx[2],OCL.parallel.command_queue);

	// double run_time = wtime() - start_time;
	// cout<<run_time<<" seconds."<<endl;
	clReleaseMemObject(LJmem);
}


// read data from cl_mem to 3d array
void readbufferdouble(double ***the3darray, cl_mem &thebuffer, int xdim, int ydim,
                      int zdim, cl_command_queue &command_queue)
{
   cl_int err;
   double *thearray = new double[(xdim+1)*(ydim+1)*(zdim+1)];

   err = clEnqueueReadBuffer(command_queue,thebuffer,CL_TRUE,0,
                             (zdim+1)*(ydim+1)*(xdim+1)*sizeof(double),thearray,0,NULL,
                             NULL);
   errchk("readbufferdouble:clEnqueueReadBuffer",err);
   for (int k = 0; k <= zdim; k++)
      for (int j = 0; j <= ydim; j++)
         for (int i = 0; i <= xdim; i++)
            the3darray[i][j][k] = thearray[(ydim+1)*(xdim+1)*k+(xdim+1)*j+i];

   delete [] thearray;
}


// helper for parallel GB, read buffer to 4darray CFA component
void readbufferdouble4d(double ****the4darray, cl_mem &thebuffer, int xdim, int ydim,
							 int zdim, cl_command_queue &command_queue)
{
	cl_int err;
	double *thearray = new double[(xdim+1)*(ydim+1)*(zdim+1)*3];

	err = clEnqueueReadBuffer(command_queue,thebuffer,CL_TRUE,0,
									  3*(zdim+1)*(ydim+1)*(xdim+1)*sizeof(double),thearray,0,NULL,
									  NULL);

	errchk("readbufferdouble4d:clEnqueueReadBuffer",err);
	for (int k = 0; k <= zdim; k++)
		for (int j = 0; j <= ydim; j++)
			for (int i = 0; i <= xdim; i++)
				for (int d = 0; d < 3; d++)
					the4darray[i][j][k][d] = thearray[ d*(xdim+1)*(ydim+1)*(zdim+1)+ (ydim+1)*(xdim+1)*k+ (xdim+1)*j + i];

	delete [] thearray;
}

// write solute data as 1 D, average with water
double* solute1d(SoluteStruct &psol, bool average){
	const int dim = 3;
	int N = psol.dnum;
	double *temp = new double[(dim+3)*N];//[x 1-N, y 1-N, z 1-N, epsilon 1-N, sigma 1-N, Q 1-N]
	for (int i = 0; i < N; i++){
		for (int j = 0; j < dim; j++){
			temp[j*N+i] = psol.d[i][j];
		}
	}
	double sigma, epsilon;
	for (int i = 0; i < N; i++)
	{
		if(average){
			sigma = (psol.sigma[i]+WATERSIG)/2;
			epsilon = sqrt(psol.epsilon[i]*WATEREPS);	
		}else{
			sigma = psol.sigma[i];
			epsilon = psol.epsilon[i];	
		}

	  temp[dim*N+i] = epsilon;
	  temp[(dim+1)*N+i] = sigma;
	  temp[(dim+2)*N+i] = psol.Q[i];
	}
	return temp;
}

// for each atom in large protein, calculate it's LJ with all atom of small protein in parallel
// psol is the smaller one, N is num of psol
double ocl_LJsolsol(SoluteStruct &sol,SoluteStruct &psol){
	cl_int err;
	const int dim = 3;
	int N = psol.dnum;

	OCL.parallel.dim = 1;
	OCL.parallel.global_work_size[0] = N;
	OCL.parallel.global_work_size[1] = 0;
	OCL.parallel.global_work_size[2] = 0;
	
	cl_mem psolmem = clCreateBuffer(OCL.parallel.context,CL_MEM_READ_WRITE,(dim+3)*N*sizeof(double),NULL,&err);//memery for data of psol
	errchk("clCreateBuffer", err);
	double *temp = solute1d(psol, false);

	err = clEnqueueWriteBuffer(OCL.parallel.command_queue, psolmem, CL_TRUE,0,
	                         (dim+2)*N*sizeof(double),temp,0,NULL,
	                         NULL);

	int thearg;
	cl_mem LJsolsolmem = clCreateBuffer(OCL.parallel.context, CL_MEM_READ_WRITE, N*sizeof(double),NULL,&err);
	errchk("clCreateBuffer", err);

	for (int i = 0; i < sol.dnum; i++)
	{
	  thearg = 0;
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(cl_mem),&LJsolsolmem);
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(cl_mem),&psolmem);
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(int),&i);
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(double),&(sol.d[i][0]));
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(double),&(sol.d[i][1]));
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(double),&(sol.d[i][2]));
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(double),&(sol.epsilon[i]));
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(double),&(sol.sigma[i]));
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clSetKernelArg(OCL.getLJsolcl,thearg,sizeof(int),&N);
	  thearg++;
	  errchk("clSetKernelArg",err);
	  err = clEnqueueNDRangeKernel(OCL.parallel.command_queue,OCL.getLJsolcl,OCL.parallel.dim, NULL,
	  							OCL.parallel.global_work_size,
	                               NULL,0,NULL,NULL);
	  errchk("clEnqueueNDRangeKernel",err);
	}

	double *thesum = new double[N];

	err = clEnqueueReadBuffer(OCL.parallel.command_queue,LJsolsolmem,CL_TRUE,0,
	                         N*sizeof(double),thesum,0,NULL,NULL);
	errchk("clEnqueueReadBuffer",err);
	double value = 0.0;
	for (int i = 0; i < N; i++)
	  value += thesum[i];

	delete [] thesum;
	delete [] temp;
	clReleaseMemObject(psolmem);
	clReleaseMemObject(LJsolsolmem);
	return value;
}

// parallel GB
void parCFA(SoluteStruct &sol, double*** GB, GridData &grid){
	// double start_time = wtime();
	cl_int err;
	
	double ****CFAcomp = matrix4d<double>(grid.nx[0],grid.nx[1],grid.nx[2]);
	parCFAcomp2(sol, CFAcomp, grid);

	for (int i = 0; i <= grid.nx[0]; i++){
		for (int j = 0; j <= grid.nx[1]; j++){
			for (int k = 0; k <= grid.nx[2]; k++){
				GB[i][j][k] = 
				1.0/(32.0*M_PI*M_PI*sol.epsilon0)*(1.0/sol.epsilonex-1.0/sol.epsilonin)*
				(CFAcomp[i][j][k][0]*CFAcomp[i][j][k][0] +CFAcomp[i][j][k][1]*CFAcomp[i][j][k][1] + CFAcomp[i][j][k][2]*CFAcomp[i][j][k][2]);
			}
		}
	}

	free_matrix4d<double>(CFAcomp,grid.nx[0],grid.nx[1],grid.nx[2]);
}

// parallel CFA component
void parCFAcomp(SoluteStruct &sol, double**** CFAcomp, GridData &grid){
	// double start_time = wtime();
	cl_int err;
	setparallel(OCL.parallel, grid.nx, grid.dim);
	// v is a vecotr field as 4D array, 
	// v = sum of Qi (x-xi)/|x-xi|^3, v[0/1/2]=component in x/y/z
	// GB = 1/(32 pi^2 epsilon0) (1/epsilonex - 1/epsilonin) v^2
	cl_mem v = clCreateBuffer(OCL.parallel.context, CL_MEM_READ_WRITE,
							 3*(grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)*sizeof(double),
							 NULL,&err);
	errchk("parCFA:clCreateBuffer", err);
	//for each particle, calculate LJ in parallel on the grid
	// work-item: LJ on 1 grid point
	for (int i = 0; i < sol.dnum; i++)//for each particle, 
	{
		int thearg = 0;
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(cl_mem), &v);
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(int), &i);
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(sol.d[i][0]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(sol.d[i][1]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(sol.d[i][2]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(sol.Q[i]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(int), &(grid.nx[0]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(int), &(grid.nx[1]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(int), &(grid.nx[2]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(grid.dx[0]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(grid.dx[1]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(grid.dx[2]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(grid.a[0]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(grid.a[1]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clSetKernelArg(OCL.getGBcomponent, thearg, sizeof(double), &(grid.a[2]));
		thearg++;
		errchk("clSetKernelArg",err);
		err = clEnqueueNDRangeKernel(OCL.parallel.command_queue, OCL.getGBcomponent, grid.dim, NULL,
			OCL.parallel.global_work_size, OCL.parallel.local_work_size,
			0,NULL,NULL);
		errchk("clEnqueueNDRangeKernel",err);
		clFinish(OCL.parallel.command_queue);
	}
	// read data from buffer to  master
	readbufferdouble4d(CFAcomp, v, grid.nx[0],grid.nx[1],grid.nx[2],OCL.parallel.command_queue);
	// double run_time = wtime() - start_time;
	// cout<<run_time<<" seconds."<<endl;
	// print_sol(sol,grid);
	clReleaseMemObject(v);
}

//initialize parallel solute
void parLJm2(SoluteStruct &sol, double*** LJ, GridData &grid){

	const int dim = 3;
	cl_int err;
	int N = sol.dnum;

	OCL.parallel.global_work_size[0] = grid.nx[0]+1;
	OCL.parallel.global_work_size[1] = grid.nx[1]+1;
	OCL.parallel.global_work_size[2] = grid.nx[2]+1;
	// OCL.parallel.local_work_size = NULL;	
	
	cl_mem solmem = clCreateBuffer(OCL.parallel.context, CL_MEM_READ_WRITE,(dim+3)*N*sizeof(double),NULL,&err);//memery for data of sol
	errchk("parLJm2:clCreateBuffer", err);
	double *sol1d = solute1d(sol,true);
	err = clEnqueueWriteBuffer(OCL.parallel.command_queue, solmem, CL_TRUE,0,
	                         (dim+3)*N*sizeof(double),sol1d,0,NULL,
	                         NULL);
	//opencl memory object, 1D array for LJ field
	cl_mem LJmem = clCreateBuffer(OCL.parallel.context,CL_MEM_READ_WRITE,
							  (grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)*sizeof(double),
							  NULL,&err);
	errchk("parLJm2:clCreateBuffer", err);
	
	int thearg = 0;
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(cl_mem),&LJmem);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(cl_mem),&solmem);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(int),&N);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(int),&(grid.nx[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(int),&(grid.nx[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(int),&(grid.nx[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.dx[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.dx[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.dx[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.a[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.a[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.a[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clEnqueueNDRangeKernel(OCL.parallel.command_queue, OCL.getLJm2cl, dim, NULL,
								OCL.parallel.global_work_size,
								NULL,
								0, NULL, NULL);
	errchk("parLJm2:clEnqueueNDRangeKernel",err);
	clFinish(OCL.parallel.command_queue);

	readbufferdouble(LJ, LJmem, grid.nx[0], grid.nx[1], grid.nx[2], OCL.parallel.command_queue);

	// double run_time = wtime() - start_time;
	// cout<<run_time<<" seconds."<<endl;
	delete [] sol1d;
	clReleaseMemObject(solmem);
	clReleaseMemObject(LJmem);

}


//initialize parallel solute
void ocl_initsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid, bool tightfit){

	const int dim = 3;
	cl_int err;
	int N = sol.dnum;

	OCL.parallel.global_work_size[0] = grid.nx[0]+1;
	OCL.parallel.global_work_size[1] = grid.nx[1]+1;
	OCL.parallel.global_work_size[2] = grid.nx[2]+1;
	// OCL.parallel.local_work_size = NULL;
	
	cl_mem solmem = clCreateBuffer(OCL.parallel.context,CL_MEM_READ_WRITE,(dim+3)*N*sizeof(double),NULL,&err);//memery for data of sol
	errchk("ocl_initsurf:clCreateBuffer", err);
	double *sol1d = solute1d(sol,true);
	err = clEnqueueWriteBuffer(OCL.parallel.command_queue, solmem, CL_TRUE,0,
	                         (dim+3)*N*sizeof(double),sol1d,0,NULL,
	                         NULL);
	//opencl memory object, 1D array for LJ field
	cl_mem phimem = clCreateBuffer(OCL.parallel.context,CL_MEM_READ_WRITE,
							  (grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)*sizeof(double),
							  NULL,&err);
	errchk("ocl_initsurf:clCreateBuffer", err);
	int tightfitint = (tightfit?1:0);
	int thearg = 0;
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(cl_mem),&phimem);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(cl_mem),&solmem);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(int),&tightfitint);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(int),&N);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(int),&(grid.nx[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(int),&(grid.nx[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(int),&(grid.nx[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(double),&(grid.dx[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(double),&(grid.dx[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(double),&(grid.dx[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(double),&(grid.a[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(double),&(grid.a[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getphi2cl,thearg,sizeof(double),&(grid.a[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clEnqueueNDRangeKernel(OCL.parallel.command_queue, OCL.getphi2cl, grid.dim, NULL, OCL.parallel.global_work_size, NULL, 0, NULL, NULL);
	errchk("ocl_initsurf:clEnqueueNDRangeKernel",err);
	clFinish(OCL.parallel.command_queue);

	readbufferdouble(surf.phi, phimem, grid.nx[0], grid.nx[1], grid.nx[2], OCL.parallel.command_queue);

	// double run_time = wtime() - start_time;
	// cout<<run_time<<" seconds."<<endl;
	delete [] sol1d;
	clReleaseMemObject(solmem);
	clReleaseMemObject(phimem);


}


// parallel CFA component
void parCFAcomp2(SoluteStruct &sol, double**** CFAcomp, GridData &grid){
	OCL.parallel.global_work_size[0] = grid.nx[0]+1;
	OCL.parallel.global_work_size[1] = grid.nx[1]+1;
	OCL.parallel.global_work_size[2] = grid.nx[2]+1;
	// OCL.parallel.local_work_size = NULL;
	
	const int dim = 3;
	cl_int err;
	int N = sol.dnum;
	cl_mem solmem = clCreateBuffer(OCL.parallel.context,CL_MEM_READ_WRITE,6*N*sizeof(double),NULL,&err);//memery for data of sol
	errchk("parCFAcompm2:clCreateBuffer", err);

	double *sol1d = solute1d(sol,true);
	err = clEnqueueWriteBuffer(OCL.parallel.command_queue, solmem, CL_TRUE,0,6*N*sizeof(double), sol1d, 0, NULL, NULL);
	errchk("parCFAcompm2:clEnqueueWriteBuffer", err);

	cl_mem cl_CFAcomp = clCreateBuffer(OCL.parallel.context, CL_MEM_READ_WRITE,
							 3*(grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)*sizeof(double),
							 NULL,&err);
	errchk("parCFAcompm2:clCreateBuffer", err);

	int thearg = 0;
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(cl_mem), &cl_CFAcomp);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(cl_mem), &solmem);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp,thearg, sizeof(int),&N);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(int), &(grid.nx[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(int), &(grid.nx[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(int), &(grid.nx[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(double), &(grid.dx[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(double), &(grid.dx[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(double), &(grid.dx[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(double), &(grid.a[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(double), &(grid.a[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getCFAcomp, thearg, sizeof(double), &(grid.a[2]));
	thearg++;
	errchk("clSetKernelArg",err);

	err = clEnqueueNDRangeKernel(OCL.parallel.command_queue, OCL.getCFAcomp, grid.dim, NULL, OCL.parallel.global_work_size, NULL,0,NULL,NULL);
	errchk("parCFAcompm2:clEnqueueNDRangeKernel",err);
	clFinish(OCL.parallel.command_queue);
	
	readbufferdouble4d(CFAcomp, cl_CFAcomp, grid.nx[0], grid.nx[1], grid.nx[2],OCL.parallel.command_queue);
	delete [] sol1d;
	clReleaseMemObject(solmem);
	clReleaseMemObject(cl_CFAcomp);
}


// test integration of 1/rho^4 on ball radius 1 outside
double paroutsidebox(int ntheta, int nphi, int nr, GridData &grid){

	int N = (grid.nx[0])*(grid.nx[1])*(grid.nx[2]);
	size_t global_work_size[3];
	global_work_size[0] = ntheta;
	global_work_size[1] = nphi;
	global_work_size[2] = nr;
	cl_int err;
	cl_mem cl_field = clCreateBuffer(OCL.parallel.context, CL_MEM_READ_WRITE,
							 N*sizeof(double),
							 NULL,&err);
	errchk("paroutsidebox:clCreateBuffer", err);
	cl_kernel integrate = getKernel(OCL.program,"ousidebox_integrate");

	err = clSetKernelArg(integrate, 0, sizeof(cl_mem), &cl_field);
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 1, sizeof(int), &(grid.nx[0]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 2, sizeof(int), &(grid.nx[1]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 3, sizeof(int), &(grid.nx[2]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 4, sizeof(double), &(grid.dx[0]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 5, sizeof(double), &(grid.dx[1]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 6, sizeof(double), &(grid.dx[2]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 7, sizeof(double), &(grid.a[0]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 8, sizeof(double), &(grid.a[1]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(integrate, 9, sizeof(double), &(grid.a[2]));
	errchk("clSetKernelArg",err);

	err = clEnqueueNDRangeKernel(OCL.parallel.command_queue, integrate, grid.dim, NULL,
		global_work_size, NULL, // default local work size
		0,NULL,NULL);
	errchk("paroutsidebox:clEnqueueNDRangeKernel",err);
	clFinish(OCL.parallel.command_queue);
	
	double *thesum = new double[N];
	err = clEnqueueReadBuffer(OCL.parallel.command_queue, cl_field, CL_TRUE, 0,
	                         N*sizeof(double),thesum,0,NULL,NULL);
	errchk("paroutsidebox:clEnqueueReadBuffer",err);
	double value = 0.0;
	for (int i = 0; i < N; i++){
		// cout<<thesum[i]<<endl;
		value += thesum[i];
	}
	delete [] thesum;

	clReleaseMemObject(cl_field);
	return value;
}

// test integration of 1/rho^4 on ball radius 1 outside
double ocl_outsidebox_cfa(SoluteStruct &sol, int ntheta, int nphi, int nr, GridData &boxgrid){

	GridData sphgrid(ntheta, nphi, nr, 0, 0, 0, 2*M_PI, M_PI, 1.0/boxgrid.b[0]);
	int N = (sphgrid.nx[0])*(sphgrid.nx[1])*(sphgrid.nx[2]);
	size_t global_work_size[3];
	global_work_size[0] = ntheta;
	global_work_size[1] = nphi;
	global_work_size[2] = nr;
	
	double boxax, boxay, boxaz;// boundary of box 
	boxax = boxgrid.a[0] - boxgrid.dx[0]/2.0;
	boxay = boxgrid.a[1] - boxgrid.dx[1]/2.0;
	boxaz = boxgrid.a[2] - boxgrid.dx[2]/2.0;

	cl_int err;
	cl_mem solmem = clCreateBuffer(OCL.parallel.context,CL_MEM_READ_WRITE,6*sol.dnum*sizeof(double),NULL,&err);//memery for data of sol
	errchk("ocl_outsidebox_cfa:clCreateBuffer", err);

	double *sol1d = solute1d(sol,true);
	err = clEnqueueWriteBuffer(OCL.parallel.command_queue, solmem, CL_TRUE,0,6*sol.dnum*sizeof(double), sol1d, 0, NULL, NULL);
	errchk("ocl_outsidebox_cfa:clEnqueueWriteBuffer", err);


	cl_mem cl_field = clCreateBuffer(OCL.parallel.context, CL_MEM_READ_WRITE,
							 N*sizeof(double),
							 NULL,&err);
	errchk("ocl_outsidebox_cfa:clCreateBuffer", err);
	

	err = clSetKernelArg(OCL.getoutsideGBkernel, 0, sizeof(cl_mem), &cl_field);
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 1, sizeof(cl_mem), &solmem);
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 2, sizeof(int), &sol.dnum);
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 3, sizeof(int), &(sphgrid.nx[0]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 4, sizeof(int), &(sphgrid.nx[1]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 5, sizeof(int), &(sphgrid.nx[2]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 6, sizeof(double), &(sphgrid.dx[0]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 7, sizeof(double), &(sphgrid.dx[1]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 8, sizeof(double), &(sphgrid.dx[2]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 9, sizeof(double), &(sphgrid.a[0]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 10, sizeof(double), &(sphgrid.a[1]));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 11, sizeof(double), &(sphgrid.a[2]));
	errchk("clSetKernelArg",err);	
	err = clSetKernelArg(OCL.getoutsideGBkernel, 12, sizeof(double), &(boxax));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 13, sizeof(double), &(boxay));
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getoutsideGBkernel, 14, sizeof(double), &(boxaz));
	errchk("clSetKernelArg",err);

	err = clEnqueueNDRangeKernel(OCL.parallel.command_queue, OCL.getoutsideGBkernel, sphgrid.dim, NULL,
		global_work_size, NULL, // default local work size
		0,NULL,NULL);
	errchk("ocl_outsidebox_cfa:clEnqueueNDRangeKernel",err);
	clFinish(OCL.parallel.command_queue);
	
	double coeff = 1.0/(32.0*M_PI*M_PI*sol.epsilon0)*(1.0/sol.epsilonex-1.0/sol.epsilonin);
	double value = 0.0;

	double *thesum = new double[N];
	err = clEnqueueReadBuffer(OCL.parallel.command_queue, cl_field, CL_TRUE, 0,
	                         N*sizeof(double),thesum,0,NULL,NULL);
	errchk("ocl_outsidebox_cfa:clEnqueueReadBuffer",err);
	
	for (int i = 0; i < N; i++){
		// cout<<thesum[i]<<endl;
		value += thesum[i];
	}
	value *= coeff;

	delete[] thesum;
	clReleaseMemObject(solmem);
	clReleaseMemObject(cl_field);
	return value;
}

// overload 
double ocl_outsidebox_cfa(SoluteStruct &sol, GridData &boxgrid){
	const int ntheta = 50;
	const int nphi = 50;
	const int nr = 50;
	return ocl_outsidebox_cfa(sol, ntheta, nphi, nr, boxgrid);
}


//initialize parallel solute
void parLJ1d(SoluteStruct &sol, double* LJ1d, GridData &grid){
	// cout << "parallel FORMING LJ" << endl;
	// double start_time = wtime();
	OCL.parallel.global_work_size[0] = grid.nx[0]+1;
	OCL.parallel.global_work_size[1] = grid.nx[1]+1;
	OCL.parallel.global_work_size[2] = grid.nx[2]+1;
	const int dim = 3;
	cl_int err;
	int N = sol.dnum;
	
	cl_mem solmem = clCreateBuffer(OCL.parallel.context,CL_MEM_READ_ONLY,(dim+3)*N*sizeof(double),NULL,&err);//memery for data of sol
	errchk("parLJm2:clCreateBuffer", err);
	double *sol1d = solute1d(sol,true);
	err = clEnqueueWriteBuffer(OCL.parallel.command_queue, solmem, CL_TRUE,0,
	                         (dim+2)*N*sizeof(double),sol1d,0,NULL,
	                         NULL);
	//opencl memory object, 1D array for LJ field
	cl_mem LJmem = clCreateBuffer(OCL.parallel.context,CL_MEM_WRITE_ONLY,
							  (grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)*sizeof(double),
							  NULL,&err);
	errchk("parLJm2:clCreateBuffer", err);
	
	int thearg = 0;
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(cl_mem),&LJmem);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(cl_mem),&solmem);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(int),&N);
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(int),&(grid.nx[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(int),&(grid.nx[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(int),&(grid.nx[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.dx[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.dx[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.dx[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.a[0]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.a[1]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clSetKernelArg(OCL.getLJm2cl,thearg,sizeof(double),&(grid.a[2]));
	thearg++;
	errchk("clSetKernelArg",err);
	err = clEnqueueNDRangeKernel(OCL.parallel.command_queue, OCL.getLJm2cl, grid.dim, NULL,
								OCL.parallel.global_work_size,NULL,0, NULL, NULL);
	errchk("parLJm2:clEnqueueNDRangeKernel",err);
	clFinish(OCL.parallel.command_queue);

   err = clEnqueueReadBuffer(OCL.parallel.command_queue,LJmem,CL_TRUE,0,(grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)*sizeof(double),LJ1d,0,NULL,NULL);
   errchk("readbufferdouble:clEnqueueReadBuffer",err);
 
	delete [] sol1d;
	clReleaseMemObject(solmem);
	clReleaseMemObject(LJmem);
}
