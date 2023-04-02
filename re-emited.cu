/* 
# Made by:
# Sergio Mendoza <sergio@mendozza.org>
# Milton Santibañez <msantibanez@astro.unam.mx>
# Gustavo Magallanes-Guijon <gustavo.magallanes.guijon@ciencias.unam.mx>
# Instituto de Astronomia UNAM
# Ciudad Universitaria
# Ciudad de Mexico
# Mexico
# Fri 21 Oct 2020 05:40:36 PM UTC
*/

//includes
#include "content/Headers.h"

// defines
#define COLUMNS 3 
#define LADO 400
#define ROWS LADO*LADO

#define MATRIX COLUMNS*ROWS

//Grid
#define THREADS 128 
#define BLOCKS (int)ceil(ROWS/THREADS) + 1


/*datafile to be open for write*/
FILE *datafilewrite;

/*datafile to be open for read*/
FILE *datafileread;

//////////////////////////////////////////////////////////////////////////////
////////////////////////  K E R N E L  ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
__global__ 
void remited(int k, 
	     int steps, 
	     double *dev_vector_emisor, 
	     double *dev_vector_receptor, 
	     double d_min, 
	     double l_min, 
	     double d_max, 
	     double l_max)
{

    double x_grid, y_grid, z_grid,  diff_y, diff_z;
  
    double Px, Py, Pz, Tau_x, Tau_y, Tau_z;
    double r_evol, theta_evol, phi_evol, r_pre, theta_pre, phi_pre, e_r, e_theta, e_phi, e_r_pre, e_theta_pre, e_phi_pre, dl;
    double distancia, r_orbita, dt, theta_0;
    int sentinel, fila;
  
    double epsilon = 0.1;
  
    sentinel = 0;
    distancia = 1000.0;
    theta_0 = 0.5*M_PI;
  
    dl = 0.01;
    r_orbita = 100.0;
    dt = M_PI/steps;

    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    //points of the grid which will be used like initial points in the evolution
    Px = r_orbita*cosf((0.5*M_PI)+k*dt);
    Py = dev_vector_emisor[(idx)*COLUMNS + 0] + r_orbita*sinf((0.5*M_PI)+k*dt);
    Pz = dev_vector_emisor[(idx)*COLUMNS + 1];

    //we obtain the initial values of Tau
    Tau_x = sinf(theta_0) * cosf(0.0);
    Tau_y = sinf(theta_0) * sinf(0.0);
    Tau_z = cosf(theta_0);

    //There are the initialization of the coordenates and the vector e for de leap frog algoritm
    r_evol = sqrt( Px*Px + Py*Py + Pz*Pz );
    if(Pz > 0.0)theta_evol = atanf(sqrt( Px*Px + Py*Py )/Pz);
    if(Pz == 0.0)theta_evol = 0.5*M_PI;
    if(Pz < 0.0)theta_evol = M_PI + atanf(sqrt( Px*Px + Py*Py )/Pz);
    if(Px > 0.0 && Py >= 0.0)phi_evol = atanf(Py/Px);
    if(Px > 0.0 && Py < 0.0)phi_evol = 2.0*M_PI + atanf(Py/Px);
    if(Px == 0.0)phi_evol = copysign(0.5*M_PI,Py);
    if(Px < 0.0)phi_evol = M_PI + atanf(Py/Px);

    e_r = Tau_x*sinf(theta_evol)*cosf(phi_evol) + Tau_y*sinf(theta_evol)*sinf(phi_evol) + Tau_z*cosf(theta_evol);
    e_theta = Tau_x*cosf(theta_evol)*cosf(phi_evol)/r_evol + Tau_y*cosf(theta_evol)*sinf(phi_evol)/r_evol - Tau_z*sinf(theta_evol)/r_evol;
    e_phi = -1.0*Tau_x*sinf(phi_evol)/(r_evol*sinf(theta_evol)) + Tau_y*cosf(phi_evol)/(r_evol*sinf(theta_evol));

    #pragma unroll 
    for (;;){
        //We rename the variables for the next calculus
        r_pre = r_evol; 
        theta_pre = theta_evol;
        phi_pre = phi_evol;
  
        e_r_pre = e_r;
        e_theta_pre = e_theta;
        e_phi_pre = e_phi;
  
        //We evolve the coordenates and the components of e, with the variables of the past step
        r_evol = r_pre + e_r_pre*dl;
        theta_evol = theta_pre + e_theta_pre*dl;
        phi_evol = phi_pre + e_phi_pre*dl;
  
        e_r = e_r_pre + dl*de_r(r_pre,theta_pre,phi_pre,e_r_pre,e_theta_pre,e_phi_pre);
        e_theta = e_theta_pre + dl*de_theta(r_pre,theta_pre,phi_pre,e_r_pre,e_theta_pre,e_phi_pre);
        e_phi = e_phi_pre + dl*de_phi(r_pre,theta_pre,phi_pre,e_r_pre,e_theta_pre,e_phi_pre); 
  
        //We write the next point in cartesian coordenates
        Px = r_evol*sinf(theta_evol)*cosf(phi_evol);
        Py = r_evol*sinf(theta_evol)*sinf(phi_evol);
        Pz = r_evol*cosf(theta_evol);
  
        x_grid = Px - r_orbita*cosf((0.5*M_PI)+k*dt);
        y_grid = Py - r_orbita*sinf((0.5*M_PI)+k*dt);
        z_grid = Pz;
  

	if( 0.5*distancia - epsilon <= x_grid && x_grid <= 0.5*distancia + epsilon
  		       && d_min <= y_grid && y_grid <= d_max
  		       && l_min <= z_grid && z_grid <= l_max && sentinel==0 ){
             sentinel = 1;
             diff_y = 1.0*LADO;
             diff_z = 1.0*LADO;

             fila = 0;
             #pragma unroll
	     for(int m = 0 ; m < ROWS ; m++){
          	  if(fabsf(dev_vector_receptor[m*COLUMNS + 0] - y_grid) <= diff_y && 
		     fabsf(dev_vector_receptor[m*COLUMNS + 1] - z_grid) <= diff_z ){
              	          diff_y = fabsf(dev_vector_receptor[m*COLUMNS + 0] - y_grid);
  		          diff_z = fabsf(dev_vector_receptor[m*COLUMNS + 1] - z_grid);
        	          fila = m;
            	  }//END IF
             }//END FOR

	atomicAdd(&dev_vector_receptor[fila*COLUMNS + 2], dev_vector_emisor[idx*COLUMNS + 2]);

__syncthreads();
	}//END IF

    if(!( r_evol<=(1.1*distancia) && r_evol>1.1 && sentinel==0)) break;
    }//END FOR 
    sentinel = 0;   
}

//////////////////////////////////////////////////////////////////////////////
////////////////////////////////  D E V I C E  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void onDevice(  double* hst_vector_emisor, 
		double* hst_vector_receptor, 
		double d_min, 
		double l_min, 
		double d_max, 
		double l_max  )
{
 // configuration Grids and Theads
 // dim3 GridBlocks(50,100);
 // dim3 ThreadsBlocks(8,4);


    //INIT KERNEL LOOP
    int steps = 180;

    //NAME OF FILE
    char nombre[20];

    // start timer
    GpuTimer timer;
    timer.Start();

    #pragma unroll steps
    for(int k=0; k<360; k++){

        double *dev_vector_emisor, *dev_vector_receptor;

       //cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

        // memory in the device
        HANDLER_ERROR_ERR(cudaMalloc((void**)&dev_vector_emisor, MATRIX * sizeof(double)));
        HANDLER_ERROR_ERR(cudaMalloc((void**)&dev_vector_receptor, MATRIX * sizeof(double)));
        
        // copy the data to the device
        HANDLER_ERROR_ERR(cudaMemcpy(dev_vector_emisor, hst_vector_emisor, MATRIX * sizeof(double), 
				cudaMemcpyHostToDevice));

	for(int i = 0; i < ROWS; i++){
	   hst_vector_receptor[i*3 + 2] = 0.0;
	  }

        // copy the data to the device
        HANDLER_ERROR_ERR(cudaMemcpy(dev_vector_receptor, hst_vector_receptor, MATRIX * sizeof(double), 
				cudaMemcpyHostToDevice));

        sprintf(nombre,"second_display_image%d.dat",k);
        
        datafilewrite = fopen ( nombre , "w" );

        ///////////////////////////////////////////////  K E R N E L  /////////////////////////////////////////////////////
        //remited<<<GridBlocks,ThreadsBlocks>>>(k, steps, dev_vector_emisor, dev_vector_receptor, d_min, l_min, d_max, l_max);
        remited<<<BLOCKS,THREADS>>>(k, steps, dev_vector_emisor, dev_vector_receptor, d_min, l_min, d_max, l_max);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        HANDLER_ERROR_MSG("Kernel Panic!!!");
        
        // data from the  host to device
        HANDLER_ERROR_ERR(cudaMemcpy(hst_vector_receptor, dev_vector_receptor, MATRIX  * sizeof(double), 
				cudaMemcpyDeviceToHost));

	#pragma unroll ROWS
        for (int i = 0; i<ROWS; i++){
            if(hst_vector_receptor[(i)*COLUMNS + 0] != 0 && hst_vector_receptor[(i)*COLUMNS + 1] != 0){
          		fprintf(datafilewrite, "%lE \t %lE \t %lE \n", hst_vector_receptor[(i)*COLUMNS + 0]* 5.029643510174233, \
           						               hst_vector_receptor[(i)*COLUMNS + 1]* 5.029643510174233, \
           							       hst_vector_receptor[(i)*COLUMNS + 2]);
            }
        }

        fclose(datafilewrite);

        // liberamos memoria del device
        HANDLER_ERROR_ERR(cudaFree( dev_vector_receptor ));
        HANDLER_ERROR_ERR(cudaFree( dev_vector_emisor ));
    } //FOR 
  
    timer.Stop();

    // print time
//    printf("Time :  %f ms\n", timer.Elapsed());
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////  H O S T ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void onHost()
{
  
    double y, z, intensity;
  
    double d_min, l_min, d_max, l_max;
  
    d_min = 0.0;
    d_max = 0.0;
    l_min = 0.0;
    l_max = 0.0;
  
    // declaration of vectors
    double *hst_vector_emisor, *hst_vector_receptor;
    
    // memory in the  host
    hst_vector_emisor = (double*)malloc(MATRIX * sizeof(double));
    hst_vector_receptor = (double*)malloc(MATRIX * sizeof(double));
    
    // init the input file
    datafileread = fopen ( "emited_image.dat" , "r" );
  
    #pragma unroll ROWS
    for(int i = 0; i < ROWS; i++){
        fscanf(datafileread,"%lE\t%lE\t%lE\n", &y, &z, &intensity);
        hst_vector_emisor[(i)*COLUMNS + 0] = y / 5.029643510174233;
        hst_vector_emisor[(i)*COLUMNS + 1] = z / 5.029643510174233;
        hst_vector_emisor[(i)*COLUMNS + 2] = intensity;
        hst_vector_receptor[(i)*COLUMNS + 0] = hst_vector_emisor[(i)*COLUMNS + 0];
        hst_vector_receptor[(i)*COLUMNS + 1] = hst_vector_emisor[(i)*COLUMNS + 1];

        if(d_max < hst_vector_emisor[(i)*COLUMNS + 0]) d_max = hst_vector_emisor[(i)*COLUMNS + 0];
        if(d_min > hst_vector_emisor[(i)*COLUMNS + 0]) d_min = hst_vector_emisor[(i)*COLUMNS + 0];
        if(l_max < hst_vector_emisor[(i)*COLUMNS + 1]) l_max = hst_vector_emisor[(i)*COLUMNS + 1];
        if(l_min > hst_vector_emisor[(i)*COLUMNS + 1]) l_min = hst_vector_emisor[(i)*COLUMNS + 1];
    }

    fclose(datafileread);

    // start timer
    CpuTimer timer;
    timer.Start();
  
    ///////////////////////////////  D E V I C E  ////////////////////////////////
    onDevice(hst_vector_emisor, hst_vector_receptor,  d_min, l_min, d_max, l_max);
    //////////////////////////////////////////////////////////////////////////////
  
    // stop timer
    timer.Stop();
    // print time
    //printf("CPU Time :  %f ms\n", timer.Elapsed());
  
  
    //salida del programa
    free(hst_vector_receptor);
    free(hst_vector_emisor);
//    printf("-: successful execution :-\n");
}

//////////////////////////////////////////////////////////////////////////////
/////////////////////////  M A I N ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int main() 
{
    onHost();
    return 0;
}

