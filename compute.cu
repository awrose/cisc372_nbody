
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "vector.h"
#include "config.h"

__global__ void computePairwiseAccels(vector3** accels, vector3* hPos, double* mass){
	int i = (blockDim.y * blockIdx.y) + threadIdx.y;
	int j = (blockDim.x * blockIdx.x) + threadIdx.x;

	if(i < NUMENTITIES && j < NUMENTITIES){
		if(i == j){
			FILL_VECTOR(accels[i][j],0,0,0);
		}else{
			vector3 distance;
			for (int k=0;k<3;k++) distance[k]=hPos[i][k]-hPos[j][k];
			double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
			double magnitude=sqrt(magnitude_sq);
			double accelmag=-1*GRAV_CONSTANT*mass[j]/magnitude_sq;
			FILL_VECTOR(accels[i][j],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);			
		}
	}
}

__global__ void computeSum(vector3** accels, vector3* hPos, vector3* hVel){
	int in = blockIdx.x;
	if(in < NUMENTITIES){
		vector3 accel_sum={0,0,0};
		for (int j=0;j<NUMENTITIES;j++){
			for (int k=0;k<3;k++)
				//printf("screaming crying throwing up");
				accel_sum[k]+=accels[in][j][k];
		}
		//compute the new velocity based on the acceleration and time interval
		//compute the new position based on the velocity and time interval
		for (int k=0;k<3;k++){
			hVel[in][k]+=accel_sum[k]*INTERVAL;
			hPos[in][k]+=hVel[in][k]*INTERVAL;
		}
	}

}

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
void compute(){

	dim3 block_size(16,16);
	dim3 block_count((NUMENTITIES+15) / block_size.x, (NUMENTITIES+15) / block_size.y);
	computePairwiseAccels<<<block_count, block_size>>>(accels, d_hPos, d_mass);
	cudaDeviceSynchronize();

	dim3 grid_dims (NUMENTITIES, 1, 1);
	computeSum<<<grid_dims, 3>>>(accels, d_hPos, d_hVel);
	cudaDeviceSynchronize();

}