#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

extern "C" {
	#include <athread.h>
    void slave_setHall(void);
}

#define MAX_LEN 30
#define MAX_VERT 1000
#define MAX_PARA 1000
#define mat_0_len 100
#define mat_1_len 100
#define mat_2_len 100
#define halo_size 1

double local_mat_old_old[local_size_x+2*halo_size][mat_1_len+2*halo_size][mat_2_len+2*halo_size]; 
double local_mat_old_new[local_size_x+2*halo_size][mat_1_len+2*halo_size][mat_2_len+2*halo_size]; ;  //　watch out!both new and old are with halo!
double local_H[local_size_x+2*halo_size][mat_1_len+2*halo_size][mat_2_len+2*halo_size];
double local_xgrad[local_size_x][mat_1_len][mat_2_len];
double local_ygrad[local_size_x][mat_1_len][mat_2_len];
double local_zgrad[local_size_x][mat_1_len][mat_2_len];
int local_size_x=mat_0_len/numprocs; //local_size without halo
int local_size_y=mat_1_len;
int local_size_z=mat_2_len; 
int vert[MAX_VERT];
int para[MAX_PARA];
int vert_len ;
int para_len ;

int set_obstacle(int x,int y,int z,int target){ //the input is mat index,transform to local_mat index width halo and change it
	int id=x/local_size_x;
	if(rank==id)
		local_mat_old[x%local_size_x+halo_size][y+halo_size][z+halo_size] = target;
}
double min(double x1, double x2) {  
	if (x1 > x2) return x2;
	else return x1;
}

double getGrad(double ***GV, double x, double y, double z)
bool crossObstacle(double x1, double y1, double z1, double x2, double y2, double z2, int ***local_mat_old, int c)
double findPath(int x1, int y1, int z1, int x2, int y2, int z2, int *vert, int cons) {
	if(rank==0){
		int change = 0;
		int lastpop, nowpop = 0;
        double dist = 0; 
	}
	int icount = 0;          
	double rat = 0.1;
	double speed = 0.5;    
	set_obstacle(x1,y1,z1,1);
	double* p1= &local_mat_old[0][0][0];
	double* p2= &local_mat_new[0][0][0];
	double* buffer[2]={p1,p2};
	while(change<20){
		icount=icount+1;
		double* a0 = buffer[(icount-1) % 2]; 
		double* a1 = buffer[icount%2]; 
		for(int i=0;i<mat_1_len;i++){ 
				MPI_Request request1,request2;
				MPI_Status ierr1,ierr2;
				if(rank<15){  
						MPI_Irecv(&(a0[local_size_x+halo_size][i+halo_size][halo_size]),
						local_size_z,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&request1);
				}   
				if(rank>0){
						MPI_Irecv(&(a0[0][i+halo_size][halo_size]),local_size_z,
						MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&request2);
				}   
				if(rank<15){                                  
						MPI_Send(&(a0[local_size_x][i+halo_size][halo_size]),
						local_size_z,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
				}   
				if(rank>0){
						MPI_Send(&(a0[halo_size][i+halo_size][halo_size]),
						local_size_z,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);
				}    
				if(rank<16-1){
						MPI_Wait(&request1,&ierr1);    
				}   
				if(rank>0){
						MPI_Wait(&request2,&ierr2);    
				}   
			}    
			MPI_Barrier(MPI_COMM_WORLD);
        }   
		for (int i = halo_size; i < local_size_x+halo_size; i++)
		{
			for (int j = halo_size; j < local_size_y+halo_size; j++)
			{
				for (int k = halo_size; k < local_size_z+halo_size; k++)
				{
					a1[i][j][k] \
					= a0[i][j][k] * (1 + rat) + (a0[i - 1][j][k] + a0[i + 1][j][k] +\
					 a0[i][j - 1][k] +a0[i][j + 1][k] + a0[i][j][k - 1] + a0[i][j][k + 1]) * rat;
					a1[i][j][k]=a1[i][j][k] * a0[i][j][k]; //mutiply the mat with obstacle
					if (a1[i][j][k] > 1 && local_H[i][j][k] == 0) local_H[i][j][k] = icount;
				}
			}
		}
		int local_nowpop=0;
		for (int i = halo_size; i < local_size_x+halo_size; i++)
		{
			for (int j = halo_size; j < mat_1_len+halo_size; j++)
			{
				for (int k = halo_size; k < mat_2_len+halo_size; k++)
				{
					if (a1[i][j][k] > 1) local_nowpop = local_nowpop + 1;
				}
			}
		}            
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&local_nowpop,nowpop,1,int,MPI_SUM,0,MPI_COMM_WORLD);
		if(rank==0){
			if (nowpop - lastpop < 1) change = change + 1;
			last_pop=now_pop;
		}
	}

	int id=x2/local_size_x;
	if(rank==id){
		if(local_mat_old[x2%local_size_x+halo_size][y2+halo_size][z2+halo_size] == 0)
			printf("cant spread to the end point \n");
	}

	if(rank==0){
			double ***xgrad = NULL;
			xgrad = (double***)malloc(sizeof(double)*mat_0_len);
			for (int i = 0; i < mat_0_len; i++)
			{
				xgrad[i] = (double**)malloc(sizeof(double)*mat_1_len);
				for (int k = 0; k < mat_1_len; k++)
				{
					xgrad[i][k] = (double*)malloc(sizeof(double)*mat_2_len);
				}
			}
			double ***ygrad = NULL;
			ygrad = (double***)malloc(sizeof(double)*mat_0_len);
			for (int i = 0; i < mat_0_len; i++)
			{
				ygrad[i] = (double**)malloc(sizeof(double)*mat_1_len);
				for (int k = 0; k < mat_1_len; k++)
				{
					ygrad[i][k] = (double*)malloc(sizeof(double)*mat_2_len);
				}
			}
			double ***zgrad = NULL;
			zgrad = (double***)malloc(sizeof(double)*mat_0_len);
			for (int i = 0; i < mat_0_len; i++)
			{
				zgrad[i] = (double**)malloc(sizeof(double)*mat_1_len);
				for (int k = 0; k < mat_1_len; k++)
				{
					zgrad[i][k] = (double*)malloc(sizeof(double)*mat_2_len);
				}
			}
		}
	}

	//update the halo between local_H
	for(int i=0;i<mat_1_len;i++){ 
		MPI_Request request1,request2;
		MPI_Status ierr1,ierr2;
		if(rank<15){  
				MPI_Irecv(&(H[local_size_x+halo_size][i+halo_size][halo_size]),
				local_size_z,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&request1);
		}   
		if(rank>0){
				MPI_Irecv(&(H[0][i+halo_size][halo_size]),local_size_z,
				MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&request2);
		}   
		if(rank<15){                                  
				MPI_Send(&(H[local_size_x][i+halo_size][halo_size]),
				local_size_z,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
		}   
		if(rank>0){
				//ÏòÏÂ·¢
				MPI_Send(&(H[halo_size][i+halo_size][halo_size]),
				local_size_z,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);
		}    
		if(rank<16-1){
				MPI_Wait(&request1,&ierr1);    
		}   
		if(rank>0){
				MPI_Wait(&request2,&ierr2);    
		}   
	}    
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < local_size_x; i++)
	{
		for (int j = 0; j < mat_1_len; j++)
		{
			for (int k = 0; k < mat_2_len; k++)
			{
				local_x_grad[i][j][k] = -((local_H[i + 1+halo_size][j+halo_size][k+halo_size] - local_H[i - 1+halo_size][j+halo_size][k+halo_size]) / 2);       
				local_y_grad[i][j][k] = -((local_H[i+halo_size][j + 1+halo_size][k+halo_size] - local_H[i+halo_size][j - 1+halo_size][k+halo_size]) / 2);     
				local_z_grad[i][j][k] = -((local_H[i+halo_size][j+halo_size][k + 1+halo_size] - local_H[i+halo_size][j+halo_size][k - 1+halo_size]) / 2);    
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&local_xgrad[0][0][0],local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,&xgrad[0][0][0],mat_0_len*mat_1_len*mat_2_len,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&local_ygrad[0][0][0],local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,&ygrad[0][0][0],mat_0_len*mat_1_len*mat_2_len,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&local_ygrad[0][0][0],local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,&ygrad[0][0][0],mat_0_len*mat_1_len*mat_2_len,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if(rank==0){
		double *r = (double*)malloc(sizeof(double) * mat_0_len * mat_1_len * mat_2_len);
		double *c = (double*)malloc(sizeof(double) * mat_0_len * mat_1_len * mat_2_len);
		double *h = (double*)malloc(sizeof(double) * mat_0_len * mat_1_len * mat_2_len);
		double rmax = mat_0_len;
		double cmax = mat_1_len;
		double hmax = mat_2_len;

		r[0] = x2;
		c[0] = y2;
		h[0] = z2;
		int iter = 0;
		bool continflag = true;
		while (continflag)
		{
			double dr, dc, dh, dv;
			dr = getGrad(xgrad, r[iter], c[iter], h[iter]);
			dc = getGrad(ygrad, r[iter], c[iter], h[iter]);
			dh = getGrad(zgrad, r[iter], c[iter], h[iter]);
			dv = sqrt(dr * dr + dc * dc + dh * dh);
			dc = speed * dc / dv;
			dh = speed * dh / dv;
			r[iter + 1] = r[iter] + dr;
			c[iter + 1] = c[iter] + dc;
			h[iter + 1] = h[iter] + dh;
			if (r[iter + 1] < 0)
			{
				r[iter + 1] = 0;
			}
			if (c[iter + 1] < 0)
			{
				c[iter + 1] = 0;
			}
			if (h[iter + 1] < 0)
			{
				h[iter + 1] = 0;
			}
			if (r[iter + 1] > rmax - 1)
			{
				r[iter + 1] = rmax - 1;
			}
			if (c[iter + 1] > cmax - 1)
			{
				c[iter + 1] = cmax - 1;
			}
			if (h[iter + 1] > hmax - 1)
			{
				h[iter + 1] = hmax - 1;
			}                               
			if (sqrt(pow((r[iter + 1] - x1), 2) + pow((c[iter + 1] - y1), 2) + pow((h[iter + 1] - z1), 2)) < 1)
			{
				r[iter + 1] = x1;
				c[iter + 1] = y1;
				h[iter + 1] = z1;
				continflag = false;
			}            
			iter = iter + 1;
		}              

		int infr_Length = 20;
		double infr[20] = { 0 };
		double infc[20] = { 0 };
		double infh[20] = { 0 };       
		infr[0] = r[iter];
		infc[0] = c[iter];
		infh[0] = h[iter];
		int g = 0;                                             
		bool con = true;
		bool ju = true;
		while (con)
		{
			ju = true;
			for (int i = 0; i <= iter; i++)
			{
				if (crossObstacle(r[i], c[i], h[i], r[iter], c[iter], h[iter], local_mat_old, cons) && ju)
				{
					ju = false;
					++g;
					infr[g] = r[i];
					infc[g] = c[i];
					infh[g] = h[i];
					iter = i;
					if (i == 0)
						con = false;
				}

			}
		}
		double s[MAX_VERT] = { 0 };           
		for (int i = 1; i < infr_Length - 1; i++)
		{
			if (infr[i + 1] != 0 || infc[i + 1] != 0 || infh[i + 1] != 0)
			{
				double testx, testy, testz, minest;
				if (infh[i + 1] > infh[i - 1])
					testz = fabs(sqrt(pow((infr[i] - infr[i - 1]), 2) + pow((infc[i] - infc[i - 1]), 2)) / (sqrt(pow((infr[i] - infr[i - 1]), 2) + pow((infc[i] - infc[i - 1]), 2)) + sqrt(pow((infr[i + 1] - infr[i]), 2) + pow((infc[i + 1] - infc[i]), 2))) * fabs(infh[i + 1] - infh[i - 1]) + infh[i - 1] - infh[i]);
				else
					testz = fabs(infh[i - 1] - sqrt(pow((infr[i] - infr[i - 1]), 2) + pow((infc[i] - infc[i - 1]), 2)) / (sqrt(pow((infr[i] - infr[i - 1]), 2) + pow((infc[i] - infc[i - 1]), 2)) + sqrt(pow((infr[i + 1] - infr[i]), 2) + pow((infc[i + 1] - infc[i]), 2))) * fabs(infh[i + 1] - infh[i - 1]) - infh[i]);

				testx = fabs((infc[i] - infc[i + 1]) / (infc[i - 1] - infc[i + 1]) * (infr[i - 1] - infr[i + 1]) + infr[i + 1] - infr[i]);
				testy = fabs((infr[i] - infr[i + 1]) / (infr[i - 1] - infr[i + 1]) * (infc[i - 1] - infc[i + 1]) + infc[i + 1] - infc[i]);
				minest = min(min(testx, testy), testz);
				if (testz == minest)
				{
					for (int j = 0; j < MAX_VERT - 5; j++)
					{
						if (j % 6 == 0)
						{
							s[j] = sqrt(pow((infr[i] - vert[j]), 2) + pow((infc[i] - vert[j + 1]), 2));
							s[j + 2] = sqrt(pow((infr[i] - vert[j + 2]), 2) + pow((infc[i] - vert[j + 3]), 2));
						}
					}
					double mins = s[0];
					int index = 0;
					for (int k = 0; k < MAX_VERT; k++)
					{
						if (s[k] != 0)
						{
							if (s[k] == min(s[k], mins))
							{
								mins = s[k];
								index = k;
							}
						}
					}
					infr[i] = vert[index];
					infc[i] = vert[index + 1];
				}
				else
				{

				}
			}
		}                       
		double dist1[20] = { 0 };
		int dist1_Length = 20;
		for (int i = 0; i < infr_Length - 1; i++)
		{
			if (infr[i + 1] != 0 || infc[i + 1] != 0 || infh[i + 1] != 0)
				dist1[i] = sqrt(pow((infr[i] - infr[i + 1]), 2) + pow((infc[i] - infc[i + 1]), 2) + pow((infh[i] - infh[i + 1]), 2));
		}
		for (int i = 0; i < dist1_Length; i++)
			printf("dist1[%d]=%f\n", i, dist1[i]);
		for (int i = 0; i < infr_Length; i++)
			dist = dist + dist1[i];
		return dist;
	}
}

int main(int argc,char **argv)
{
	int rank,numprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	if(rank==0){
		char str[MAX_LEN] = { 0 };
		FILE* sr;

		sr = fopen("vert.txt", "r");
		if (sr == NULL) {
			puts("can not open file!");
			return 0;
		}
		fgets(str, MAX_LEN, sr);
		while (!feof(sr)) {
			fgets(str, MAX_LEN, sr);
			char * token = strtok(str, " ");
			while (token != NULL) {
				vert[vert_len++] = atoi(token);
				token = strtok(NULL, " ");
			}
		}
		fclose(sr);

		sr = fopen("para.txt", "r");
		if (sr == NULL) {
			puts("can not open file!");
			return 0;
		}
		fgets(str, MAX_LEN, sr);
		while (!feof(sr)) {
			fgets(str, MAX_LEN, sr);
			char * token = strtok(str, " ");
			while (token != NULL) {
				para[para_len++] = atoi(token);
				token = strtok(NULL, " ");
			}
		}
		fclose(sr);


		double end,start;
		start=MPI_Wtime();         
	
		athread_init();
		athread_set_num_threads(nums);
		/*
		__real_athread_spawn((void *)func,0);	
		*/
		MPI_Bcast(verta,vert_len,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&vert_len,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(para,vert_len,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&para_len,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	for (int i = 0; i < vert_len - 5; i++)
	{
		if (i % 6 == 0)
		{
			int h;
			h = vert[i + 5] - vert[i + 4] + 1;
			if (vert[i + 1] == vert[i + 3])
			{
				if (vert[i] < vert[i + 2])
				{
					for (int j = vert[i + 4]; j < h; j++)
					{
						for (int k = vert[i]; k <= vert[i + 2]; k++)
						{
							set_obstacle(k,vert[i+1],j,1);
						}
					}
				}
				else
				{
					for (int j = vert[i + 4]; j < h; j++)
					{
						for (int k = vert[i + 2]; k <= vert[i]; k++)
						{
						    set_obstacle(k,vert[i+1],j,1);
						}
					}
				}
			}
			else if (vert[i] == vert[i + 2])
			{
				if (vert[i + 1] < vert[i + 3])
				{
					for (int j = vert[i + 4]; j < h; j++)
					{
						for (int k = vert[i + 1]; k <= vert[i + 3]; k++)
						{
							set_obstacle(vert[i],k,j,1);
						}
					}
				}
				else
				{
					for (int j = vert[i + 4]; j < h; j++)
					{
						for (int k = vert[i + 3]; k <= vert[i + 1]; k++)
						{
							set_obstacle(vert[i],k,j,1);
						}
					}
				}
			}
			else
			{
				double a = (double)(vert[i + 1] - vert[i + 3]) / (vert[i] - vert[i + 2]);
				double b = (double)(vert[i + 1] - a * vert[i]);
				if (fabs(vert[i + 1] - vert[i + 3]) <= fabs(vert[i] - vert[i + 2]))
				{
					if (vert[i] < vert[i + 2])
					{
						for (int j = vert[i + 4]; j < h; j++)
						{
							for (int k = vert[i]; k <= vert[i + 2]; k++)
							{
								set_obstacle(k,(int)(a * k + b + 0.5),j,1)
								set_obstacle(k+1,(int)(a*k+b+0.5),j,1);
							}
						}
					}
					else
					{
						for (int j = vert[i + 4]; j < h; j++)
						{
							for (int k = vert[i + 2]; k <= vert[i]; k++)
							{
								set_obstacle(k,(int)(a * k + b + 0.5),j,1)
								set_obstacle(k+1,(int)(a*k+b+0.5),j,1);
							}
						}
					}
				}
				else
				{
					if (vert[i + 1] < vert[i + 3])
					{
						for (int j = vert[i + 4]; j < h; j++)
						{
							for (int k = vert[i + 1]; k <= vert[i + 3]; k++)
							{
								set_obstacle((int)((k - b) / a + 0.5),k,j,1);
								set_obstacle((int)((k - b) / a + 0.5),k+1,j,1);
							}
						}
					}
					else
					{
						for (int j = vert[i + 4]; j < h; j++)
						{
							for (int k = vert[i + 3]; k <= vert[i + 1]; k++)
							{
								set_obstacle((int)((k - b) / a + 0.5),k,j,1);
								set_obstacle((int)((k - b) / a + 0.5),k+1,j,1);
							}
						}
					}
				}

			}
		}
	}

	for (int i = 0; i < para_len - 5; i++)
	{
		if (para[i] == para[i + 2])
		{
			if (para[i + 1] < para[i + 3])
			{
				for (int k = para[i] - para[i + 5] / 2; k <= para[i] + para[i + 5] / 2; k++)
				{
					for (int j = para[i + 1]; j <= para[i + 3]; j++)
					{
						set_obstacle(k,j,para[i + 4],1);
					}
				}
			}
			else
			{
				for (int k = para[i] - para[i + 5] / 2; k <= para[i] + para[i + 5] / 2; k++)
				{
					for (int j = para[i + 3]; j <= para[i + 1]; j++)
					{
						set_obstacle(k,j,para[i + 4],1);
					}
				}

			}

		}
	}
	for (int i = 0; i < mat_0_len; i++)
	{
		for (int j = 0; j < mat_1_len; j++)
		{
			for (int k = 0; k < mat_2_len; k++)
			{
				if (local_mat_old[i][j][k] == 1)
				{
					set_obstacle(i,j,k,0);
					c++;
				}
				else
				{
					set_obstacle(i,j,k,1);
				}
			}
		}
	}

	double dist = findPath(10, 5, 5, 15, 80, 5, vert, c);
	if(myid==0){
		printf("dist = %lf", dist);
		end=MPI_Wtime();
		printf("\ttime:%lf\n");
	}

	athread_halt();
	MPI_Finalize();
	return 0;
}













double getGrad(double ***GV, double x, double y, double z)
{
	double value;
	double x1, x2, y1, y2, z1, z2;
	if ((int)(x) != x && (int)(y) != y && (int)(z) != z)
	{
		x1 = (ceil(x) - x) * GV[(int)(x)][(int)(y)][(int)(z)] + (x - (int)(x)) * GV[(int)ceil(x)][(int)(y)][(int)(z)];
		x2 = (ceil(x) - x) * GV[(int)(x)][(int)ceil(y)][(int)(z)] + (x - (int)(x)) * GV[(int)ceil(x)][(int)ceil(y)][(int)(z)];
		z1 = (ceil(y) - y) * x1 + (y - (int)(y)) * x2;
		y1 = (ceil(x) - x) * GV[(int)(x)][(int)(y)][(int)ceil(z)] + (x - (int)(x)) * GV[(int)ceil(x)][(int)(y)][(int)ceil(z)];
		y2 = (ceil(x) - x) * GV[(int)(x)][(int)ceil(y)][(int)ceil(z)] + (x - (int)(x)) * GV[(int)ceil(x)][(int)ceil(y)][(int)ceil(z)];
		z2 = (ceil(y) - y) * y1 + (y - (int)(y)) * y2;
		value = (ceil(z) - z) * z1 + (z - (int)(z)) * z2;
	}
	else if ((int)(x) == x && (int)(y) != y && (int)(z) != z)
	{
		x1 = GV[(int)x][(int)(y)][(int)(z)];
		x2 = GV[(int)x][(int)ceil(y)][(int)(z)];
		z1 = (ceil(y) - y) * x1 + (y - (int)(y)) * x2;
		y1 = GV[(int)x][(int)(y)][(int)ceil(z)];
		y2 = GV[(int)x][(int)ceil(y)][(int)ceil(z)];
		z2 = (ceil(y) - y) * y1 + (y - (int)(y)) * y2;
		value = (ceil(z) - z) * z1 + (z - (int)(z)) * z2;
	}
	else if ((int)(x) != x && (int)(y) == y && (int)(z) != z)
	{
		x1 = GV[(int)(x)][(int)y][(int)(z)];
		x2 = GV[(int)ceil(x)][(int)y][(int)(z)];
		z1 = (ceil(x) - x) * x1 + (x - (int)(x)) * x2;
		y1 = GV[(int)(x)][(int)y][(int)ceil(z)];
		y2 = GV[(int)ceil(x)][(int)y][(int)ceil(z)];
		z2 = (ceil(y) - y) * y1 + (y - (int)(y)) * y2;
		value = (ceil(z) - z) * z1 + (z - (int)(z)) * z2;

	}
	else if ((int)(x) != x && (int)(y) != y && (int)(z) == z)
	{
		x1 = (ceil(x) - x) * GV[(int)(x)][(int)(y)][(int)z] + (x - (int)(x)) * GV[(int)ceil(x)][(int)(y)][(int)z];
		x2 = (ceil(x) - x) * GV[(int)(x)][(int)ceil(y)][(int)z] + (x - (int)(x)) * GV[(int)ceil(x)][(int)ceil(y)][(int)z];
		value = (ceil(y) - y) * x1 + (y - (int)(y)) * x2;
	}
	else if ((int)(x) != x && (int)(y) == y && (int)(z) == z)
	{
		value = (ceil(x) - x) * GV[(int)(x)][(int)y][(int)z] + (x - (int)(x)) * GV[(int)ceil(x)][(int)y][(int)z];
	}
	else if ((int)(x) == x && (int)(y) != y && (int)(z) == z)
	{
		value = (ceil(y) - y) * GV[(int)x][(int)(y)][(int)z] + (y - (int)(y)) * GV[(int)x][(int)ceil(y)][(int)z];
	}
	else if ((int)(x) == x && (int)(y) == y && (int)(z) != z)
	{
		value = (ceil(z) - z) * GV[(int)x][(int)y][(int)(z)] + (z - (int)(z)) * GV[(int)x][(int)y][(int)ceil(z)];
	}
	else
	{
		value = GV[(int)x][(int)y][(int)z];
	}

	return value;
}

bool crossObstacle(double x1, double y1, double z1, double x2, double y2, double z2, int ***local_mat_old, int c)
{

	if (x1 < x2)
	{
		for (int x = (int)floor(x1); x < (int)ceil(x2); x++)
		{
			int y = (int)round((x - x1) / (x2 - x1) * (y2 - y1) + y1);
			int z = (int)round((x - x1) / (x2 - x1) * (z2 - z1) + z1);
			if (y > 0 && y < mat_1_len && z > 0 && z < mat_2_len&&local_mat_old[x][y][z]==0) return false;
		}             
	}
	else if (x1 == x2);
	else
	{
		for (int x = (int)floor(x2); x < (int)ceil(x1); x++)
		{
			int y = (int)round((x - x1) / (x2 - x1) * (y2 - y1) + y1);
			int z = (int)round((x - x1) / (x2 - x1) * (z2 - z1) + z1);
			if (y > 0 && y < mat_1_len && z > 0 && z < mat_2_len&&local_mat_old[x][y][z]==0) return false;
		}               
	}
	if (y1 < y2)
	{
		for (int y = (int)floor(y1); y < (int)ceil(y2); y++)
		{
			int x = (int)round((y - y1) / (y2 - y1) * (x2 - x1) + x1);
			int z = (int)round((y - y1) / (y2 - y1) * (z2 - z1) + z1);
			if (y > 0 && y < mat_1_len && z > 0 && z < mat_2_len&&local_mat_old[x][y][z]==0) return false;
		}                  
	}
	else if (y1 == y2);
	else
	{
		for (int y = (int)floor(y2); y < (int)ceil(y1); y++)
		{
			int x = (int)round((y - y1) / (y2 - y1) * (x2 - x1) + x1);
			int z = (int)round((y - y1) / (y2 - y1) * (z2 - z1) + z1);
			if (y > 0 && y < mat_1_len && z > 0 && z < mat_2_len&&local_mat_old[x][y][z]==0) return false;
		}                 
	}
	if (z1 < z2)
	{
		for (int z = (int)floor(z1); z < (int)ceil(z2); z++)
		{
			int x = (int)round((z - z1) / (z2 - z1) * (x2 - x1) + x1);
			int y = (int)round((z - z1) / (z2 - z1) * (y2 - y1) + y1);
			if (y > 0 && y < mat_1_len && z > 0 && z < mat_2_len&&local_mat_old[x][y][z]==0) return false;
		}                 
	}
	else if (z1 == z2);
	else
	{
		for (int z = (int)floor(z2); z < (int)ceil(z1); z++)
		{
			int x = (int)round((z - z1) / (z2 - z1) * (x2 - x1) + x1);
			int y = (int)round((z - z1) / (z2 - z1) * (y2 - y1) + y1);
			if (y > 0 && y < mat_1_len && z > 0 && z < mat_2_len&&local_mat_old[x][y][z]==0) return false;
		}               
	}
	return true;
}