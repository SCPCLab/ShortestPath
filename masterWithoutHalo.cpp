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
    void slave_func(void);
}

#define MAX_LEN 30
#define MAX_VERT 1000
#define MAX_PARA 1000
#define mat_0_len 100
#define mat_1_len 100
#define mat_2_len 100
#define halo_size 1
#define local_size_x 10
#define local_size_y 100
#define local_size_z 100
int rank,numprocs;


#define INDEX(x,y,z,ldy,ldz)  ( ((y)+(x)*(ldy))*(ldz) +(z) )

int *local_mat_3d=NULL;
double *local_gradmat=NULL;
double *local_gradmatnew=NULL;
double *local_H=NULL;
double *local_xgrad=NULL;
double *local_ygrad=NULL;
double *local_zgrad=NULL;

int *mat3d=NULL;
double *test=NULL;

int vert[MAX_VERT];
int para[MAX_PARA];	
void set_obstacle(int x,int y,int z,int target){ //输入是全局矩阵坐标，需要转换为带环晕的局部矩阵并设置障碍物
	int id=x/local_size_x;
	if(rank==id){ 
		local_mat_3d[INDEX(x%local_size_x,y,z,local_size_y,local_size_z)] = target;
	}
}
double min(double x1, double x2) {  
	if (x1 > x2) return x2;
	else return x1;
}

double getGrad(double *GV, double x, double y, double z);
bool crossObstacle(double x1, double y1, double z1, double x2, double y2, double z2, int *mat_3d, int c, double *test);

double findPath(int x1, int y1, int z1, int x2, int y2, int z2, int *vert, int cons) {
	int change = 0;
	int lastpop=0, nowpop = 0;
       	double dist = 0; 
	int icount = 0;          
	double rat = 0.1;
	double speed = 0.5;    
	int ldy=local_size_y;
	int ldz=local_size_z;
	if(x1/local_size_x == rank){
		local_gradmat[INDEX(x1%local_size_x,y1,z1,ldy,ldz)] = 1;
	}
	
	double *upper=NULL;
	double *down=NULL;
	if(rank<numprocs-1){  
		upper=(double*)malloc(sizeof(double)*ldy*ldz);
	}
	if(rank>0){
		down=(double*)malloc(sizeof(double)*ldy*ldz);
	}
	

	
	
	while(change<20){
		icount=icount+1; 
		for(int i=0;i<mat_1_len;i++){ 
				MPI_Request request1,request2;
				MPI_Status ierr1,ierr2;
				if(rank<numprocs-1){  
					MPI_Irecv(&(upper[i*ldz]),local_size_z,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&request1);
					
				}   
				
				if(rank>0){	
					MPI_Irecv(&(down[i*ldz]),local_size_z,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&request2);	
				}
			  
				if(rank<numprocs-1){                                  
					MPI_Send(&(local_gradmat[INDEX(local_size_x-1,i,0,ldy,ldz)]),local_size_z,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
					
				}   
				if(rank>0){		
					MPI_Send(&(local_gradmat[INDEX(0,i,0,ldy,ldz)]),local_size_z,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);
					
				}  
			
				if(rank<numprocs-1){
					MPI_Wait(&request1,&ierr1);    
				}   
				if(rank>0){
					MPI_Wait(&request2,&ierr2);    
				}   
		}    
		MPI_Barrier(MPI_COMM_WORLD);
		

		for (int i = 0; i < local_size_x; i++)
		{
			if(i==0&&rank==0||i==local_size_x-1&&rank==numprocs-1) continue;
			for (int j = 1; j < local_size_y - 1; j++)
			{
				for (int k = 1; k < local_size_z - 1; k++)
				{
					if(i==0){
						local_gradmatnew[INDEX(i,j,k,ldy,ldz)] = local_gradmat[INDEX(i,j,k,ldy,ldz)] * (1 + rat) + (down[j*local_size_z+k] + local_gradmat[INDEX(i + 1,j,k,ldy,ldz)]
						+ local_gradmat[INDEX(i,j - 1,k,ldy,ldz)]+ local_gradmat[INDEX(i,j + 1,k,ldy,ldz)]+ local_gradmat[INDEX(i,j,k - 1,ldy,ldz)]+  local_gradmat[INDEX(i,j,k + 1,ldy,ldz)]) * rat;
					}else if(i==local_size_x-1){
						local_gradmatnew[INDEX(i,j,k,ldy,ldz)] = local_gradmat[INDEX(i,j,k,ldy,ldz)] * (1 + rat) + (local_gradmat[INDEX(i - 1,j,k,ldy,ldz)] + upper[j*local_size_z+k]
						+ local_gradmat[INDEX(i,j - 1,k,ldy,ldz)]+ local_gradmat[INDEX(i,j + 1,k,ldy,ldz)]+ local_gradmat[INDEX(i,j,k - 1,ldy,ldz)]+  local_gradmat[INDEX(i,j,k + 1,ldy,ldz)]) * rat;

					}else{
						local_gradmatnew[INDEX(i,j,k,ldy,ldz)] = local_gradmat[INDEX(i,j,k,ldy,ldz)] * (1 + rat) + (local_gradmat[INDEX(i - 1,j,k,ldy,ldz)] + local_gradmat[INDEX(i + 1,j,k,ldy,ldz)]
						+ local_gradmat[INDEX(i,j - 1,k,ldy,ldz)]+ local_gradmat[INDEX(i,j + 1,k,ldy,ldz)]+ local_gradmat[INDEX(i,j,k - 1,ldy,ldz)]+  local_gradmat[INDEX(i,j,k + 1,ldy,ldz)]) * rat;
					}
				}
			}
		}

		for (int i = 0; i < local_size_x; i++)
		{
			for (int j = 0; j < local_size_y; j++)
			{
				for (int k = 0; k < local_size_z; k++)
				{
					local_gradmat[INDEX(i,j,k,ldy,ldz)]= local_gradmatnew[INDEX(i,j,k,ldy,ldz)] * local_mat_3d[INDEX(i,j,k,ldy,ldz)]; //mutiply the mat with obstacle
				
					if (local_gradmat[INDEX(i,j,k,ldy,ldz)] > 1 && local_H[INDEX(i,j,k,ldy,ldz)]==0){
						local_H[INDEX(i,j,k,ldy,ldz)] = icount;	
					}
				}
			}
		}
		if(rank==0){
			lastpop=nowpop; nowpop=0;
		} 
		int local_nowpop=0;
		for (int i = 0; i < local_size_x; i++)
		{
			for (int j = 0; j < mat_1_len; j++)
			{
				for (int k = 0; k < mat_2_len; k++)
				{
					if (local_gradmat[INDEX(i,j,k,ldy,ldz)]>1) local_nowpop = local_nowpop + 1;
				}
			}
		}            
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Reduce(&local_nowpop,&nowpop,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		
		if(rank==0){
			if (nowpop - lastpop < 1) change = change + 1;
		}
		MPI_Bcast(&change,1,MPI_INT,0,MPI_COMM_WORLD);

	} 
	int id=x2/local_size_x;
	if(rank==id){
		if(local_gradmat[INDEX(x2%local_size_x,y2,z2,ldy,ldz)] == 0)
			printf("cant spread to the end point \n");
	}
	double *xgrad = NULL;
	double  *ygrad = NULL;
	double *zgrad = NULL;
	
	for(int i=0;i<local_size_x;i++){
		for(int j=0;j<local_size_y;j++){
			for(int k=0;k<local_size_z;k++){
				if(local_H[INDEX(i,j,k,ldy,ldz)]==0){
					local_H[INDEX(i,j,k,ldy,ldz)]=icount;
				}
			}
		}
	}
	//更新local_H之间的环晕部分
	for(int i=0;i<mat_1_len;i++){ 
		MPI_Request request1,request2;
		MPI_Status ierr1,ierr2;
		if(rank<numprocs-1){  
			MPI_Irecv(&(upper[i*local_size_z]),local_size_z,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&request1);
		}   
		if(rank>0){
			MPI_Irecv(&(down[i*local_size_z]),local_size_z,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&request2);	
		}   
		if(rank<numprocs-1){                                  
			MPI_Send(&(local_H[INDEX(local_size_x-1,i,0,ldy,ldz)]),local_size_z,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);		
		}   
		if(rank>0){
			MPI_Send(&(local_H[INDEX(0,i,0,ldy,ldz)]),local_size_z,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);
		}    
		if(rank<numprocs-1){
			MPI_Wait(&request1,&ierr1);    
		}   
		if(rank>0){
			MPI_Wait(&request2,&ierr2);    
		}   
	}    
	MPI_Barrier(MPI_COMM_WORLD);

	free(local_gradmatnew);
	free(local_gradmat);
	free(local_mat_3d);
 

	ldy=local_size_y;ldz=local_size_z;
	local_xgrad=(double*)malloc(sizeof(double)*local_size_x*local_size_y*local_size_z);
	local_ygrad=(double*)malloc(sizeof(double)*local_size_x*local_size_y*local_size_z);
	local_zgrad=(double*)malloc(sizeof(double)*local_size_x*local_size_y*local_size_z);
	if(local_xgrad==NULL||local_ygrad==NULL||local_zgrad==NULL) printf("malloc  failed\n");

	memset(local_xgrad,0,sizeof(double)*local_size_x*local_size_y*local_size_z);
	memset(local_ygrad,0,sizeof(double)*local_size_x*local_size_y*local_size_z);
	memset(local_zgrad,0,sizeof(double)*local_size_x*local_size_y*local_size_z);

	for (int i = 0; i < local_size_x; i++)
	{	
		for (int j = 0; j < local_size_y; j++)
		{
			for (int k = 0; k < local_size_z; k++)
			{
				if(rank==0&&i==0){
					local_xgrad[INDEX(i,j,k,ldy,ldz)] = -(local_H[INDEX(i + 1,j,k,ldy,ldz)] - local_H[INDEX(i,j,k,ldy,ldz)]) ;					       
				}else if((rank==numprocs-1)&&(i==local_size_x-1)){
					local_xgrad[INDEX(i,j,k,ldy,ldz)] = -(local_H[INDEX(i,j,k,ldy,ldz)] - local_H[INDEX(i - 1,j,k,ldy,ldz)]) ;	
				}else{
					if(i==0){
						local_xgrad[INDEX(i,j,k,ldy,ldz)] = -((local_H[INDEX(i + 1,j,k,ldy,ldz)] - down[j*local_size_z+k]) / 2);
					}else if(i==local_size_x-1){
						local_xgrad[INDEX(i,j,k,ldy,ldz)] = -((upper[j*local_size_z+k] - local_H[INDEX(i - 1,j,k,ldy,ldz)]) / 2);
					}else{
			
						local_xgrad[INDEX(i,j,k,ldy,ldz)] = -((local_H[INDEX(i + 1,j,k,ldy,ldz)] - local_H[INDEX(i - 1,j,k,ldy,ldz)]) / 2);
					}					       
				}   
				if(j==0){
					local_ygrad[INDEX(i,j,k,ldy,ldz)] = -(local_H[INDEX(i,j + 1,k,ldy,ldz)] - local_H[INDEX(i,j,k,ldy,ldz)]);	
				}else if(j==local_size_y-1){
					local_ygrad[INDEX(i,j,k,ldy,ldz)] = -(local_H[INDEX(i,j,k,ldy,ldz)] - local_H[INDEX(i,j - 1,k,ldy,ldz)]);
				}else{
					local_ygrad[INDEX(i,j,k,ldy,ldz)] = -((local_H[INDEX(i,j + 1,k,ldy,ldz)] - local_H[INDEX(i,j - 1,k,ldy,ldz)]) / 2);
				}

				if(k==0){
					local_zgrad[INDEX(i,j,k,ldy,ldz)] = -(local_H[INDEX(i,j,k + 1,ldy,ldz)] - local_H[INDEX(i,j,k,ldy,ldz)]);	
				}else if(k==local_size_z-1){
					local_zgrad[INDEX(i,j,k,ldy,ldz)] = -(local_H[INDEX(i,j,k,ldy,ldz)] - local_H[INDEX(i,j,k - 1,ldy,ldz)]);
				}else{
					local_zgrad[INDEX(i,j,k,ldy,ldz)] = -((local_H[INDEX(i,j,k + 1,ldy,ldz)] - local_H[INDEX(i,j,k - 1,ldy,ldz)]) / 2);
				}  
 
			}
		}
	}

	if(rank==0){
			xgrad = (double*)malloc(sizeof(double)*mat_0_len*mat_1_len*mat_2_len);	
			ygrad = (double*)malloc(sizeof(double)*mat_0_len*mat_1_len*mat_2_len);
			zgrad = (double*)malloc(sizeof(double)*mat_0_len*mat_1_len*mat_2_len);
			if(xgrad==NULL||ygrad==NULL||zgrad==NULL) printf("gradx y z  malloc  failed\n");			
	}
	
	MPI_Barrier(MPI_COMM_WORLD); //等待主进程分配内存
	
	MPI_Gather(local_xgrad,local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,xgrad,local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(local_ygrad,local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,ygrad,local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(local_zgrad,local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,zgrad,local_size_x*mat_1_len*mat_2_len,MPI_DOUBLE,0,MPI_COMM_WORLD);
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
		srand((time(0)));
		while (continflag)
		{
			double dr, dc, dh, dv;
			dr = getGrad(xgrad, r[iter], c[iter], h[iter]);
			dc = getGrad(ygrad, r[iter], c[iter], h[iter]);
			dh = getGrad(zgrad, r[iter], c[iter], h[iter]);
			
			dv = sqrt(dr * dr + dc * dc + dh * dh);

			dr = speed * dr / dv;
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
				if (crossObstacle(r[i], c[i], h[i], r[iter], c[iter], h[iter], mat3d, cons, test) && ju)
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
	double start,end;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

	local_mat_3d=(int*)malloc(sizeof(int)*local_size_x*local_size_y*local_size_z);

	local_gradmat=(double*)malloc(sizeof(double)*local_size_x*local_size_y*local_size_z);
	local_gradmatnew=(double*)malloc(sizeof(double)*local_size_x*local_size_y*local_size_z);
	local_H=(double*)malloc(sizeof(double)*local_size_x*local_size_y*local_size_z);
	
	test=(double*)malloc(sizeof(double)*mat_0_len*mat_1_len*mat_2_len);
	mat3d=(int*)malloc(sizeof(int)*mat_0_len*mat_1_len*mat_2_len);

	memset(local_mat_3d,0,sizeof(int)*local_size_x*local_size_y*local_size_z);
	memset(local_gradmat,0,sizeof(double)*local_size_x*local_size_y*local_size_z);
	memset(local_gradmatnew,0,sizeof(double)*local_size_x*local_size_y*local_size_z);
	memset(local_H,0,sizeof(double)*local_size_x*local_size_y*local_size_z);	

	int vert_len=0;
	int para_len=0;
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
	}
	start=MPI_Wtime();
	MPI_Bcast(&para_len,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&vert_len,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Bcast(vert,vert_len,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(para,para_len,MPI_INT,0,MPI_COMM_WORLD);
	
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
								set_obstacle(k,(int)(a * k + b + 0.5),j,1);
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
								set_obstacle(k,(int)(a * k + b + 0.5),j,1);
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
	int c=0;
	int sum=0;	
	for (int i = 0; i < local_size_x; i++)
	{
		for (int j = 0; j < local_size_y; j++)
		{
			for (int k = 0; k < local_size_z; k++)
			{	
				if (local_mat_3d[INDEX(i,j,k,local_size_y,local_size_z)]==1.0)
				{					
					local_mat_3d[INDEX(i,j,k,local_size_y,local_size_z)] = 0;
					c++;
				}
				else
				{
					local_mat_3d[INDEX(i,j,k,local_size_y,local_size_z)] = 1;
				}
			}
		}
	}
	MPI_Gather(local_mat_3d,local_size_x*mat_1_len*mat_2_len,MPI_INT,mat3d,local_size_x*mat_1_len*mat_2_len,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Reduce(&c,&sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	double dist = findPath(10, 5, 5, 15, 80, 5, vert, sum);
	end=MPI_Wtime();
	if(rank==0){
		printf("dist = %lf", dist);
		printf("\ttime:%lf\n",end-start);
	}

	MPI_Finalize();
	return 0;
}













double getGrad(double *GV, double x, double y, double z)
{
	double value;
	double x1, x2, y1, y2, z1, z2;
	if ((int)(x) != x && (int)(y) != y && (int)(z) != z)     							
	{
		x1 = (ceil(x) - x) * GV[INDEX((int)(x),(int)(y),(int)(z),local_size_y,local_size_z)] + (x - (int)(x)) * GV[INDEX((int)ceil(x),(int)(y),(int)(z),local_size_y,local_size_z)];
		x2 = (ceil(x) - x) * GV[INDEX((int)(x),(int)ceil(y),(int)(z),local_size_y,local_size_z)] + (x - (int)(x)) * GV[INDEX((int)ceil(x),(int)ceil(y),(int)(z),local_size_y,local_size_z)];
		z1 = (ceil(y) - y) * x1 + (y - (int)(y)) * x2;
		y1 = (ceil(x) - x) * GV[INDEX((int)(x),(int)(y),(int)ceil(z),local_size_y,local_size_z)] + (x - (int)(x)) * GV[INDEX((int)ceil(x),(int)(y),(int)ceil(z),local_size_y,local_size_z)];
		y2 = (ceil(x) - x) * GV[INDEX((int)(x),(int)ceil(y),(int)ceil(z),local_size_y,local_size_z)] + (x - (int)(x)) * GV[INDEX((int)ceil(x),(int)ceil(y),(int)ceil(z),local_size_y,local_size_z)];
		z2 = (ceil(y) - y) * y1 + (y - (int)(y)) * y2;
		value = (ceil(z) - z) * z1 + (z - (int)(z)) * z2;
	}
	else if ((int)(x) == x && (int)(y) != y && (int)(z) != z)
	{
		x1 = GV[INDEX((int)(x),(int)(y),(int)(z),local_size_y,local_size_z)];
		x2 = GV[INDEX((int)(x),(int)ceil(y),(int)(z),local_size_y,local_size_z)];
		z1 = (ceil(y) - y) * x1 + (y - (int)(y)) * x2;
		y1 = GV[INDEX((int)(x),(int)(y),(int)ceil(z),local_size_y,local_size_z)];
		y2 = GV[INDEX((int)(x),(int)ceil(y),(int)ceil(z),local_size_y,local_size_z)];
		z2 = (ceil(y) - y) * y1 + (y - (int)(y)) * y2;
		value = (ceil(z) - z) * z1 + (z - (int)(z)) * z2;
	}
	else if ((int)(x) != x && (int)(y) == y && (int)(z) != z)
	{
		x1 = GV[INDEX((int)(x),(int)(y),(int)(z),local_size_y,local_size_z)];
		x2 = GV[INDEX((int)ceil(x),(int)(y),(int)(z),local_size_y,local_size_z)];
		z1 = (ceil(x) - x) * x1 + (x - (int)(x)) * x2;
		y1 = GV[INDEX((int)(x),(int)(y),(int)ceil(z),local_size_y,local_size_z)];
		y2 = GV[INDEX((int)ceil(x),(int)(y),(int)ceil(z),local_size_y,local_size_z)];
		z2 = (ceil(y) - y) * y1 + (y - (int)(y)) * y2;
		value = (ceil(z) - z) * z1 + (z - (int)(z)) * z2;

	}
	else if ((int)(x) != x && (int)(y) != y && (int)(z) == z)
	{
		x1 = (ceil(x) - x) * GV[INDEX((int)(x),(int)(y),(int)(z),local_size_y,local_size_z)] + (x - (int)(x)) * GV[INDEX((int)ceil(x),(int)(y),(int)(z),local_size_y,local_size_z)];
		x2 = (ceil(x) - x) * GV[INDEX((int)(x),(int)ceil(y),(int)(z),local_size_y,local_size_z)] + (x - (int)(x)) * GV[INDEX((int)ceil(x),(int)ceil(y),(int)(z),local_size_y,local_size_z)];
		value = (ceil(y) - y) * x1 + (y - (int)(y)) * x2;
	}
	else if ((int)(x) != x && (int)(y) == y && (int)(z) == z)
	{
		value = (ceil(x) - x) * GV[INDEX((int)(x),(int)(y),(int)(z),local_size_y,local_size_z)] + (x - (int)(x)) * GV[INDEX((int)ceil(x),(int)(y),(int)(z),local_size_y,local_size_z)];
	}
	else if ((int)(x) == x && (int)(y) != y && (int)(z) == z)
	{
		value = (ceil(y) - y) * GV[INDEX((int)(x),(int)(y),(int)(z),local_size_y,local_size_z)] + (y - (int)(y)) * GV[INDEX((int)(x),(int)ceil(y),(int)(z),local_size_y,local_size_z)];
	}
	else if ((int)(x) == x && (int)(y) == y && (int)(z) != z)
	{
		value = (ceil(z) - z) * GV[INDEX((int)(x),(int)(y),(int)(z),local_size_y,local_size_z)] + (z - (int)(z)) * GV[INDEX((int)(x),(int)(y),(int)ceil(z),local_size_y,local_size_z)];
	}
	else
	{
		value = GV[INDEX((int)(x),(int)(y),(int)(z),local_size_y,local_size_z)];
	}

	return value;
}

bool crossObstacle(double x1, double y1, double z1, double x2, double y2, double z2, int *mat_3d, int c, double *test)
{
	for (int i = 0; i < mat_0_len; i++)
	{
		for (int j = 0; j < mat_1_len; j++)
		{
			for (int k = 0; k < mat_2_len; k++)
			{
				test[INDEX(i,j,k,local_size_y,local_size_z)] = mat_3d[INDEX(i,j,k,local_size_y,local_size_z)];
			}
		}
	}

	bool cross;
	int c1 = 0;

	if (x1 < x2)
	{
		for (int x = (int)floor(x1); x < (int)ceil(x2); x++)
		{
			int y = (int)round((x - x1) / (x2 - x1) * (y2 - y1) + y1);
			int z = (int)round((x - x1) / (x2 - x1) * (z2 - z1) + z1);
			if (y > 0 && y < mat_1_len && z > 0 && z < mat_2_len)
			{
				test[INDEX(x,y,z,local_size_y,local_size_z)] = 1;
			}

		}//空间中直线方程              
	}
	else if (x1 == x2);
	else
	{
		for (int x = (int)floor(x2); x < (int)ceil(x1); x++)
		{
			int y = (int)round((x - x1) / (x2 - x1) * (y2 - y1) + y1);
			int z = (int)round((x - x1) / (x2 - x1) * (z2 - z1) + z1);
			if (y > 0 && y < mat_1_len && z > 0 && z < mat_2_len)
			{
				test[INDEX(x,y,z,local_size_y,local_size_z)] = 1;
			}

		}//空间中直线方程          
	}
	if (y1 < y2)
	{
		for (int y = (int)floor(y1); y < (int)ceil(y2); y++)
		{
			int x = (int)round((y - y1) / (y2 - y1) * (x2 - x1) + x1);
			int z = (int)round((y - y1) / (y2 - y1) * (z2 - z1) + z1);
			if (x > 0 && x < mat_0_len && z > 0 && z < mat_2_len)
			{
				test[INDEX(x,y,z,local_size_y,local_size_z)] = 1;
			}

		}//空间中直线方程              
	}
	else if (y1 == y2);
	else
	{
		for (int y = (int)floor(y2); y < (int)ceil(y1); y++)
		{
			int x = (int)round((y - y1) / (y2 - y1) * (x2 - x1) + x1);
			int z = (int)round((y - y1) / (y2 - y1) * (z2 - z1) + z1);
			if (x > 0 && x < mat_0_len && z > 0 && z < mat_2_len)
			{
				test[INDEX(x,y,z,local_size_y,local_size_z)] = 1;
			}

		}//空间中直线方程            
	}
	if (z1 < z2)
	{
		for (int z = (int)floor(z1); z < (int)ceil(z2); z++)
		{
			int x = (int)round((z - z1) / (z2 - z1) * (x2 - x1) + x1);
			int y = (int)round((z - z1) / (z2 - z1) * (y2 - y1) + y1);
			if (x > 0 && x < mat_0_len && y > 0 && y < mat_1_len)
			{
				test[INDEX(x,y,z,local_size_y,local_size_z)] = 1;
			}

		}//空间中直线方程              
	}
	else if (z1 == z2);
	else
	{
		for (int z = (int)floor(z2); z < (int)ceil(z1); z++)
		{
			int x = (int)round((z - z1) / (z2 - z1) * (x2 - x1) + x1);
			int y = (int)round((z - z1) / (z2 - z1) * (y2 - y1) + y1);
			if (x > 0 && x < mat_0_len && y > 0 && y < mat_1_len)
			{
				test[INDEX(x,y,z,local_size_y,local_size_z)] = 1;
			}

		}//空间中直线方程               
	}
	for (int i = 0; i < mat_0_len; i++)
	{
		for (int j = 0; j < mat_1_len; j++)
		{
			for (int k = 0; k < mat_2_len; k++)
			{
				if (test[INDEX(i,j,k,local_size_y,local_size_z)] == 0)
				{
					c1++;//记录现在矩阵中的0个数
				}
			}
		}
	}
	if (c == c1)
	{
		cross = true;
	}
	else
	{
		cross = false;
	}
	return cross;
}//判断两点之间是否可见，即两点连线是否穿过障碍物
