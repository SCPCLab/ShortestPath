#include<slave.h>
#include<math.h>
#include<stdio.h>
#define size_x 25
#define size_y 25
#define size_z 100
#define local_size_x 5
#define local_size_y 5
#define local_size_z 50
#define halo_size 1

#define INDEX(i,j,k,ldy,ldz)  ( (i)*(ldy)*(ldz) + (j)*(ldz)  + (k) )

extern double* local_gradmat;
extern double* local_H;
extern double* local_gradmatnew;
extern int* local_mat_3d;

__thread_local double local_old[(local_size_x+2*halo_size)*(local_size_y+2*halo_size)*(local_size_z+2*halo_size)];
__thread_local double local_new[local_size_x*local_size_y*local_size_z];
__thread_local double llocal_H[local_size_x*local_size_y*local_size_z];
__thread_local int mat_3d[local_size_x*local_size_y*local_size_z]; 
void readMat(){
	volatile int get_reply;
	volatile int put_reply;
    	int t_id = athread_get_id(-1);
	int cx,cy,cz;
	int ldy=local_size_y;   int ldz=local_size_z; 
    	cx=size_x/local_size_x; //x方向迭代次数
    	cy=size_y/local_size_y; //y方向迭代次数
    	cz=size_z/local_size_z;
	//start 10 5 5   5 4 100
    	//按z,y,x顺序得到进程对应的坐标
   	const int x=t_id/(cy*cz);//线程x坐标
    	const int y=(t_id/cz)%cy;//线程y坐标
    	const int z=t_id%cz;
	for(int i=0;i<local_size_x;i++){
              	  	for(int j=0;j<local_size_y;j++){
                        		get_reply=0;
                       		 athread_get(PE_MODE,&(local_mat_3d[INDEX(x*local_size_x+i,y*local_size_y+j,z*local_size_z,size_y,size_z)]), 
                        		&(mat_3d[INDEX(i,j,0,local_size_y,local_size_z)]),local_size_z*sizeof(int),
                        		&get_reply, 0, 0, 0); 
                        		while(get_reply!=1);
                		}   
	}
}
void fun1(int data[3]){
    	volatile int get_reply;
	volatile int put_reply;
    	int t_id = athread_get_id(-1);
    	double rat=0.1;

    	int icount=data[2];
    	int gx=data[0],gy=data[1];

    	int cx,cy,cz;
	int ldy=local_size_y+2*halo_size;   int ldz=local_size_z+2*halo_size; 
    	cx=size_x/local_size_x; //x方向迭代次数
    	cy=size_y/local_size_y; //y方向迭代次数
    	cz=size_z/local_size_z;
	//start 10 5 5   5 4 100
    	//按z,y,x顺序得到进程对应的坐标
   	const int x=t_id/(cy*cz);//线程x坐标
    	const int y=(t_id/cz)%cy;//线程y坐标
    	const int z=t_id%cz;
	//单个从核一次计算量
	for(int i=0;i<local_size_x+2*halo_size;i++){
		for(int j=0;j<local_size_y+2*halo_size;j++){
			//每次传1*1*(50+2)
			get_reply=0;
			athread_get(PE_MODE,&(local_gradmat[INDEX(x*local_size_x+i,y*local_size_y+j,z*local_size_z,size_y+2*halo_size,size_z+2*halo_size)]),
			&(local_old[INDEX(i,j,0,ldy,ldz)]),(local_size_z+2*halo_size)*sizeof(double),
			&get_reply, 0, 0, 0); 
			while(get_reply!=1);
		}   
	}  	

	

	//计算部分
	for(int i=0;i<local_size_x;i++){ 
		if(x==0&&gx==0&&i==0||x==cx-1&&gx==3&&i==local_size_x-1) continue; //如果位于原始数据块的边界,略过
		for(int j=0;j<local_size_y;j++){
			if(y==0&&gy==0&&j==0||y==cy-1&&gy==3&&j==local_size_y-1) continue; //同理
			for(int k=0;k<local_size_z;k++){
				if(z==0&&k==0||z==cz-1&&k==local_size_z-1) continue;
				local_new[INDEX(i,j,k,local_size_y,local_size_z)] = 
				local_old[INDEX(i+halo_size,j+halo_size,k+halo_size,ldy,ldz)] * (1 + rat)
				+ (local_old[INDEX(i - 1+halo_size,j+halo_size,k+halo_size,ldy,ldz)]
				+ local_old[INDEX(i + 1+halo_size,j+halo_size,k+halo_size,ldy,ldz)]
				+ local_old[INDEX(i+halo_size,j - 1+halo_size,k+halo_size,ldy,ldz)]
				+ local_old[INDEX(i+halo_size,j + 1 + halo_size,k+halo_size,ldy,ldz)]
				+ local_old[INDEX(i+halo_size,j+halo_size,k - 1+halo_size,ldy,ldz)]
				+ local_old[INDEX(i+halo_size,j+halo_size,k + 1+halo_size,ldy,ldz)]) * rat; 
			}
		}
	}
	for(int i=0;i<local_size_x;i++){ 
                for(int j=0;j<local_size_y;j++){
                        for(int k=0;k<local_size_z;k++){
				local_old[INDEX(i+halo_size,j+halo_size,k+halo_size,ldy,ldz)] = local_new[INDEX(i,j,k,local_size_y,local_size_z)]*mat_3d[INDEX(i,j,k,local_size_y,local_size_z)];
				
				if(local_old[INDEX(i+halo_size,j+halo_size,k+halo_size,ldy,ldz)]>1&&llocal_H[INDEX(i,j,k,local_size_y,local_size_z)]==0){
					llocal_H[INDEX(i,j,k,local_size_y,local_size_z)]=icount;
				}

                        }   
                }   
        }  
	
	//写回数据
	for(int i=0;i<local_size_x;i++){
		for(int j=0;j<local_size_y;j++){
			put_reply = 0;
			//写回数据,每次写回(1*1*50)
			athread_put(PE_MODE,&(local_old[INDEX(i+1,j+1,1,ldy,ldz)]),
			&(local_gradmat[INDEX(x*local_size_x+i+halo_size,y*local_size_y+j+halo_size,z*local_size_z+halo_size,size_y+2*halo_size,size_z+2*halo_size)]),
			local_size_z*sizeof(double),&put_reply, 0, 0);
			while (put_reply != 1);
		}
	}
}
void putH(int icount){
	volatile int put_reply;
    	int t_id = athread_get_id(-1);
    	int cx,cy,cz;
	int ldy=local_size_y;   int ldz=local_size_z; 
    	cx=size_x/local_size_x; //x方向迭代次数
    	cy=size_y/local_size_y; //y方向迭代次数
    	cz=size_z/local_size_z;
   	const int x=t_id/(cy*cz);//线程x坐标
    	const int y=(t_id/cz)%cy;//线程y坐标
    	const int z=t_id%cz;
	for(int i=0;i<local_size_x;i++){
		for(int j=0;j<local_size_y;j++){
			for(int k=0;k<local_size_z;k++){
				if(llocal_H[INDEX(i,j,k,local_size_y,local_size_z)]==0) 
					llocal_H[INDEX(i,j,k,local_size_y,local_size_z)]=icount;
			}
		}
	}
	for(int i=0;i<local_size_x;i++){
		for(int j=0;j<local_size_y;j++){
			put_reply = 0;
			//写回数据,每次写回(1*1*50)
			athread_put(PE_MODE,&(llocal_H[INDEX(i,j,0,ldy,ldz)]),
			&(local_H[INDEX(x*local_size_x+i+halo_size,y*local_size_y+j+halo_size,z*local_size_z+halo_size,size_y+2*halo_size,size_z+2*halo_size)]),
			local_size_z*sizeof(double),&put_reply, 0, 0);
			while (put_reply != 1);
		}
	}
}
