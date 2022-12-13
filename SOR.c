#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// SOR 迭代法 求解 Ax=b, 
//	停机准则：|x_{k+1}-x_K|<\epsilon 或 迭代次数大于 K 
//	本次作业停机准则按 \infty-范数 
int SORSolver(double A[], double b[], double x[], double x0[], int n, double epsilon, int K, double omega){
	int i, j, k=0;
	double err=1, temp;
	
	while(err>epsilon && k<K){// 停机准则 
		for(i=0; i<n; i++){
			temp=0;
			for(j=0; j<=i-1; j++){
				temp=temp+A[i*n+j]*x[j];
			}
			
			for(j=i+1; j<n; j++){
				temp=temp+A[i*n+j]*x0[j];
			}
			
			x[i]=(1-omega)*x0[i]-omega*(temp-b[i])/A[i*n+i];
		}
		
		k=k+1;
		
		// 计算残差：\infty-范数
		err=fabs(x[0]-x0[0])/(1+fabs(x[0]));
		for(i=1; i<n; i++){
			temp=fabs(x[i]-x0[i])/(1+fabs(x[i]));
			if(temp>err){ err=temp;}
		} 
		
		// 更新
		for(i=0; i<n; i++){	x0[i]=x[i]; } 
	}
	
	return k;
} 


// JOR 算法求解 Ax=b 
//		停机准则：|x_{k+1}-x_K|<\epsilon 或 迭代次数大于 K 
//		本次作业停机准则按 \infty-范数 
int JORSolver(double A[], double b[], double x[], double x0[], int n, double epsilon, int K, double omega){
	int i, j, k=0;
	double err=1, temp;
	
	while(err>epsilon && k<K){// 停机准则 
		for(i=0; i<n; i++){
			temp=0;
			for(j=0; j<=i-1; j++){
				temp=temp+A[i*n+j]*x0[j];
			}
			
			for(j=i+1; j<n; j++){
				temp=temp+A[i*n+j]*x0[j];
			}
			
			x[i]=(1-omega)*x0[i]-omega*(temp-b[i])/A[i*n+i];
		}
		
		k=k+1;
		
		// 计算残差：\infty-范数
		err=fabs(x[0]-x0[0])/(1+fabs(x[0]));
		for(i=1; i<n; i++){
			temp=fabs(x[i]-x0[i])/(1+fabs(x[i]));
			if(temp>err){ err=temp;}
		} 
		
		// 更新
		for(i=0; i<n; i++){	x0[i]=x[i]; } 
	}
	
	return k;
} 


int main(){
	int i, j;
	
	// 测试数据
	double A[16]={13, 9, 9, 11, 9, 11, 6, 9, 9, 6, 12, 7, 11, 9, 7, 10};
	double b[4]={42,35,34,37};
	double epsilon=0.000001;
	int K=5000;


	// 作业
	double omega=0.5;
	int k;
	double x[4], x0[4]={0, 0, 0, 0};
	
	// SOR 迭代法 
	k=SORSolver(A, b, x, x0, 4, epsilon, K, omega);
	
	// JOR 迭代法
	//k=JORSolver(A, b, x, x0, 4, epsilon, K, omega);
		 
	printf("step = %d\n", k);
	printf("x =\n");
	for(i=0; i<4; i++){
		printf("%lf\n", x[i]);
	}	 



/*
	// 检查不同松弛因子对 SOR 迭代效率的影响 
	double omega[100];
	int steps[100];
	for(i=0; i<100; i++){ omega[i]=0.02*(i+1);	}

	for(i=0; i<100; i++){
		// SOR 迭代法并检验不同松弛因子对迭代效率的影响 
		double x[4], x0[4]={0, 0, 0, 0};	
		steps[i]=SORSolver(A, b, x, x0, 4, epsilon, K, omega[i]);
	} 
	
	// 输出
	for(i=0; i<100; i++){
		printf("omega = %lf \t step = %d \n", omega[i], steps[i]);
	}
*/

/*	
	// 检查不同松弛因子对 JOR 迭代效率的影响 
	double omega[100];
	int steps[100];
	for(i=0; i<100; i++){ omega[i]=0.01*(i+1);	}

	for(i=0; i<100; i++){
		// SOR 迭代法并检验不同松弛因子对迭代效率的影响 
		double x[4], x0[4]={0, 0, 0, 0};	
		steps[i]=JORSolver(A, b, x, x0, 4, epsilon, K, omega[i]);
	} 
	
	// 输出
	for(i=0; i<100; i++){
		printf("omega = %lf \t step = %d \n", omega[i], steps[i]);
	}
*/
		
	return 0;
}




