#include<stdio.h>
#include<math.h>

int main(){
	int SORSolver(double A[], double b[], double x[], double x0[], int n, double epsilon, int K,double w);
	int i, j;
	
	// 测试数据
	double A[16]={13,9,9,11,9,11,6,9,9,6,12,7,11,9,7,10},b[4]={42,35,34,37};
	double x[4],x0[4];
	double epsilon=0.000001,w=0.5;
	int K=2000, k;	
	
	// Jacobi 迭代法 
	k=SORSolver(A, b, x, x0, 4, epsilon, K, w);
	
	// 输出
	printf("step = %d \n", k);
	printf("x =\n");
	for(i=0; i<4; i++){
		printf("%lf\n", x[i]);
	} 
	
	return 0;
}

//  SOR迭代法 求解 Ax=b, 
//	停机标准 max(fabs(x[i]-x0[i])/(1+fabs(x[i])))<精度 
int SORSolver(double A[], double b[], double x[], double x0[], int n, double epsilon, int K,double w){
	int i,j,k;
	double err,temp;
	
	for(k=0,err=1; err>epsilon && k<K; k++){// 停机准则 
		for(i=0; i<n; i++){
			temp=0;
			for(j=0; j<=i-1; j++){
				temp=temp+A[i*n+j]*x[j];
			}
			
			for(j=i; j<n; j++){
				temp=temp+A[i*n+j]*x0[j];
			}
			
			x[i]=x0[i]-w*(temp-b[i])/A[i*n+i];
		}
		// 计算  停机标准 max(fabs(x[i]-x0[i])/(1+fabs(x[i])))<精度  
		err=fabs(x[0]-x0[0])/(1+fabs(x[0]));
		for(i=1; i<n; i++){
			temp=fabs(x[i]-x0[i])/(1+fabs(x[i]));
			if(err<temp){ err=temp;}
		} 
		// 更新
		for(i=0; i<n; i++){	x0[i]=x[i]; } 	
		
	}
	return k;
} 

