#include<stdio.h>
#include<math.h>

// Gauss-Seidel 迭代法 求解 Ax=b, 
int main(){

	// 相关函数声明 	
	int GaussSeidelSolver(double A[], double b[], double x[], double x0[], int n, double epsilon, int K);	
	
	int i, j;
	// 测试数据
	double A[9]={20,2,3,1,8,1,2,-3,15},b[3]={25,10,14},x[3],x0[3]={0,0,1};
	double epsilon=0.000001;
	int K=1000, k;
	
	
	// G-S 迭代法 
	k=GaussSeidelSolver(A, b, x, x0, 3, epsilon, K);
	/*	
	// 输出
	printf("step = %d \n", k);
	printf("x=\n");
	for(i=0; i<3; i++){
		printf("%lf\n", x[i]);
	} 
	*/
	return 0;
}

//	停机准则：|x_{k+1}-x_K|<\epsilon 或 迭代次数大于 K 
//	本次作业停机准则按 \infty-范数 
int GaussSeidelSolver(double A[], double b[], double x[], double x0[], int n, double epsilon, int K){
	int i, j, k;
	double err, temp;
	
	for(k=0,err=1; err>epsilon && k<K; k++){// 停机准则 
		for(i=0; i<n; i++){
			temp=0;
			for(j=0; j<=i-1; j++){
				temp=temp+A[i*n+j]*x[j];
			}
			
			for(j=i+1; j<n; j++){
				temp=temp+A[i*n+j]*x0[j];
			}
			
			x[i]=-(temp-b[i])/A[i*n+i];
		}
		
		// 计算残差：\infty-范数
		err=fabs(x[0]-x0[0]);
		for(i=1; i<n; i++){
			temp=fabs(x[i]-x0[i]);
			if(temp<err){ err=temp;}
		} 
		
		// 更新
		for(i=0; i<n; i++){	x0[i]=x[i]; } 
		/*
		printf("x_(%d)\n",k+1);
		for(i=0; i<n; i++){
			printf("%lf\t",x[i]);
		}
		printf("\n");
		*/
	}
	
	return k;
} 

