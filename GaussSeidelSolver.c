#include<stdio.h>
#include<math.h>

// Gauss-Seidel 迭代法 求解 Ax=b, 
int main(){

	// 相关函数声明 	
	int GaussSeidelSolver(double A[], double b[], double x[], double x0[], int n, double epsilon, int K);	
	//int GaussSeidelSolverA(double A[], double b[], double x[], double x0[], int n, double epsilon, int K);
	//int GaussSeidelSolverB(double A[], double b[], double x[], double x0[], int n, double epsilon, int K);	
	
	int i, j;
	// 测试数据
	//double A[9]={8, -1, 1, 2, 10, -1, 1, 1, -5}, b[3]={1, 4, 3}, x[3], x0[3]={0.125, 0.4, -0.6};
	double A[16]={12,9,4,3,2,11,1,5,7,2,13,1,-5,2,-7,-10};
	double b[4]={10,20,1,15};
	double x[4],x0[4];
	double epsilon=0.000001;
	int K=1000, k;
	
	
	// G-S 迭代法 
	k=GaussSeidelSolver(A, b, x, x0, 4, epsilon, K);
	//k=GaussSeidelSolverA(A, b, x, x0, 3, epsilon, K);	
	//k=GaussSeidelSolverB(A, b, x, x0, 3, epsilon, K);
		
	// 输出
	printf("step = %d \n", k);
	printf("x=\n");
	for(i=0; i<4; i++){
		printf("%lf\n", x[i]);
	} 
	
	return 0;
}

// 下三角矩阵求解 
void LowTriSys(double A[], double x[], double y[], int n){
	int i, j;
	for(i=0; i<n; i++){
		x[i]=y[i]/A[i*n+i];
		for(j=i+1; j<n; j++){
			y[j]=y[j]-x[i]*A[j*n+i];
		}
	} 
}

// 下三角矩阵的逆
void InverseL(double A[], double IL[], int n){
	int i, j, k;
	double temp;
	
	// 求 IL 的第 i 列元素
	for(i=0; i<n; i++){
		IL[i*n+i]=1/A[i*n+i];
		// 求 IL 的第 k 行第 i 列元素
		for(k=i+1; k<n; k++){
			temp=0;
			for(j=i; j<k; j++){
				temp=temp+A[k*n+j]*IL[j*n+i];
			}
			IL[k*n+i]=-temp/A[k*n+k];
		}	
	} 
} 


//	停机准则：|x_{k+1}-x_K|<\epsilon 或 迭代次数大于 K 
//	本次作业停机准则按 \infty-范数 
int GaussSeidelSolver(double A[], double b[], double x[], double x0[], int n, double epsilon, int K){
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
			
			x[i]=-(temp-b[i])/A[i*n+i];
		}
		
		k=k+1;
		
		// 计算残差：\infty-范数
		err=fabs(x[0]-x0[0]);
		for(i=1; i<n; i++){
			//temp=fabs(x[i]-x0[i]);
			if(fabs(x[i]-x0[i])>err){ err=fabs(x[i]-x0[i]);}
		} 
		
		// 更新
		for(i=0; i<n; i++){	x0[i]=x[i]; } 
	}
	
	return k;
} 

//	停机准则：|x_{k+1}-x_K|<\epsilon 或 迭代次数大于 K 
//	本次作业停机准则按 \infty-范数 
int GaussSeidelSolverA(double A[], double b[], double x[], double x0[], int n, double epsilon, int K){
	int i, j, k=0;
	double err=1, temp;
	double tempb[n];
	
	while(err>epsilon && k<K){// 停机准则 
		// 计算 tempb=Ux+b 
		for(i=0; i<n; i++){
			temp=0;
			
			for(j=i+1; j<n; j++){
				temp=temp+A[i*n+j]*x0[j];
			}
			
			tempb[i]=-temp+b[i];
		}
		
		// 求解 (D-L)x=tempb 
		LowTriSys(A, x, tempb, n);
		 
		k=k+1;
		
		// 计算残差：\infty-范数
		err=fabs(x[0]-x0[0]);
		for(i=1; i<n; i++){
			temp=fabs(x[i]-x0[i]);
			if(temp<err){ err=temp;}
		} 
		
		// 更新
		for(i=0; i<n; i++){	x0[i]=x[i]; } 
	}
	
	return k;	
}

//	停机准则：|x_{k+1}-x_K|<\epsilon 或 迭代次数大于 K 
//	停机准则按 \infty-范数 
//	做预处理 (D-L)^{-1}U, (D-L)^{-1}b 
int GaussSeidelSolverB(double A[], double b[], double x[], double x0[], int n, double epsilon, int K){
	int i, j, k=0;
	double err=1, temp;
	double tempb[n];
	
	// 预处理
	double IDL[n*n], DLU[n*n], DLb[n]; 
	InverseL(A, IDL, n);
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			temp=0;
			for(k=0; k<fmin(i+1, j); k++){
				temp=temp+IDL[i*n+k]*A[k*n+j];
			}
			DLU[i*n+j]=temp;
		}
		
		temp=0;
		for(j=0; j<=i; j++){
			temp=temp+IDL[i*n+j]*b[j];
		}
		DLb[i]=temp;
	}
	
	k=0;
	while(err>epsilon && k<K){// 停机准则 
		// 计算 tempb=Ux+b 
		for(i=0; i<n; i++){
			temp=0;
			for(j=0; j<n; j++){
				temp=temp+DLU[i*n+j]*x0[j];
			}
			
			x[i]=-temp+DLb[i];
		}
		 
		k=k+1;
		
		// 计算残差：\infty-范数
		err=fabs(x[0]-x0[0]);
		for(i=1; i<n; i++){
			temp=fabs(x[i]-x0[i]);
			if(temp<err){ err=temp;}
		} 
		
		// 更新
		for(i=0; i<n; i++){	x0[i]=x[i]; } 
	}
	
	return k;	
}

