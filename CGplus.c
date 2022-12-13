#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// 计算矩阵向量相乘I b=Aa
void MatrixVectorMultiI(double A[], double a[], double b[], int n, int m){
	int i, j;
	double c;
	for(i=0; i<n; i++){
		c=0;
		for(j=0; j<m; j++){
			c=c+A[i*m+j]*a[j];
		}
		b[i]=c;
	}
} 

// 求向量的无穷范数
double VecNorm_inf(double a[], int n){
	int i;
	double temp, norm=fabs(a[0]);
	
	for(i=1; i<n; i++){
		temp=fabs(a[i]);
		if(norm<temp){	norm=temp; }
	}
	
	return norm;
} 

// 求向量的二范数 
double VecNorm_2(double a[], int n){
	int i;
	double temp, norm=0;
	
	for(i=0; i<n; i++){
		norm=norm+a[i]*a[i];
	}
	
	return norm;
} 


// 计算向量乘法I c=ab^T 
double VectorsMultiI(double a[], double b[],  int n){
	int i;
	double c=0;
	for(i=0; i<n; i++){
		c=c+a[i]*b[i];
	}
	
	return c;
}

// 最速下降法求解 Ax=b, A 为对称正定矩阵
int SteepestDecentSolver(double A[], double b[], double x[], int n, double epsilon, int K){
	double r[n], Ar[n];
	double ak;
	int i, j, k=0;
	double err;
	
	// 计算初始残差 r0=b-Ax0 
	MatrixVectorMultiI(A, x, Ar, n, n);
	for(i=0; i<n; i++){
		r[i]=b[i]-Ar[i];
	}
	err=VecNorm_inf(r, n);
	
	while(err>epsilon && k<K ){
		// 更新 Ar
		MatrixVectorMultiI(A, r, Ar, n, n); 
		
		// 计算 alpha_k
		ak=VectorsMultiI(r, r, n)/VectorsMultiI(Ar, r, n); 
		
		// 更新迭代解 
		for(i=0; i<n; i++){	x[i]=x[i]+ak*r[i];	}
		
		// 更新残差
		for(i=0; i<n; i++){	r[i]=r[i]-ak*Ar[i];	} 
		
		k=k+1;
		err=VecNorm_inf(r, n);
	}
	
	return k;	
}

// 极小残量法求解 Ax=b, A 为对称正定矩阵
int MinResidualSolver(double A[], double b[], double x[], int n, double epsilon, int K){
	double r[n], Ar[n];
	double ak;
	int i, j, k=0;
	double err;
	
	// 计算初始残差 r0=b-Ax0 
	MatrixVectorMultiI(A, x, Ar, n, n);
	for(i=0; i<n; i++){
		r[i]=b[i]-Ar[i];
	}
	err=VecNorm_inf(r, n);
	
	while(err>epsilon && k<K ){
		// 更新 Ar
		MatrixVectorMultiI(A, r, Ar, n, n); 
		
		// 计算 alpha_k
		ak=VectorsMultiI(Ar, r, n)/VectorsMultiI(Ar, Ar, n); 
		
		// 更新迭代解 
		for(i=0; i<n; i++){	x[i]=x[i]+ak*r[i];	}
		
		// 更新残差
		for(i=0; i<n; i++){	r[i]=r[i]-ak*Ar[i];	} 
		
		k=k+1;
		err=VecNorm_inf(r, n);
	}
	
	return k;	
}

// 极共轭梯度求解 Ax=b, A 为对称正定矩阵
int CGSolver(double A[], double b[], double x[], int n, double epsilon, int K){
	double r[n], p[n], Ap[n];
	double ak, betak;
	int i, j, k=0;
	double err, temp;
	
	// 计算初始残差 r0=b-Ax0 
	MatrixVectorMultiI(A, x, Ap, n, n);
	for(i=0; i<n; i++){
		r[i]=b[i]-Ap[i];
		p[i]=r[i];
	}
	err=VecNorm_inf(r, n);
	//err=VecNorm_2(r, n);	
	
	while(err>epsilon && k<K ){
		// 更新 Ap
		MatrixVectorMultiI(A, p, Ap, n, n); 
		
		// 计算 Ap_k
		temp=VectorsMultiI(Ap, p, n);
		
		// 计算 alpha_k
		ak=VectorsMultiI(r, p, n)/temp; 
		
		// 更新迭代解 
		for(i=0; i<n; i++){	x[i]=x[i]+ak*p[i];	}
		
		// 更新残差
		for(i=0; i<n; i++){	r[i]=r[i]-ak*Ap[i];	} 

		// 计算 beta_k
		betak=-VectorsMultiI(r, Ap, n)/temp; 

		// 更新 p
		for(i=0; i<n; i++){	p[i]=r[i]+betak*p[i];	} 
						
		k=k+1;
		err=VecNorm_inf(r, n);
		//err=VecNorm_2(r, n);	
	}
	
	return k;	
}
 
int main(){
	int i, j;
	
	// 测试数据
	double A[16]={13, 9, 9, 11, 9, 11, 6, 9, 9, 6, 12, 7, 11, 9, 7, 10};
	double b[4]={42,35,34,37};
	double epsilon=0.000001;
	int K=2000;

	// 作业
	int k;
	double x[4]={0, 0, 0, 0};
	
	// 最速下降法 
	//k=SteepestDecentSolver(A, b, x, 4, epsilon, K);
	
	// 极小残量法 
	//k=MinResidualSolver(A, b, x, 4, epsilon, K);	
	
	// 共轭梯度法
	k=CGSolver(A, b, x, 4, epsilon, K);
	
	// 输出结果 
	printf("step = %d\n", k);
	printf("x =\n");
	for(i=0; i<4; i++){
		printf("%lf\n", x[i]);
	}	 
		
	return 0;
}

