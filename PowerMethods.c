#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
	
	// 函数声明
	void MatrixVectorMultiI(double A[], double a[], double b[], int n, int m);
	void Initialization(double A[], int n, int m);
	int MaxComponent(double A[], int n);
	void VecScalarMul(double A[], double B[], double a, int n);
	double VecDiff(double A[], double B[], int n);
	double PowerMethodA(double A[], double z[], double VecY[], int n, double eps, int K);
 
	int i, j, k, n, K;

	// 输入矩阵的阶 
	scanf("%d", &n);
	double A[n*n], z[n], y[n], lambda, eps;
	Initialization(A, n, n);
	Initialization(z, n, 1);	
	Initialization(y, n, 1);	

	// 输入矩阵元素
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			scanf("%lf", &A[i*n+j]);
		}
	} 
	
	// 输入初始向量
	for(i=0; i<n; i++){
		scanf("%lf", &z[i]);
	} 

	// 输入计算精度及最大迭代步数
	scanf("%lf %d", &eps, &K); 
	
	// 调用基本乘幂法 
	lambda=PowerMethodA(A, z, y, n, eps, K);
	
	// 输出结果
	printf("%.4lf\n", lambda);
	for(i=0; i<n; i++){
		printf("%.4lf ", y[i]);
	}

	return 0;
}

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

// 初始化 
void Initialization(double A[], int n, int m){
	int i;
	for(i=0; i<n*m; i++){ A[i]=0.; }	
}

// 找到向量 A 绝对值最大的分量 
int MaxComponent(double A[], int n){
	int i, k;
	double temp;
	temp=fabs(A[0]);
	k=0;
	
	for(i=1; i<n; i++){
		if(temp<fabs(A[i])){
			k=i;
			temp=fabs(A[i]);
		}
	}
	
	return k;
}


// 向量标量乘法
void VecScalarMul(double A[], double B[], double a, int n){
	int i;
	for(i=0; i<n; i++){
		B[i]=A[i]*a;
	}
} 

// 计算向量差的范数
double VecDiff(double A[], double B[], int n){
	int i;
	double temp=0;
	for(i=0; i<n; i++){
		temp=fabs(A[i]-B[i]);
	}
	
	return temp;
} 


// 基本乘幂法：lambda_1>lambda_2
//	返回按模最大特征值 lambda
//	输入：
//		A[]: n-by-n 实矩阵 
//		z[]: 初始向量及计算所得特征向量 
//		eps: 计算精度 
//		K: 最大迭代步数 
double PowerMethodA(double A[], double z[], double VecY[], int n, double eps, int K){
	int i, j, k=1;
	double err=1.;
	double m, temp;
	
	// y=Az
	MatrixVectorMultiI(A, z, VecY, n, n);
		
	// m=max(y_k)
	j=MaxComponent(VecY, n);
	m=VecY[j];
	temp=m;
			
	// z=y/m
	VecScalarMul(VecY, z, 1./m, n);
				
	while(k<K && err>eps){
		// y=Az
		MatrixVectorMultiI(A, z, VecY, n, n);
	
		// m=max(y_k)
		j=MaxComponent(VecY, n);
		m=VecY[j];
	
		// 计算误差: 方式一 
		err=fabs(m-temp);
				
		// 更新：z=y/m
		VecScalarMul(VecY, z, 1./m, n);
		
		// 更新 
		temp=m;	
		
		k=k+1;	
	}
	
	return m;
}



