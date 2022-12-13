#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// 原点位移反幂法：
int main(){
	
	// 函数声明
	void MatrixVectorMultiI(double A[], double a[], double b[], int n, int m);
	double VectorsMultiI(double a[], double b[],  int n);
	void Initialization(double A[], int n, int m);
	int MaxComponent(double A[], int n);
	void VecScalarMul(double A[], double B[], double a, int n);
	double VecDiff(double A[], double B[], int n);
	void DoolittleP(double A[], double P[], int n);
	void MatrixLower(double A[], double L[], int n);
	void MatrixUpper(double A[], double U[], int n);
	void LowTriSys(double A[], double x[], double y[], int n);
	void UpTriSys(double A[], double x[], double y[], int n);
	double InvPowerMethodB(double A[], double z[], double VecY[], int n, double eps, int K, double lambda0);
 		
	int i, j, k, n, K;

	// 输入矩阵的阶 
	scanf("%d", &n);
	double A[n*n], z[n], y[n], lambda, lambda0, eps;
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
		
	// 输入 lambda0	
	scanf("%lf", &lambda0);
	
	// 输入计算精度及最大迭代步数a
	scanf("%d %lf", &K, &eps); 

	// 调用反幂法 
	lambda=InvPowerMethodB(A, z, y, n, eps, K, lambda0);
	printf("%lf\n", lambda);
	for(i=0; i<n; i++){
		printf("%lf ", z[i]);
	}
	printf("\n");
	
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

// 计算向量乘法I c=ab^T 
double VectorsMultiI(double a[], double b[],  int n){
	int i;
	double c=0;
	for(i=0; i<n; i++){
		c=c+a[i]*b[i];
	}
	
	return c;
}

// 数组初始化 
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


// 带有行交换 Doolittle 算法做 LU 分解
void DoolittleP(double A[], double P[], int n){
	int i, j, k, m;
	double temp;
	int VecP[n];
	
	// 对 VecP 初始化
	for(i=0; i<n; i++){ VecP[i]=i;}  
	
	for(i=0; i<n; i++){
		// 选主元
		m=i; temp=fabs(A[i*n+i]);
		for(k=i+1; k<n; k++){
			if(temp<fabs(A[k*n+i])){
				temp=fabs(A[k*n+i]);
				m=k;
			}
		} 
		
		// 换主元
		for(k=0; k<n; k++){
			temp=A[m*n+k]; A[m*n+k]=A[i*n+k];
			A[i*n+k]=temp;
		} 
		// 更新 P 
		temp=VecP[m]; VecP[m]=VecP[i]; VecP[i]=temp; 	
	

		// 生成矩阵 U 并存储到 A 的上三角部分 
		for(k=i; k<n; k++){
			for(j=0; j<=i-1; j++){
				A[i*n+k]=A[i*n+k]-A[i*n+j]*A[j*n+k]; 
			}
		}
	
		// 生成矩阵 L 并存储到 A 的下三角部分
		for(k=i+1; k<n; k++){
			for(j=0; j<=i-1; j++){
				A[k*n+i]=A[k*n+i]-A[k*n+j]*A[j*n+i];
			}
			A[k*n+i]=A[k*n+i]/A[i*n+i];
		}	 
	} 
	
	// 生成行变换矩阵 P 
	for(i=0; i<n; i++){
		P[i*n+VecP[i]]=1;
	} 
} 

// 获取单位下三角部分
void MatrixLower(double A[], double L[], int n){
	int i, j;
	
	for(i=0; i<n; i++){	L[i*n+i]=1;}
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++){
			L[i*n+j]=0;
			L[j*n+i]=A[j*n+i];
		}
	}
}

// 获取上三角部分 
void MatrixUpper(double A[], double U[], int n){
	int i, j;
	
	for(i=0; i<n; i++){	U[i*n+i]=A[i*n+i];}
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++){
			U[j*n+i]=0;
			U[i*n+j]=A[i*n+j];
		}
	}
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


// 上三角矩阵求解 
void UpTriSys(double A[], double x[], double y[], int n){
	int i, j;
	for(i=n-1; i>=0; i--){
		x[i]=y[i]/A[i*n+i];
		for(j=0; j<=i-1; j++){
			y[j]=y[j]-x[i]*A[j*n+i];
		}
	} 
}

// 原点位移反幂法：lambda_N<lambda_N-1
//	利用带行变换的 Doolittle 分解完成预处理 
//	返回按模最大特征值 lambda
//	输入：
//		A[]: n-by-n 实矩阵 
//		z[]: 初始向量
//		VecY[]: 计算所得特征向量 
//		eps: 计算精度 
//		K: 最大迭代步数 
//      lambda0: 初始近似 
double InvPowerMethodB(double A[], double z[], double VecY[], int n, double eps, int K, double lambda0){
	int i, j, k=1;
	double err=1.;
	double m, temp, Pz[n], tempx[n];
	double P[n*n], L[n*n], U[n*n];
	Initialization(P, n, n);
	Initialization(L, n, n);
	Initialization(U, n, n);	

	// 基于行变换的 Doolittle 分解 A 的预处理		
	for(i=0; i<n; i++){	A[i*n+i]=A[i*n+i]-lambda0; }
	DoolittleP(A, P, n);
	MatrixLower(A, L, n);
	MatrixUpper(A, U, n);

	// z=Ay
	MatrixVectorMultiI(P, z, Pz, n, n);	
	LowTriSys(L, tempx, Pz, n);
	UpTriSys(U, VecY, tempx, n);
		
	// m=max(y_k)
	j=MaxComponent(VecY, n);
	m=VecY[j];
	temp=m;

	k=k+1;
	
	// z=y/m
	VecScalarMul(VecY, z, 1./m, n);
				
	while(k<K && err>eps){
		// z=Ay
		MatrixVectorMultiI(P, z, Pz, n, n);	
		LowTriSys(L, tempx, Pz, n);
		UpTriSys(U, VecY, tempx, n);
		
		// m=max(y_k)
		j=MaxComponent(VecY, n);
		m=VecY[j];

		// 计算误差: 方式一 
		err=fabs(1./m-1./temp); 	
				
		// 更新：z=y/m
		VecScalarMul(VecY, z, 1./m, n);
		
		// 更新 
		temp=m;	
		
		k=k+1;	
	}
	
	return 1/m+lambda0;
}



