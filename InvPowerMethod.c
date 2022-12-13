#include<stdio.h>
#include<math.h>

int main(){
	double InvPowerMethodr1(double A[],double y[],double z[],int n,int K,double epsilon,double p);
	double InvPowerMethodr2(double A[],double y[],double z[],int n,int K,double epsilon,double p);
	void MatPrintf(double Q[],int m,int n);
	void MatCon(double Q[],double z[],int n);
	
	int n,i,K;
	//scanf("%d",&n);
	//printf("%d\n",n);
	
	//double A[n*n],z[n],y[n];
	double A[9]={2,-1,0,0,2,-1,0,-1,2},z[3]={0,0,1},y[3];
	
	double epsilon,m,p;
	n=3;  K=50;  epsilon=0.01;  p=2.93;
	/*
	//输入 
	for(i=0;i<n*n;i++){
		scanf("%lf",A+i);
	}
	//MatPrintf(A,n,n);
	for(i=0;i<n;i++){
		scanf("%lf",z+i);
	}
	//MatPrintf(z,1,n); 
	scanf("%lf",&p);
	//printf("%.8lf\n",p);	
	scanf("%d",&K);
	//printf("%d\n",K);
	scanf("%lf",&epsilon);
	//printf("%.8lf\n",epsilon);
	*/
	//输出
	double A0[n*n],z0[n],y0[n];
	MatCon(A,A0,n*n);
	MatCon(z,z0,n);
	
	m=InvPowerMethodr2(A,y,z,n,K,epsilon,p);
	m=1/m+p;
	printf("%lf\n",m);
	
	m=InvPowerMethodr1(A0,y0,z0,n,K,epsilon,p);
	for(i=0;i<3;i++){
		printf("%lf ",z0[i]);
	}
	
}

//反幂法  
//r=fabs(m-t)
double InvPowerMethodr1(double A[],double y[],double z[],int n,int K,double epsilon,double p){
	int Vecmax(double v[],int n);//向量的绝对值最大的分量
	void Doolittle(double A[], int n);// Doolittle 算法做 LU 分解
	void MatrixUpper(double A[], double U[], int n);// 获取上三角部分 
	void MatrixLower(double A[], double L[], int n);// 获取单位下三角部分
	void LUSolver(double L[],double U[],double x[],double b[],int n);//求解 LUx=b 
	void MatPrintf(double Q[],int m,int n);
	
	int i,k,j;
	double m,r,t;
	double L[n*n],U[n*n];
	
	//初始化
	for(i=0;i<n;i++){
		A[i*n+i]-=p;
	} 
	
	Doolittle(A,n);
	MatrixUpper(A,U,n);
	MatrixLower(A,L,n);
	//MatPrintf(L,n,n);
	//MatPrintf(U,n,n);
	LUSolver(L,U,y,z,n);
	j=Vecmax(y,n);
	m=y[j];
	t=m;
	
	for(k=2;k<=K;k++){
		
		for(i=0;i<n;i++){
			z[i]=y[i]/m;
		}//更新 z
		
		LUSolver(L,U,y,z,n);//求解Ay=z
		
		j=Vecmax(y,n);
		m=y[j];//绝对值最大的分量
		
		r=fabs(m-t);//残差 
		
		if(r<epsilon){
			break;
		}
		
		t=m;//更新t
	}
	
	return m; 
}

//反幂法
//r=fabs(1/m-1/t)
double InvPowerMethodr2(double A[],double y[],double z[],int n,int K,double epsilon,double p){
	int Vecmax(double v[],int n);//向量的绝对值最大的分量
	void Doolittle(double A[], int n);// Doolittle 算法做 LU 分解
	void MatrixUpper(double A[], double U[], int n);// 获取上三角部分 
	void MatrixLower(double A[], double L[], int n);// 获取单位下三角部分
	void LUSolver(double L[],double U[],double x[],double b[],int n);//求解 LUx=b 
	void MatPrintf(double Q[],int m,int n);
	void MatCon(double Q[],double z[],int n);
	
	int i,k,j;
	double m,r,t;
	double L[n*n],U[n*n],P[n*n];
	
	//初始化
	MatCon(A,P,n*n);
	for(i=0;i<n;i++){
		A[i*n+i]-=p;
	} 
	
	Doolittle(A,n);
	MatrixUpper(A,U,n);
	MatrixLower(A,L,n);
	//MatPrintf(L,n,n);
	//MatPrintf(U,n,n);
	LUSolver(L,U,y,z,n);
	j=Vecmax(y,n);
	m=y[j];
	t=m;
	
	for(k=2;k<=K;k++){
		
		for(i=0;i<n;i++){
			z[i]=y[i]/m;
		}//更新 z
		
		LUSolver(L,U,y,z,n);//求解Ay=z
		
		j=Vecmax(y,n);
		m=y[j];//绝对值最大的分量
		
		r=fabs(1./m-1./t);//残差 
		
		if(r<epsilon){
			break;
		}
		
		t=m;//更新t
	}
	MatCon(P,A,n*n);
	return m; 
}

// Doolittle 算法做 LU 分解
void Doolittle(double A[], int n){
	int i, j, k;
	
	for(i=0; i<n; i++){
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
} 

//求解 LUx=b 
void LUSolver(double L[],double U[],double x[],double b[],int n){
	void UpTriSys(double U[],double x[],double y[],int n);
	void DwTriSys(double L[],double x[],double y[],int n);
	double p[n];
	int i;
	DwTriSys(L,x,b,n);
	for(i=0;i<n;i++){
		p[i]=x[i];
	}
	UpTriSys(U,x,p,n);	
}

//U为上三角矩阵，求解 Ux=y 
void UpTriSys(double U[],double x[],double y[],int n){
	int i,j,m;
	double p[n];
	for(j=0;j<n;j++)
		p[j]=y[j];
	for(j=n;j>=1;j--)
	{
		m=(j-1)*(n+1);
		x[j-1]=p[j-1]/U[m];
		for(i=j-1;i>=1;i--)
			p[i-1]=p[i-1]-U[(i-1)*n+j-1]*x[j-1];
	}
}

//L为下三角矩阵，求解 Lx=y 
void DwTriSys(double L[],double x[],double y[],int n){
	int i,j,m;
	double p[n];
	for(j=0;j<n;j++)
		p[j]=y[j];
	for(j=1;j<=n;j++)
	{
		m=(j-1)*(n+1);
		x[j-1]=p[j-1]/L[m];
		for(i=j+1;i<=n;i++)
			p[i-1]=p[i-1]-L[(i-1)*n+j-1]*x[j-1];
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
 
//矩阵乘法
void MatMul(double A[],double B[],double y[],int m,int n,int p){
	int i,j,k;
	for(i=1;i<=m;i++)
		for(j=1;j<=p;j++)
			y[(i-1)*p+j-1]=0;
	for(i=1;i<=m;i++)
		for(j=1;j<=p;j++)
			for(k=1;k<=n;k++)
				y[(i-1)*p+j-1]+=A[(i-1)*n+k-1]*B[(k-1)*p+j-1];
}

//向量的绝对值最大的分量
int Vecmax(double v[],int n){
	double t;
	int i,j;
	j=0;
	t=fabs(v[0]);
	for(i=0;i<n;i++){	
		if(t<fabs(v[i])){
			t=fabs(v[i]);
			j=i;
		}
	}
	
	return j;
} 

//矩阵输出
void MatPrintf(double Q[],int m,int n){
	int i,j;
	for(i=1;i<=m;i++)
	{
		for(j=1;j<=n;j++){
			printf("%lf\t",Q[(i-1)*n+j-1]);
		}
		printf("\n");
	}
} 

//矩阵继承
void MatCon(double Q[],double z[],int n){
	int i;
	for(i=0;i<n;i++)
		z[i]=Q[i];
} 




