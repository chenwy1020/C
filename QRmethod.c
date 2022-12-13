#include<stdio.h>
#include<math.h>

int main(){
	//函数声明 
	int QRMethod(double A[],double lambda[],int n,int K,double epsilon);
	void MatPrintf(double Q[],int m,int n);
	double InvPowerMethodr1(double A[],double y[],double z[],int n,int K,double epsilon,double p);
	
	int n=4,k,i;
	double A[16]={13,9,9,11,9,11,6,9,9,6,12,7,11,9,7,10};
	double z[4]={1,0,0,0},y[4];
	double lambda[4];
	double epsilon=0.000001,K=10000;
	
	k=QRMethod(A,lambda,n,K,epsilon);
	printf("QRMethod迭代步数：%d\n",k);
	MatPrintf(lambda,1,n);
	printf("\n"); 
	
	for(i=0;i<n;i++){
		printf("反幂法第 %d 数值特征值：%lf\n",i+1,lambda[i]);
		InvPowerMethodr1(A,y,z,n,K,epsilon,lambda[i]);
		printf("反幂法第 %d 数值特向量：",i+1);
		MatPrintf(z,1,n);
		printf("\n"); 
	}
	
}

//QR基本算法 
int QRMethod(double A[],double lambda[],int n,int K,double epsilon){
	void MatMul(double A[],double B[],double y[],int m,int n,int p);
	void Householder(double A[],double Q[],double R[],int n);
	void MatPrintf(double Q[],int m,int n);
	void MatCon(double Q[],double z[],int n);
	
	double t,r;
	double Q[n*n],R[n*n],A0[n*n];
	int i,k;
	//初始化
	MatCon(A,A0,n*n); 
	for(i=0;i<n;i++){
		lambda[i]=A[i*n+i];
	}

	for(k=2,r=1;k<K && r>epsilon;k++){
		//A=QR
		Householder(A,Q,R,n);
		
		//计算 A_k+1=R*Q 
		MatMul(R,Q,A,n,n,n);
		
		//计算残差
		r=0;
		for(i=1;i<=n;i++){
			t=A[(i-1)*n+i-1];
			r+=fabs(t-lambda[i-1]);
			lambda[i-1]=t;		
		}
		if(r<=epsilon) break;		
	}
	MatCon(A0,A,n*n);
	return k;
	

}

//Householder变换 
void Householder(double A[],double Q[],double R[],int n){
	void MatT(double Q[],int n);
	void MatCon(double Q[],double z[],int n);
	
	int i,j,k,i1,j1;
	double t,t1=0,t2=0,t3=0,v[n];
	double P[n*n];
	
	//初始化 
	MatCon(A,P,n*n); 
	for(i1=1;i1<=n;i1++)
		for(j1=1;j1<=n;j1++)
		{
			Q[(i1-1)*n+j1-1]=0;
			Q[(i1-1)*(n+1)]=1; 
			R[(i1-1)*n+j1-1]=0;
		}
	for(i=1;i<=n;i++)
		v[i-1]=0;		
		
	for(i=1;i<=n-1;i++)
	{
		
		//重置单位矩阵，清除旧的househorder矩阵信息 
		for(i1=1;i1<=n;i1++)
			for(j1=1;j1<=n;j1++)
			{
				P[(i1-1)*n+j1-1]=0;
				P[(i1-1)*(n+1)]=1; 
			}		
		
		//一.生成householder矩阵 
		
		
		t=0;
		for(k=i;k<=n;k++)
			t+=A[(k-1)*n+i-1]*A[(k-1)*n+i-1];
		t=sqrt(t);
		if(A[(i-1)*(n+1)]<0) t=-t;
		
	
		v[i-1]=A[(i-1)*(n+1)]+t;
		t1=v[i-1]*v[i-1];
		for(k=i+1;k<=n;k++)
		{
			v[k-1]=A[(k-1)*n+i-1];
			t1+=v[k-1]*v[k-1];
		}
		t1=2/t1;
		
		//将A转化为R ，先对第i列赋值，在对后n-i列更新 
		A[(i-1)*(n+1)]=-t;
		for(k=i+1;k<=n;k++)
			A[(k-1)*n+i-1]=0;
		for(j=i+1;j<=n;j++)
		{
			t2=0;  //对j每次循环后要将t2还原为0
			for(k=i;k<=n;k++)
				t2+=v[k-1]*A[(k-1)*n+j-1];
			t3=t1*t2;
			for(k=i;k<=n;k++)
				A[(k-1)*n+j-1]-=t3*v[k-1];
		} 
		
		//计算Q
		for(j=1;j<=n;j++)
		{
			t2=0;  //对j每次循环后要将t2还原为0
			for(k=i;k<=n;k++)
				t2+=v[k-1]*Q[(k-1)*n+j-1];
			t3=t1*t2;
			for(k=i;k<=n;k++)
				Q[(k-1)*n+j-1]-=t3*v[k-1];
		}	
		
	}
	//完成了A->R，I->Q; 
	MatT(Q,n);
	MatCon(A,R,n*n);
	MatCon(P,A,n*n);
	
}

//反幂法  输出离 p 最近的特征向量 y 
//r=fabs(m-t)
double InvPowerMethodr1(double A[],double y[],double z[],int n,int K,double epsilon,double p){
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
		y[i]=0;
		z[i]=0;
		z[0]=1;
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
	MatCon(P,A,n*n);
	//return m; 
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

//矩阵转置 
void MatT(double Q[],int n){
	double t;
	int p,q,i,j;
	for(i=1;i<=n;i++)
		for(j=i+1;j<=n;j++)
		{
			p=(i-1)*n+j-1;
			q=(j-1)*n+i-1;
			t=Q[p];
			Q[p]=Q[q];
			Q[q]=t;
		}
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

// 数组初始化 
void Initialization(double A[], int n, int m){
	int i;
	for(i=0; i<n*m; i++){ A[i]=0.; }	
}

