#include<stdio.h>
#include<math.h>

int main(){
	//函数声明 
	void MatPrintf(double Q[],int m,int n);
	void GMRES(double A[],double x0[],double b[],double x[],int m,int n,double epsilon);
	void MatMul(double A[],double B[],double y[],int m,int n,int p);
	int n=4,i,m=4;
	double A[16]={13,9,9,11,9,11,6,9,9,6,12,7,11,9,7,10},b[4]={42,35,34,37};
	double x0[4]={0,0,0,0},x[4];
	double epsilon=0.000001;
	
	GMRES(A,x0,b,x,m,n,epsilon);
	printf("\nxm = \n");
	MatPrintf(x,n,1);
}

//GMRES基本算法 
void GMRES(double A[],double x0[],double b[],double x[],int m,int n,double epsilon){
	//函数声明 
	double InnerProduct(double A[],double B[],int n);				//内积 
	double Norm2(double A[],int n);					 				//2范数
	void MatMul(double A[],double B[],double y[],int m,int n,int p);//矩阵乘法
	void Givens(double A[],double Q[],double R[],double y[],int n);	//Given变换 
	void Givensij(double A[],double c,double s,int i,int j,int n);	//对 A 的i行和 j行做的Givens变换 
	void MatConVi(double V[],double vj[],int m,int n,int j);		//继承 V的第 j列 
	void UpTriSys(double A[],double x[],double y[],int n);			//A为上三角矩阵，求解 Ax=y 
	void MatConvi(double V[],double vj[],int m,int n,int j);		//将 vj 放置在V的第 j列
	void MatPrintf(double Q[],int m,int n);							//矩阵输出 
	
	int i,j,k;
	double V[n*m],H[(m+1)*m],tempb[m+1],CS[2*m];
	double w[n],v[n],h[m+1],xm[n];
	double t,p;
	double c,s;
	
	//初始化 ，计算v1，将 e1放置在H的 第m+1列 
	MatMul(A,x0,v,n,n,1);
	for(k=0;k<n;k++){
		v[k]=b[k]-v[k];
	}
	p=Norm2(v,n);
	for(k=0;k<n;k++){
		v[k]=v[k]/p;
		V[k*m]=v[k];
	}	
	
	printf("第 0 步残差为 %lf\n",p);
	tempb[0]=p;
	
	for(j=1;j<=m;j++){
		
		/* GMRES算法 生成 H和 V */
		MatConVi(V,v,n,m,j);  //v继承 V的第 j列

		MatMul(A,v,w,n,n,1);//计算wj = Avj
		
		for(i=1;i<=j;i++){
			MatConVi(V,v,n,m,i);
			t=InnerProduct(w,v,n);//计算内积 (w,v_j) 
			H[(i-1)*m+j-1]=t; 	  //对H_ij赋值 
			for(k=1;k<=n;k++){w[k-1]-=t*v[k-1];}//更新 W		
		}
		
		t=Norm2(w,n);//计算 w的2范数 
		H[j*m+j-1]=t;//H_(j+1)j 赋值 
		

		
		/*Givens变换求残差*/
		MatConVi(H,h,m+1,m,j);//h继承 H的第 j列
		for(i=1;i<=j-1;i++){
			Givensij(h,CS[i-1],CS[m+i-1],i,i+1,1);
		}	
		
		//更新 c, s,并记录在 CS矩阵中， H_j+1,j  H_j,j 
		t=sqrt(h[j-1]*h[j-1]+h[j]*h[j]);
		c=h[j-1]/t;   s=h[j]/t; 
		CS[j-1]=c;    CS[m+j-1]=s;
		h[j-1]=t;     h[j]=0;

		//更新右端向量
		Givensij(tempb,CS[j-1],CS[m+j-1],j,j+1,1);
		
		t=H[j*m+j-1];//留存 H_j+1,j的信息 ，用于更新 V的第 j+1 列 
		
		//更新 H
		for(i=1;i<=j+1;i++){ 
			H[(i-1)*m+j-1]=h[i-1]; 
		} 
	
		//计算残差
		p=tempb[j];
		printf("第 %d 的残差为 %lf\n",j,p);
		 
		if(fabs(p)<epsilon){
			break;
		}
		if(fabs(t)<epsilon){
			break;
		}
	
		//更新 V的第 j+1 列 
		for(k=1;k<=n ;k++){ 
			V[(k-1)*m+j] = w[k-1]/t;
		}
		

	} 
	
	/*求解上三角系统*/		
	//从 H 中提出 上三角部分 ，bj
	double R[j*j],bj[j],y[j];
	for(i=1;i<=j;i++){
		for(k=1;k<=j;k++)
			R[(i-1)*j+k-1]=H[(i-1)*m+k-1];
		bj[i-1]=tempb[i-1];
	}
	UpTriSys(R,y,bj,j);
	
	//更新解
	for(i=1;i<=n;i++){
		t=0;
		for(k=1;k<=j;k++){
			t+=V[(i-1)*m+k-1]*y[k-1];
		}
		x[i-1]=x0[i-1]+t;
	} 
	
	/*
	//得到 Qm+1的转置得第一行的前 m项 
	for(i=1;i<=m;i++)
		b0[i-1]=Q[i-1]; 
	if(p<epsilon){
		UpTriSys(R,x,b0,m);
		MatMul(V,x,xm,n,m,1);
		for(i=1;i<=n;i++)xm[i]+=x0[i];
		//MatPrintf(xm,n,1);
		for(i=1;i<=n;i++){printf("%lf\n",xm[i-1]);}
	}
	*/
	
	
}
 
//Given变换
void Givens(double A[],double Q[],double R[],double y[],int n){
	//函数声明 
	void MatT(double Q[],int n);
	void MatCon(double Q[],double z[],int n);
	
	int i,j,k;
	int p,q;
	double c,s,t;
	double A0[n*n];
	MatCon(A,A0,n*n);
	
	//Q的初值 
	for(i=1;i<=n*n;i++)
		Q[i-1]=0;
	for(i=1;i<=n;i++)
		Q[(i-1)*(n+1)]=1;
	for(i=1;i<=n-1;i++) 
	{
		for(j=i+1;j<=n;j++)
		{
			//Step1: 生成Givens 矩阵（所需参数，并对Aii和Aji赋值） 
			p=(i-1)*(n+1);
			q=(j-1)*n+i-1;
			c=A[p]/sqrt(A[p]*A[p]+A[q]*A[q]);
			s=A[q]/sqrt(A[p]*A[p]+A[q]*A[q]);
			A[p]=sqrt(A[p]*A[p]+A[q]*A[q]);
			A[q]=0;
			
			//Step2: 对A第i, j 行的其他元素进行变换得到矩阵R
			for(k=i+1;k<=n;k++)
			{
				p=(i-1)*n+k-1;
				q=(j-1)*n+k-1;
				t=A[p];
				A[p]= c*t+s*A[q];
				A[q]=-s*t+c*A[q];
			}
			//Step3: 对右端项第i, j 元素进行变换
			t=y[i-1];
			y[i-1]= c*t+s*y[j-1];
			y[j-1]=-s*t+c*y[j-1];
			
			//t=e1[i-1];
			//e1[i-1]= c*t+s*e1[j-1];
			//e1[j-1]=-s*t+c*e1[j-1];
			
			//Step4: 对正交矩阵 Q 的第i, j 行元素进行变换
			for(k=1;k<=n;k++)
			{
				p=(i-1)*n+k-1;
				q=(j-1)*n+k-1;
				t=Q[p];
				Q[p]= c*t+s*Q[q];
				Q[q]=-s*t+c*Q[q];
			}	
		}
	}
	//对Q做转置 
	MatT(Q,n);	
	MatCon(A,R,n*n);
	MatCon(A0,A,n*n);
}

//对 A 的i行和 j行做的Givens变换 
void Givensij(double A[],double c,double s,int i,int j,int n){
	int k;
	double p,q;
	for(k=1;k<=n;k++){
		p=c*A[(i-1)*n+k-1]+s*A[(j-1)*n+k-1];
		q=c*A[(j-1)*n+k-1]-s*A[(i-1)*n+k-1];
		//printf("p=%lf\n",p);
		//printf("q=%lf\n",q);
		A[(i-1)*n+k-1]=p;
		A[(j-1)*n+k-1]=q;
	}
}

//Householder变换
void Householder(double A[],double Q[],double R[],double y[],int n){
	//函数声明 
	void MatPrintf(double Q[],int m,int n);
	void MatT(double Q[],int n);
	void MatCon(double Q[],double z[],int n);
	
	int i,j,k,i1,j1;
	double t,t1=0,t2=0,t3=0,v[n];
	double A0[n*n];
	MatCon(A,A0,n*n);
	//double P[n*n];
	
	//初始化 
	for(i1=1;i1<=n;i1++)
		for(j1=1;j1<=n;j1++)
		{
			Q[(i1-1)*n+j1-1]=0;
			Q[(i1-1)*(n+1)]=1; 
		}
	for(i=1;i<=n;i++)
		v[i-1]=0;		
		
	for(i=1;i<=n-1;i++)
	{
		/*
		//重置单位矩阵，清除旧的househorder矩阵信息 
		for(i1=1;i1<=n;i1++)
			for(j1=1;j1<=n;j1++)
			{
				P[(i1-1)*n+j1-1]=0;
				P[(i1-1)*(n+1)]=1; 
			}	
		*/	
		
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
			v[k-1]=A[(k-1)*n+i-1];//计算v = x + αe1 和  x + αe1的2范数 
			t1+=v[k-1]*v[k-1];
		}
		t1=2/t1;
		
		
		//将A转化为 R ，先对第i列赋值，在对后n-i列更新 
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
		
		//计算y
		t2=0;  //对j每次循环后要将t2还原为0
		for(k=i;k<=n;k++)
			t2+=v[k-1]*y[k-1];
		t3=t1*t2;
		for(k=i;k<=n;k++)
			y[k-1]-=t3*v[k-1];
		
		//计算 Q
		for(j=1;j<=n;j++)
		{
			t2=0;  //对j每次循环后要将t2还原为0
			for(k=i;k<=n;k++)
				t2+=v[k-1]*Q[(k-1)*n+j-1];
			t3=t1*t2;
			for(k=i;k<=n;k++)
				Q[(k-1)*n+j-1]-=t3*v[k-1];
		}
		MatT(Q,n); 	
		/* 
		//用单位矩阵来记录househorder矩阵的每一次信息 
		for(j=1;j<=n;j++)
		{
			t2=0;  //对j每次循环后要将t2还原为0
			for(k=i;k<=n;k++)
				t2+=v[k-1]*P[(k-1)*n+j-1];
			t3=t1*t2;
			for(k=i;k<=n;k++)
				P[(k-1)*n+j-1]-=t3*v[k-1];
		}
		printf("H(%d) = \n",i-1);
		MatPrintf(P,n,n);
		*/
	}
	//完成A->R，I->Q,y->y(new);
	MatCon(A,R,n*n);
	MatCon(A0,A,n*n);
	
}

//A为上三角矩阵，求解 Ax=y 
void UpTriSys(double A[],double x[],double y[],int n){
	int i,j,m;
	double p[n];
	for(j=0;j<n;j++)
		p[j]=y[j];
	for(j=n;j>=1;j--)
	{
		m=(j-1)*(n+1);
		x[j-1]=p[j-1]/A[m];
		for(i=j-1;i>=1;i--)
			p[i-1]=p[i-1]-A[(i-1)*n+j-1]*x[j-1];
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

//内积 
double InnerProduct(double A[],double B[],int n){
	double p=0;
	int i;
	for(i=0;i<n;i++){p+=A[i]*B[i];}
	return p;
}

//2范数 
double Norm2(double A[],int n){
	double p=0;
	int i;
	for(i=0;i<n;i++){p+=A[i]*A[i];}
	return sqrt(p);
}

//矩阵继承
void MatCon(double Q[],double z[],int n){
	int i;
	for(i=0;i<n;i++)
		z[i]=Q[i];
} 

//继承 V的第 j 列 
void MatConVi(double V[],double vj[],int m,int n,int j){
	int k;
	for(k=1;k<=m;k++){
		vj[k-1]=V[(k-1)*n+j-1];
	}
} 

//将 vj 放置在V的第 j列 
void MatConvi(double V[],double vj[],int m,int n,int j){
	int k;
	for(k=1;k<=m;k++){
		V[(k-1)*n+j-1]=vj[k-1];
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
