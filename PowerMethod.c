#include<stdio.h>
#include<math.h>

int main(){
	double PowerMethod(double A[],double y[],double z[],int n,int K,double epsilon);
	void MatPrintf(double Q[],int m,int n);
	
	int n,i,K;
	
	scanf("%d",&n);
	//printf("%d\n",n);
	
	double A[n*n],z[n],y[n];
	//double A[9]={2,-1,0,0,2,-1,0,-1,2},z[3]={0,0,1},y[3];
	double epsilon,m;
	//n=3; K=10000;  epsilon=1e-8;
	
	//输入 
	for(i=0;i<n*n;i++){
		scanf("%lf",A+i);
	}
	//MatPrintf(A,n,n);
	for(i=0;i<n;i++){
		scanf("%lf",z+i);
	}
	//MatPrintf(z,1,n); 
	scanf("%lf",&epsilon);
	//printf("%.8lf\n",epsilon);
	scanf("%d",&K);
	//printf("%d\n",K);

	//输出 
	m=PowerMethod(A,y,z,n,K,epsilon);
	printf("%.4lf\n",m);
	for(i=0;i<3;i++){
		printf("%.4lf ",y[i]);
	}
	
}

//乘幂法
double PowerMethod(double A[],double y[],double z[],int n,int K,double epsilon){
	void MatMul(double A[],double B[],double y[],int m,int n,int p);//矩阵乘法
	int Vecmax(double v[],int n);
	
	int i,k,j;
	double m,r,t;
	
	MatMul(A,z,y,n,n,1);//y=Az
	j=Vecmax(y,n);
	m=y[j];
	t=m;
	for(k=0;k<n;k++){
		z[k]=y[k]/m;
	}
	
	for(k=0;k<K;k++){
		MatMul(A,z,y,n,n,1);//更新 y 
		j=Vecmax(y,n);
		m=y[j];//绝对值最大的分量 
		r=fabs(m-t);
		if(r<epsilon) break;
		t=m;//更新t 
		//更新 z 
		for(i=0;i<n;i++){
			z[i]=y[i]/m;
		}
	}
	
	return m; 
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





