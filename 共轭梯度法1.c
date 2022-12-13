#include<stdio.h>
#include<math.h>

int main(){
	int CGSolver(double A[], double b[], double x[], int n, double epsilon, int K);
	int i;
	
	// 测试数据
	double A[16]={13,9,9,11,9,11,6,9,9,6,12,7,11,9,7,10},b[4]={42,35,34,37};
	double x[4]={0,0,0,0};
	double epsilon=0.000001;
	int K=2000,k,n=4;	
	
	//共轭梯度法 
	k=CGSolver(A, b, x, n, epsilon, K);
	
	// 输出
	printf("step = %d \n", k);
	printf("x =\n");
	for(i=0; i<4; i++){
		printf("%lf\n", x[i]);
	} 
	return 0;
}


int CGSolver(double A[], double b[], double x[],int n, double epsilon, int K){
	double InnerProduct(double A[],double B[],int n);
	int i,j,k;
	double err,a,d;
	double r[n],p[n],Ap[n];

	//计算 r0,p0 
	for(i=0; i<n; i++){
		for(j=0,d=0;j<n;j++){
			d+=A[i*n+j]*x[j];
		}
		r[i]=b[i]-d;
		p[i]=r[i];
	}
	
	for(k=0,err=1; err>epsilon && k<K; k++){
		//计算 x r p 更新Ap
		for(i=0; i<n; i++){
			for(j=0,d=0;j<n;j++)
				d+=A[i*n+j]*p[j];
			Ap[i]=d;	
		}
		d=InnerProduct(r,p,n)/InnerProduct(Ap,p,n);
		for(i=0;i<n;i++)
			x[i]=x[i]+d*p[i];
		for(i=0;i<n;i++)
			r[i]=r[i]-d*Ap[i];
			
		d=-InnerProduct(r,Ap,n)/InnerProduct(Ap,p,n);
		for(i=0;i<n;i++)
			p[i]=r[i]+d*p[i];
		
		//计算残差 
		err=fabs(r[0]);
		for(i=1;i<n;i++)
			if(err<fabs(r[i]))
				err=fabs(r[i]);	
	}
	return k;
}
//内积 
double InnerProduct(double A[],double B[],int n)
{
	double p=0;
	int i;
	for(i=0;i<n;i++){p+=A[i]*B[i];}
	return p;
}
