#include<stdio.h>
#include<stdlib.h>

int main(){
 	void Lagrange(int n,int m,double *x,double *y,double *xs,double *ys);
	void output(int n,int m,double *x,double *y,double *xs,double *ys);
	int n=3;
	int m=2;
	double x[3]={1.2,1.3,1.4};
	double y[3]={1.44,1.69,1.96};
	double *xs;
	double *ys;
	
	xs=(double *)malloc(sizeof(double)*m);
	ys=(double *)malloc(sizeof(double)*m);
	
	xs[0]=1.5;
	xs[1]=1.6;
	
	Lagrange(n,m,x,y,xs,ys);
	
	output(n,m,x,y,xs,ys);
	
	free(xs);
	free(ys);
	return 0;
} 

void Lagrange(int n,int m,double *x,double *y,double *xs,double *ys){
	int i,j,k;
	double l;

	for(k=0;k<m;k++){
		ys[k]=0;
		for(i=0;i<n;i++){
			l=1;
			for(j=0;j<n;j++){
				if(i!=j){
					l*=(xs[k]-x[j])/(x[i]-x[j]);
				}
			}
			ys[k]+=y[i]*l;
		}	
	}
	
}

void output(int n,int m,double *x,double *y,double *xs,double *ys){
	int i;
	printf("插值结点为：");
	for(i=0;i<n;i++){
		printf("%lf,",x[i]);
	}
	printf("\n");
	
	printf("插值结点函数值为：");
	for(i=0;i<n;i++){
		printf("%lf,",y[i]);
	}
	printf("\n");

	printf("待求点为：");
	for(i=0;i<m;i++){
		printf("%lf,",xs[i]);
	}
	printf("\n");
	
	
	printf("待求点Lagrange插值为：");
	for(i=0;i<m;i++){
		printf("%lf,",ys[i]);
	}
	
	
}













