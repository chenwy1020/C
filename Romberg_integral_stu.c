#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

double f(double x);
double composite_integral(int n, double left, double right);
double trapezoid(double left, double right);
void uniform_partition(int cells, double left, double right, double *xs);
void Romberg(int m, int n0, double left, double right, double *int_f);

int main(){
	
	double left = 0.0, right = 1.0;
	double *int_f;
	int m = 4;
	int n0 = 4;
	int i;
	
	int_f = (double *)malloc(sizeof(double)*m);
	
	Romberg(m, n0, left, right, int_f);
	
	printf("精确值为 %.16f\n", sin(1.0));
	
	for(i=0;i<m;i++){
		printf("\n小区间数为 %d\n", n0*(int)(pow(2,i)));
		printf("Romberg积分法结果为 %.16f\n", int_f[i]);
		printf("误差为 %e\n", fabs(int_f[i]-sin(1.0)));
	}
	
	free(int_f);
	
	return 0;
}

double f(double x){

	double y; 
	
	y = cos(x);
	
	return y;
}

double composite_integral(int n, double left, double right){
	int i;
	double x[n+1];
	double temp;
	uniform_partition(n,left,right,x);
	
	temp=0.0;
	for(i=0;i<n;i++){
		temp+=trapezoid(x[i],x[i+1]);
	}
	return temp;
	
}
 
double trapezoid(double left, double right){
	double temp;
	
	temp=(right-left)*(f(left)+f(right))/2.0;
	
	return temp;
}
 
void uniform_partition(int cells, double left, double right, double *xs){
	/*
	该程序用于对区间进行均匀剖分, 所得剖分点包括两个端点 
	
	m 表示剖分点个数 
	left 表示区间左端点
	right 表示区间右端点
	xs 表示剖分点坐标 
	*/	
	int i;
	double h;  
	
	//h表示剖分步长 
	h = (right-left) / cells;
	
	//每个剖分点坐标即为左端点+步长*个数 
	for(i=0;i<cells+1;i++){
		xs[i] = left + h*i;
	}
	
	return;
}

void Romberg(int m, int n0, double left, double right, double *int_f){
	
	double **int_table;
	int i,j,n;
	
	int_table = (double **)malloc(sizeof(double *)*m);
	for(i=0; i<m; i++ ){
		int_table[i]= (double*)malloc(sizeof(double)*m);
	}
	for(i=0;i<m;i++){
		for(j=0;j<m;j++){
			int_table[i][j]=0.0;
		}
	}
	
	for(i=0;i<m;i++){
		n=n0*pow(2,i);
		int_table[i][0]=composite_integral(n,left,right);
	}
	for(j=1;j<m;j++){
		for(i=j;i<m;i++){
			int_table[i][j]=(pow(4,j)*int_table[i][j-1]-int_table[i-1][j-1])/(pow(4,j)-1);
		}
	}
	for(i=0;i<m;i++){
		int_f[i]=int_table[i][i];
	}

}









