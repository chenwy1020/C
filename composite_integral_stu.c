#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

double f(double x);
void composite_integral(int n, double left, double right, double *result);
double trapezoid(double left, double right);
double Simpson(double left, double right);
double Gauss(double left, double right);
double Euler_Maclaurin(double left,double right,int n,int m);
double df(double x, int n);
void uniform_partition(int m, double left, double right, double *xs);

int main(){
	
	double left = 0.0, right = 1.0;
	double result[3]; 
	double EMm;
	//result表示3种复化积分公式的计算结果 
	// result[0]表示复化梯形公式,  result[1]表示复化Simpson公式, result[2]表示复化Gauss公式
	int n = 16;
	int m=2;

	composite_integral(n, left, right, result);
	EMm=Euler_Maclaurin(left,right,n,m);
	
	printf("精确值为 %f\n", sin(1.0));
	printf("复化梯形公式结果为 %f\n", result[0]);
	printf("误差为 %e\n", fabs(result[0]-sin(1.0)));
	printf("复化Simpson公式结果为 %f\n", result[1]);
	printf("误差为 %e\n", fabs(result[1]-sin(1.0)));
	printf("复化Gauss公式结果为 %f\n", result[2]);
	printf("误差为 %e\n", fabs(result[2]-sin(1.0)));
	printf("Euler-Maclaurin公式结果为 %f\n",result[0]-EMm);
	printf("误差为 %e\n", fabs(result[0]-EMm-sin(1.0)));
	
	
	return 0;
}


double f(double x){

	double y; 
	
	y = cos(x);
	
	return y;
}

void composite_integral(int n, double left, double right, double *result){
	int i;
	double x[n+1];
	double temp;
	uniform_partition(n,left,right,x);
	
	temp=0.0;
	for(i=0;i<n;i++){
		temp+=trapezoid(x[i],x[i+1]);
	}
	result[0]=temp;
	
	temp=0.0;
	for(i=0;i<n;i++){
		temp+=Simpson(x[i],x[i+1]);
	}
	result[1]=temp;
	
	temp=0.0;
	for(i=0;i<n;i++){
		temp+=Gauss(x[i],x[i+1]);
	}
	result[2]=temp;	
	
}
 
double trapezoid(double left, double right){
	double temp;
	
	temp=(right-left)*(f(left)+f(right))/2.0;
	
	return temp;
}
 
double Simpson(double left, double right){
	double temp;
	
	temp=(right-left)*(f(left)+4*f((left+right)/2.0)+f(right))/6.0;
	
	return temp;
}
 
double Gauss(double left, double right){
	double x0,x1,temp;
	
	x0=(1.0/2.0+pow(3,1.0/2.0)/6.0)*left+(1.0/2.0-pow(3,1.0/2.0)/6.0)*right;
	x1=(1.0/2.0-pow(3,1.0/2.0)/6.0)*left+(1.0/2.0+pow(3,1.0/2.0)/6.0)*right;
	
	temp=(right-left)*(f(x0)+f(x1))/2.0;
	
	return temp;
	
} 

double Euler_Maclaurin(double left,double right, int n,int m){
	int i,j,k,d;
	double h,temp;
	double B[7]={1.0, -0.5, 1.0/6.0, 0.0, -1.0/30.0, 0.0, 1.0/42.0};
	
	h=(right-left)/n;
	temp=0.0;
	
	for(i=1;i<=m;i++){
		for(j=1,k=1;j<=2*i;j++){
			k*=j; 
		}

		d=2*i-1;
		
		temp+=B[2*i]*pow(h,2*i)*(df(right,d)-df(left,d))/k;
	}
	return temp;
}

double df(double x, int n){
	if(n==1){	return -sin(x); }
	if(n==2){   return -cos(x); }
	if(n==3){   return  sin(x); }
	if(n==4){   return  cos(x); }
}

void uniform_partition(int m, double left, double right, double *xs){
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
	h = (right-left) / m;
	
	//每个剖分点坐标即为左端点+步长*个数 
	for(i=0;i<m+1;i++){
		xs[i] = left + h*i;
	}
	
}

