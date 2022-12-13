#include<stdio.h>
#include<math.h>

double Hermite3(double *x, double *y, double *dy, double xs);

int main(){
/*
该程序用于计算三次Hermite插值 
*/
	
/*
x 表示插值点坐标
y 表示插值点函数值
dy 表示插值点导数值 
xs 表示待求点坐标
ys 表示待求点Hermite插值 
*/ 
	double x[2] = {0.0, 1.0};
	double y[2] = {0.0, 1.0};
	double dy[2] = {0.0, 3.0};
	double xs = 0.5;
	double ys = 0.0;
	
//计算Hermite插值ys 
	ys = Hermite3(x, y, dy, xs);
	
//输出计算结果 
	printf("待求点为: %lf\n", xs);
	printf("待求点三次Hermite插值为: %lf\n", ys);
	
    return 0;
	 
} 


double Hermite3(double *x, double *y, double *dy, double xs){
	double h0,dh0,h1,dh1,ys;
	
	h0=(1-2*(xs-x[0])/(x[0]-x[1]))*pow((xs-x[1])/(x[0]-x[1]),2);
	dh0=(xs-x[0])*pow((xs-x[1])/(x[0]-x[1]),2);
	h1=(1-2*(xs-x[1])/(x[1]-x[0]))*pow((xs-x[0])/(x[0]-x[1]),2); 
	dh1=(xs-x[1])*pow((xs-x[0])/(x[0]-x[1]),2);
	ys=h0*y[0]+h1*y[1]+dh0*dy[0]+dh1*dy[1];
	
	return ys;
} 































