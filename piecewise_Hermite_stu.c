#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double piecewise_Hermite3(int n, double *x_node, double *y_node, double *dy_node, double xs, double left, double right);
void uniform_partition(int cells, double left, double right, double *xs);
void assign_function(int n, double *x_node, double *y_node, double *dy_node);
double Hermite3(double *x, double *y, double *dy, double xs);
void output(int n, double left, double right, double xs, double ys);

int main(){
/*
该程序用于计算分段三次Hermite插值 
*/ 

/*
x_node 表示所有子区间分点坐标
y_node 表示所有子区间分点处函数值
dy_node 表示所有子区间分点处导数值

left 表示区间左端点
right 表示区间右端点

xs 表示待求点
ys 表示待求点插值函数值

n 表示子区间数 
*/ 
	double *x_node, *y_node, *dy_node;
	double left = -1.0, right = 1.0;
	double xs = 1.0/3.0;
	double ys=0;
	
	int n = 80;
	int i;

//为x_node y_node dy_node xs ys分配内存空间 
	x_node = (double *)malloc(sizeof(double)*(n+1));
	y_node = (double *)malloc(sizeof(double)*(n+1));
	dy_node = (double *)malloc(sizeof(double)*(n+1));
	for(i=0;i<n+1;i++){
		x_node[i]=0;
		y_node[i]=0;
		dy_node[i]=0;
	}

//x_node为将区间[left, right]分成n_cells份后的子区间分点 
	uniform_partition(n, left, right, x_node);
	
//根据精确函数赋子区间分点处的函数值和导数值 
	assign_function(n, x_node, y_node, dy_node);

//进行分段三次Hermite插值, 计算待求点xs处的函数插值ys 
	ys = piecewise_Hermite3(n, x_node, y_node, dy_node, xs, left, right);
	
//输出计算结果 
	output(n, left, right, xs, ys);

//释放内存空间 
	free(x_node);
	free(y_node);
	free(dy_node);

	return 0;
}

double piecewise_Hermite3(int n, double *x_node, double *y_node, double *dy_node, double xs, double left, double right){
	double h,ys=0;
	double x[2],y[2],dy[2];
	int j;
	
	
	h=(right-left)/n;
	j=floor((xs-left)/h);
	
	//获取[x(j),x(j+1)]相关信息 
	x[0]=x_node[j];
	x[1]=x_node[j+1];
	y[0]=y_node[j];
	y[1]=y_node[j+1];
	dy[0]=dy_node[j];
	dy[1]=dy_node[j+1];
	 
	ys = Hermite3(x,y,dy,xs) ;
}

void assign_function(int n, double *x_node, double *y_node, double *dy_node){
	int i;
	for(i=0;i<n+1;i++){
		y_node[i] = cos(x_node[i]);
		dy_node[i] = -sin(x_node[i]);
		
	}
}

void uniform_partition(int cells, double left, double right, double *xs){
	int i;
	double h;
	h=(right-left)/cells;
	for(i=0;i<cells+1;i++){
		xs[i]=left+h*i;
	}
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

void output(int n, double left, double right, double xs, double ys){
/*
该程序用于输出结果

xs 表示待求点
ys 表示待求点插值函数值

n 表示子区间数 

h 表示步长 
*/
	double h;
	
	//计算步长 
	h = (right-left) / n;
	
	//输出结果 
	printf("待求点为: %f\n", xs);
	printf("待求点分段Hermite插值为: %lf\n", ys);
	printf("分段数为: %d\n",n);
	printf("步长为: %lf\n", h);
	printf("误差为: %e\n", fabs(ys-cos(xs)) );
	printf("误差界为: %e\n", pow(h,4)/384.0);
	
	return;
}



