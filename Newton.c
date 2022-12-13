#include<time.h>
#include<stdio.h>

// 宏定义MAX_PT表示最多可允许20个插值结点 
#define MAX_PT 20

double Newton(int n, double *x, double *y, double dt[MAX_PT][MAX_PT], double xs, double ys);//计算Newton差值 
void difference_table(int row, int col, double *x, double *y, double dt[MAX_PT][MAX_PT]);//差商表计算 
void print_table(int n, double dt[MAX_PT][MAX_PT]);//输出差值 

int main(){
	/*
	该程序用于计算Newton插值, 算法参见课本1.2节以及课程第4次PPT 
	
	参数说明:
	x 表示插值点
	y 表示插值点函数值
	xs 表示待求点
	ys 表示待求点的Newton插值
	n 表示插值点个数
	dt 表示差商表 
	*/
	int n,m=6;
	double x[MAX_PT];
	double y[MAX_PT];
	double xs[m];
	double ys[m];
	double dt[MAX_PT][MAX_PT];
	double t;//步长
	t=1.0/(m-1); 

	//将dt的初值赋为0 
	int i,j;
	for(i=0; i<MAX_PT; i++){
		for(j=0; j<MAX_PT; j++){
			dt[i][j] = 0.0;
		}
	}
	
	for(i=0;i<m;i++){
		xs[i]=t*i;
	}
	

	
	//对结点数n循环, 每次新增一个插值结点, 计算插值函数 
	for(n=1; n<MAX_PT; n++){
		
		//输入新增的插值结点坐标 
		printf("请输入第%d个插值结点: ", n);
		scanf("%lf", &x[n-1]);
		
		//输入新增的插值结点函数值 
		printf("请输入第%d个点的函数值: ", n);
		scanf("%lf", &y[n-1]);
	
	
	
	//	for(i=0;i<n;i++){
	//		difference_table(n-1,i,x,y,dt);
	//	}
		
		

		//根据新增的插值结点重新计算待求点函数值 
		for(j=0;j<m;j++){
			ys[j] = Newton(n, x, y, dt, xs[j], ys[j]);
		}
	
		
		//输出差商表 
		print_table(n, dt);
		
			//输入待求点xs 
		printf("待求点为: \t\t");
		for(i=0;i<m;i++){
			printf("%lf  ",xs[i]);
		}
		printf("\n");
		
		//输出ys 
		printf("待求点Newton插值为: \t");
		for(i=0;i<m;i++){
			printf("%lf  ", ys[i]);
		}
		printf("\n");
		printf("\n");
		
	}
	
	return 0;
	
}

double Newton(int n, double *x, double *y, double dt[MAX_PT][MAX_PT], double xs, double ys){
	int i;
	double w;
	
	//更新差商表
	for(i=0;i<n;i++){
		difference_table(n-1,i,x,y,dt);
	}
	
	w=1.0;
	//
	for(i=0;i<n-1;i++){
		w=w*(xs-x[i]);
	}
	ys=ys+w*dt[n-1][n-1];
	
	return ys;
} 

void difference_table(int row, int col, double *x, double *y, double dt[MAX_PT][MAX_PT]){
	
	
	
	if(col==0){
		dt[row][col]=y[row];
	} 
	else{
		dt[row][col]=(dt[row][col-1]-dt[row-1][col-1])/(x[row]-x[row-col]);
	}
} 

void print_table(int n, double dt[MAX_PT][MAX_PT]){
	/*
	该程序用于输出差商表 
	
	输入参数:
	n 表示插值结点个数 
	dt 表示差商表 
	*/ 
	int i,j;
	
	printf("\n第%d步差商表为:\n",n);
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			printf("%f\t",dt[i][j]);
		}
		printf("\n");
	} 
	printf("\n");
	
	return;
}

