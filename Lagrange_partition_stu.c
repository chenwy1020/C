#include <stdio.h>
//使用malloc函数需要include标准库文件stdlib.h 
#include <stdlib.h>

void Lagrange(int n, double *x, double *y, int m, double *xs, double *ys);
void output(int n, double *x, double *y, int m, double *xs, double *ys);
void uniform_partition(int m, double left, double right, double *xs); 

int main(){
	/*
	该程序用于计算Lagrange插值, 算法参见课本1.1节以及第2次课程PPT 
	
	参数说明:
	n 表示插值点个数
	x 表示插值点
	y 表示插值点函数值
	m 表示等距剖分小区间数, m+1为待求点个数 
	xs 表示待求点
	ys 表示待求点的Lagrange插值
	*/
	int n = 3;
	double x[] = {1.0, 1.3, 1.5};
	double y[] = {1.0, 1.69, 2.25};
	int m = 5;
	double *xs;
	double *ys;

	//为xs, ys分配内存空间
	xs = (double *)malloc(sizeof(double)*(m+1));
	ys = (double *)malloc(sizeof(double)*(m+1));
	 
	//待求点坐标为将[1,2]区间等距分成5份所得的6个结点
	uniform_partition(m, 1.0, 2.0, xs); 

    //进行Lagrange插值 
	Lagrange(n, x, y, m+1, xs, ys);

	//输出结果, 将x y xs ys输出到屏幕
	output(n, x, y, m+1, xs, ys);
	
	//释放内存空间
	free(xs);
	free(ys); 

	return 0;

}

void Lagrange(int n, double *x, double *y, int m, double *xs, double *ys){
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

void output(int n, double *x, double *y, int m, double *xs, double *ys){
	/*
	该函数用于输出程序结果
	插值结点x, 插值结点函数值y, 待求点xs和待求点的Lagrange插值ys输出到屏幕
	
	输入参数:
	n 表示插值点个数
	x 表示插值点
	y 表示插值点函数值
	m 表示待求点个数
	xs 表示待求点
	ys 表示待求点的Lagrange插值
	*/
	int i;
	
	//输出插值结点 
	printf("插值结点为: "); 
	for(i=0;i<n;i++){	
		printf(" %f", x[i]);
	}
	printf("\n");
	
	//输出插值结点函数值 
	printf("插值结点函数值为: "); 
	for(i=0;i<n;i++){
		printf(" %f", y[i]);
	}
	printf("\n");

	//输出待求点 
	printf("待求点为: "); 
	for(i=0;i<m;i++){	
		printf(" %f", xs[i]);
	}
	printf("\n");
	
	//输出待求点Lagrange插值 
	printf("待求点Lagrange插值为: "); 
	for(i=0;i<m;i++){
		printf(" %f", ys[i]);
	}
	printf("\n");

	return; 
}

void uniform_partition(int m, double left, double right, double *xs){
	double h;
	int i;
	h=(right -left)/m;
	for(i=0;i<=m;i++){
		xs[i]=left+h*i;
	}
} 

 
