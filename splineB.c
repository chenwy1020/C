#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#define PI 3.1415926535898

double cubic_spline(int n, double *x, double *y, double left, double right, double xs);
double cubic_spline2(int n, double *x, double *y, double left, double right, double xs);
//a b c 赋值 
void compute_diags(int n, double *x, double *h, double *a, double *b, double *c);
void compute_diags2(int n, double *x, double *h, double *a, double *b, double *c);
void compute_diags3(int n, double *x, double *h, double *a, double *b, double *c);
//构建向量 d 
void right_side(int n, double *y, double *h, double *d);
void right_side2(int n, double *y, double *h, double *d);
void right_side3(int n, double *y, double *h, double *d);
//边界条件 
void fixed_boundary(int n, double *h, double *y, double *a, double *b, double *c, double *d);
void natural_boundary(int n, double *a, double *b, double *c, double *d);
//构建 矩阵A 
void transfer_matrix(int n, double *a, double *b, double *c, double **A);
void transfer_matrix2(int n, double *a, double *b, double *c, double **A);
void transfer_matrix3(int n, double *a, double *b, double *c, double **A);
//获取待求点的插值 
double compute_value(int n, double *x, double *y, double *h, double *M, double left, double xs);
double compute_value2(int n, double *x, double *y, double *h, double *M, double left, double xs);
double compute_value3(int n, double *x, double *y, double *h, double *M, double left, double xs);
//获取函数值向量 
void assign_function(int n, double *x, double *y);
//构建等距网格 
void uniform_partition(int n, double left, double right, double *xs);
//输出 
void output(int n, double left, double right, double xs, double ys); 
//列主元的高斯消去法 
void ColElim(int n, double **A, double *b, double *x);

int main(){
/*
该程序用于计算三次样条插值 
*/ 

/*
x 表示所有子区间分点坐标
y 表示所有子区间分点处函数值

left 表示区间左端点
right 表示区间右端点

xs 表示待求点
ys 表示待求点插值函数值

n 表示子区间数 
*/ 
	double *x, *y;
	double left = -1.0, right = 1.0;
	double xs = 1.0/3.0;
	double ys;
	int i;
	
	int n = 32;
	
	//为x y分配内存空间 并初始化 
	x = (double *)malloc(sizeof(double)*(n+1));
	y = (double *)malloc(sizeof(double)*(n+1));
	for(i=0;i<=n;i++){
		x[i]=0;
		y[i]=0;
	}
	
	//对 x和 y赋值 
	uniform_partition(n, left, right, x);
	assign_function(n, x, y);
	
	//计算待求点的函数值 
	ys = cubic_spline2(n, x, y, left, right, xs);
	
	output(n, left, right, xs, ys);
	
	free(x);
	free(y);
	
	return 0;
}

//三次样条插值 
double cubic_spline(int n, double *x, double *y, double left, double right, double xs){
/*
该程序用于计算三次样条插值 
*/

/*
h 表示每个小区间长度h[j]=x[j]-x[j-1]
a 表示下次对角线
b 表示对角线
c 表示上次对角线
d 表示右端
M 表示样条函数中的系数 
*/	
	double *h, *a, *b, *c, *d;
	double *M;
	double **A;
	double ys;
	int i,j;
	
	//初始化 
	h = (double *)malloc(sizeof(double)*(n+1));
	a = (double *)malloc(sizeof(double)*(n+1));
	b = (double *)malloc(sizeof(double)*(n+1));
	c = (double *)malloc(sizeof(double)*(n+1));
	d = (double *)malloc(sizeof(double)*(n+1));
	M = (double *)malloc(sizeof(double)*(n+1));
	c[n]=0;
	A = (double **)malloc(sizeof(double *)*(n+1));
	for(i=0; i<=n; i++ ){
		A[i] = (double*)malloc(sizeof(double)*(n+1));
	}
	for(i=0; i<=n; i++){
		for(j=0; j<=n; j++){
			A[i][j] = 0.0;
		}
	}
	
	//计算每个小区间的长度 
	for(i=1; i<=n; i++){
		h[i] = x[i] - x[i-1];
	}
	
	//计算a b c 
	compute_diags(n, x, h, a, b, c);
	//计算d 
	right_side(n, y, h, d);
	
	//根据边界条件给出矩阵和右端的第一行和最后一行 
	//fixed_boundary(n, h, y, a, b, c, d);
	natural_boundary(n, a, b, c, d);
	//cycle_boundary(n, h, y, a, b, c, d);
	//natural_boundary(n, a, b, c, d);
	
	//利用a,b,c生成矩阵A
	transfer_matrix(n, a, b, c, A); 

	/*
	printf("A=\n");
	for(i=0; i<=n; i++){
		for(j=0;j<=n;j++){
			printf("%f ",A[i][j]);
		}
		printf("\n");
	}
	
	printf("d=\n");
	for(i=0;i<=n;i++){
		printf("%f\n",d[i]);
	}
	*/


	//求解AM=d, 计算出M
	ColElim(n, A, d, M);
	
	//根据求出的系数M计算xs处的样条函数值ys 
	ys = compute_value(n, x, y, h, M, left, xs);

/*
	printf("h=\n");
	for(i=0;i<=n;i++){
		printf("%f\n",h[i]);
	}
	printf("a=\n");
	for(i=0;i<=n;i++){
		printf("%f\n",a[i]);
	}
	
	printf("b=\n");
	for(i=0;i<=n;i++){
		printf("%f\n",b[i]);
	}
	
	printf("c=\n");
	for(i=0;i<=n;i++){
		printf("%f\n",c[i]);
	}
	
	printf("M=\n");
	for(i=0;i<=n;i++){
		printf("%f\n",M[i]);
	}
*/
	
	for( i=0; i<=n; i++){
		free(A[i]);
	}
	free(A);
	
	free(h);
	free(a);
	free(b);
	free(c); 
	free(d); 
	free(M);
	
	return ys;
}

//周期边界条件的三次样条插值 
double cubic_spline2(int n, double *x, double *y, double left, double right, double xs){
/*
该程序基于周期边界条件下用于计算三次样条插值  
*/

/*
h 表示每个小区间长度h[j]=x[j]-x[j-1]
a 表示下次对角线
b 表示对角线
c 表示上次对角线
d 表示右端
M 表示样条函数中的系数 
*/	
	double *h, *a, *b, *c, *d;
	double *M;
	double **A;
	double ys;
	int i,j;
	
	//分配内存并初始化 
	h = (double *)malloc(sizeof(double)*(n+1));
	a = (double *)malloc(sizeof(double)*(n));
	b = (double *)malloc(sizeof(double)*(n));
	c = (double *)malloc(sizeof(double)*(n));
	d = (double *)malloc(sizeof(double)*(n));
	M = (double *)malloc(sizeof(double)*(n));
	A = (double **)malloc(sizeof(double *)*(n));
	for(i=0; i<n; i++ ){
		A[i] = (double*)malloc(sizeof(double)*(n));
	}
	
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			A[i][j] = 0.0;
		}
		a[i]=0; b[i]=0; c[i]=0; d[i]=0; M[i]=0; 
	}
	
	//计算每个小区间的长度 
	for(i=1; i<=n; i++){
		h[i] = x[i] - x[i-1];
	}
	h[0]=h[1];
	
	//计算 a b c 
	compute_diags3(n, x, h, a, b, c);
	//计算 d 
	right_side3(n, y, h, d);
	
	//利用a,b,c生成矩阵 A
	transfer_matrix3(n, a, b, c, A); 

	/*
	printf("A=\n");
	for(i=0; i<=n-1; i++){
		for(j=0;j<=n-1;j++){
			printf("%f ",A[i][j]);
		}
		printf("\n");
	}
	
	printf("d=\n");
	for(i=0;i<=n-1;i++){
		printf("%f\n",d[i]);
	}
	*/

	//求解AM=d, 计算出M
	ColElim(n, A, d, M);

	
	//根据求出的系数M计算xs处的样条函数值ys 
	ys = compute_value3(n, x, y, h, M, left, xs);

	/*
	printf("h=\n");
	for(i=0;i<=n;i++){
		printf("%f\n",h[i]);
	}
	printf("a=\n");
	for(i=0;i<=n-1;i++){
		printf("%f\n",a[i]);
	}
	
	printf("b=\n");
	for(i=0;i<=n-1;i++){
		printf("%f\n",b[i]);
	}
	
	printf("c=\n");
	for(i=0;i<=n-1;i++){
		printf("%f\n",c[i]);
	}
	
	printf("M=\n");
	for(i=0;i<=n-1;i++){
		printf("%f\n",M[i]);
	}
	*/
	
	for( i=0; i<n; i++){
		free(A[i]);
	}
	free(A);
	
	free(h);
	free(a);
	free(b);
	free(c); 
	free(d); 
	free(M);
	
	return ys;
}

//a b c 赋值 
void compute_diags(int n, double *x, double *h, double *a, double *b, double *c){
	int i;
	for(i=1;i<n;i++){
		a[i]=h[i]/(h[i]+h[i+1]);
		b[i]=2;
		c[i]=h[i+1]/(h[i]+h[i+1]);
	}
}

void compute_diags2(int n, double *x, double *h, double *a, double *b, double *c){
	int i;
	for(i=0;i<n;i++){
		a[i]=h[i]/(h[i]+h[i+1]);
		b[i]=2;
		c[i]=h[i+1]/(h[i]+h[i+1]);
	}
}

void compute_diags3(int n, double *x, double *h, double *a, double *b, double *c){
	int i;
	for(i=0;i<n;i++){
		a[i]=1.0/6;
		b[i]=2.0/3;
		c[i]=1.0/6;
	}
}

//构建等距网格 
void uniform_partition(int n, double left, double right, double *xs) {
	/*
	该程序用于对区间进行均匀剖分, 所得剖分点包括两个端点

	n 表示剖分小区间数
	left 表示区间左端点
	right 表示区间右端点
	xs 表示剖分点坐标
	*/
	int i;
	double h;

	//h表示剖分步长
	h = (right-left) / n;

	//每个剖分点坐标即为左端点+步长*个数
	for(i=0; i<=n; i++) {
		xs[i] = left + h*i;
	}

	return;
}

//获取函数值向量 
void assign_function(int n, double *x, double *y){
	int i;
	for(i=0;i<=n;i++){
		y[i]=cos(PI*x[i]);
	} 
} 

//构建向量 d 
void right_side(int n, double *y, double *h, double *d){
	int i;
	for(i=1;i<n;i++){
		d[i]=6*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i])/(h[i]+h[i+1]);
	}
}

void right_side2(int n, double *y, double *h, double *d){
	int i;
	for(i=0;i<n;i++){
		d[i]=6*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i])/(h[i]+h[i+1]);
	}
}

void right_side3(int n, double *y, double *h, double *d){
	int i;
	for(i=0;i<n;i++){
		d[i]=y[i];
	}
}

//M0,M1=0的边界条件 
void fixed_boundary(int n, double *h, double *y, double *a, double *b, double *c, double *d){
	double dy_left= 5.0, dy_right = 5.0;
	b[0] = 2;
	c[0] = 1;
	d[0] = 6*((y[1]-y[0])/h[1]-dy_left)/h[1];
	a[n] = 1;
	b[n] = 2;
	d[n] = 6*(dy_right-(y[n]-y[n-1])/h[n])/h[n];
	//c[n] = 0;
} 

//自然边界条件 
void natural_boundary(int n, double *a, double *b, double *c, double *d){
	b[0] = 1;
	c[0] = 0;
	d[0] = 0;
	a[n] = 0;
	b[n] = 1;
	d[n] = 0;
}

//构建 矩阵A 
void transfer_matrix(int n, double *a, double *b, double *c, double **A){
	int i;
	for(i=0; i<=n; i++){
		A[i][i]  = b[i];
	}
	for(i=1; i<=n; i++){
		A[i][i-1]= a[i];
	}
	for(i=0; i<=n-1; i++){
		A[i][i+1]= c[i];
	}
} 

void transfer_matrix2(int n, double *a, double *b, double *c, double **A){
	int i;
	for(i=0; i<=n-1; i++){
		A[i][i]  = b[i];
	}
	for(i=1; i<=n-1; i++){
		A[i][i-1]= a[i];
	}
	A[0][n-1]= a[0];
	for(i=0; i<=n-1; i++){
		A[i][i+1]= c[i];
	}
	A[n-1][0]= c[n-1];
} 

void transfer_matrix3(int n, double *a, double *b, double *c, double **A){
	int i;
	for(i=0; i<=n-1; i++){
		A[i][i]  = b[i];
	}
	for(i=1; i<=n-1; i++){
		A[i][i-1]= a[i];
	}
	A[0][n-1]= a[0];
	for(i=0; i<n-1; i++){
		A[i][i+1]= c[i];
	}
	A[n-1][0]= c[n-1];
} 

//列主元的高斯消去法 
void ColElim(int n, double **A, double *b, double *x){ 
    int i,j,k,r,kk,ik,ki,ri;
    double iMax, t, Aii;
          
    for(k=0; k<n-1; k++)    { 
       	r = k;
       	kk = k*n+k;
       	iMax = fabs(A[k][k]);
       	for(i=k+1; i<n; i++){
          	t = fabs(A[i][k]);
          	if(t>iMax){
             	r = i;
             	iMax = t;                       
          	}
       	}   // 选列主元
       	if(r!=k){
          	for(i=k; i<n; i++){
             	ki = k*n+i;
             	ri = r*n+i;             
             	t = A[k][i];
             	A[k][i] = A[r][i];
             	A[r][i] = t;
          	} // 交换矩阵 A 的 r,k 两行元素
          	t = b[k];
          	b[k] = b[r];
          	b[r] = t;  // 交换 b 的 r,k 两行元素 
       	}

	   	if(fabs(A[k][k])<1e-12){
          	printf("fail\n"); 
          	return;
       	}
       	for(i=k+1; i<n; i++){
          	ik=i*n+k;
          	A[i][k] /= A[k][k];
          	b[i] -= A[i][k]*b[k];  
          	for(j=k+1; j<n; j++){ 
             	A[i][j] -= A[i][k]*A[k][j];
            } 
          	A[i][k] = 0.0;
       	}
    }
        
    kk = k*n+k;
    if(fabs(A[k][k])<1e-12){
       	printf("fail\n"); 
       	return;
    }  
    
    Aii = A[n-1][n-1];
    if(fabs(Aii)<1e-12){
       	printf("fail\n"); 
       	return;
    }
    else{ 
        x[n-1] = b[n-1];
        if(Aii!=1.0){ 
           x[n-1] /= Aii;
        } 
    }
        
    for(i = n-2; i >= 0; i--){
        Aii = A[i][i];
        if(fabs(Aii)<1e-12){
           printf("fail\n"); 
           return;
        }
        else{
           	x[i] = 0.0;   
           	for(j = i+1; j < n; j++){   
              	x[i] += A[i][j]*x[j];
            } 
           	x[i] = b[i]-x[i];
           	if(Aii!=1.0){ 
              	x[i] /= Aii;
        	}
        }
    }
    
    return;
}

//获取待求点的插值 
double compute_value(int n, double *x, double *y, double *h, double *M, double left, double xs){
	int j;
	double ys;
	//计算xs所在的小区间[x_j-1,x_j] 
	j=ceil((xs-left)/h[1]);
	
	ys=M[j-1]*pow((x[j]-xs),3)/(6*h[j])
	   	+M[j]*pow(xs-x[j-1],3)/(6*h[j])
	   	+(y[j-1]-M[j-1]*h[j]*h[j]/6.0)*(x[j]-xs)/h[j]
	   	+(y[j]-M[j]*h[j]*h[j]/6.0)*(xs-x[j-1])/h[j]; 
	   	
	return ys; 
}

double compute_value2(int n, double *x, double *y, double *h, double *M, double left, double xs){
	int j;
	double ys;
	//计算xs所在的小区间[x_j-1,x_j] 
	j=ceil((xs-left)/h[1]);
	if(j!=n){
		ys=M[j-1]*pow((x[j]-xs),3)/(6*h[j])
	   		+M[j]*pow(xs-x[j-1],3)/(6*h[j])
	   		+(y[j-1]-M[j-1]*h[j]*h[j]/6.0)*(x[j]-xs)/h[j]
	   		+(y[j]-M[j]*h[j]*h[j]/6.0)*(xs-x[j-1])/h[j]; 
	}
	if(j==n){
		ys=M[j-1]*pow((x[j]-xs),3)/(6*h[j])
	   		+M[0]*pow(xs-x[j-1],3)/(6*h[j])
	   		+(y[j-1]-M[j-1]*h[j]*h[j]/6.0)*(x[j]-xs)/h[j]
	   		+(y[j]-M[0]*h[j]*h[j]/6.0)*(xs-x[j-1])/h[j]; 
	}

	   
	return ys; 
} 

double compute_value3(int n, double *x, double *y, double *h, double *M, double left, double xs){
	int j;
	double ys,h0;
	//计算xs所在的小区间[x_j-1,x_j] 
	j=floor((xs-left)/h[1]);
	h0=h[1];
	if(j!=0){
		ys=M[j-1]*pow((x[j+1]-xs)/h0,3)/6.0
	   		+M[j]*(pow((x[j+1]+h0-xs)/h0,3)/6.0-2*pow((x[j+1]-xs)/h0,3)/3.0 )
			+M[j+1]*(pow((xs-x[j]+h0)/h0,3)/6.0-2*pow((xs-x[j])/h0,3)/3.0 )
			+M[j+2]*pow((xs-x[j])/h0,3)/6.0; 
	}
	if(j==0){
		ys=M[n-1]*pow((x[j+1]-xs)/h0,3)/6.0
	   		+M[j]*(pow((x[j+1]+h0-xs)/h0,3)/6.0-2*pow((x[j+1]-xs)/h0,3)/3.0)
			+M[j+1]*(pow((xs-x[j]+h0)/h0,3)/6.0-2*pow((xs-x[j])/h0,3)/3.0)
			+M[j+2]*pow((xs-x[j])/h0,3)/6.0; 
	}

	return ys; 
} 

void output(int n, double left, double right, double xs, double ys){
/*
该程序用于输出结果

xs 表示待求点
ys 表示待求点插值函数值

n 表示子区间分点数

h 表示步长 
*/
	double h;
	
	//计算步长 
	h = (right-left) / n;
	
	//输出结果 
	printf("分段数为: %d\n",n);
	printf("步长为: %f\n", h);
	printf("待求点为: %f\n", xs);
	printf("待求点三次样条插值为: %f\n", ys); 
	printf("误差为: %e\n", fabs(ys-cos(PI*xs)) );
	
	return;
}

