#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#define PI 3.1415926535898

double cubic_spline(int n, double *x, double *y, double left, double right, double xs);
double cubic_spline2(int n, double *x, double *y, double left, double right, double xs);
//a b c ��ֵ 
void compute_diags(int n, double *x, double *h, double *a, double *b, double *c);
void compute_diags2(int n, double *x, double *h, double *a, double *b, double *c);
void compute_diags3(int n, double *x, double *h, double *a, double *b, double *c);
//�������� d 
void right_side(int n, double *y, double *h, double *d);
void right_side2(int n, double *y, double *h, double *d);
void right_side3(int n, double *y, double *h, double *d);
//�߽����� 
void fixed_boundary(int n, double *h, double *y, double *a, double *b, double *c, double *d);
void natural_boundary(int n, double *a, double *b, double *c, double *d);
//���� ����A 
void transfer_matrix(int n, double *a, double *b, double *c, double **A);
void transfer_matrix2(int n, double *a, double *b, double *c, double **A);
void transfer_matrix3(int n, double *a, double *b, double *c, double **A);
//��ȡ�����Ĳ�ֵ 
double compute_value(int n, double *x, double *y, double *h, double *M, double left, double xs);
double compute_value2(int n, double *x, double *y, double *h, double *M, double left, double xs);
double compute_value3(int n, double *x, double *y, double *h, double *M, double left, double xs);
//��ȡ����ֵ���� 
void assign_function(int n, double *x, double *y);
//�����Ⱦ����� 
void uniform_partition(int n, double left, double right, double *xs);
//��� 
void output(int n, double left, double right, double xs, double ys); 
//����Ԫ�ĸ�˹��ȥ�� 
void ColElim(int n, double **A, double *b, double *x);

int main(){
/*
�ó������ڼ�������������ֵ 
*/ 

/*
x ��ʾ����������ֵ�����
y ��ʾ����������ֵ㴦����ֵ

left ��ʾ������˵�
right ��ʾ�����Ҷ˵�

xs ��ʾ�����
ys ��ʾ������ֵ����ֵ

n ��ʾ�������� 
*/ 
	double *x, *y;
	double left = -1.0, right = 1.0;
	double xs = 1.0/3.0;
	double ys;
	int i;
	
	int n = 32;
	
	//Ϊx y�����ڴ�ռ� ����ʼ�� 
	x = (double *)malloc(sizeof(double)*(n+1));
	y = (double *)malloc(sizeof(double)*(n+1));
	for(i=0;i<=n;i++){
		x[i]=0;
		y[i]=0;
	}
	
	//�� x�� y��ֵ 
	uniform_partition(n, left, right, x);
	assign_function(n, x, y);
	
	//��������ĺ���ֵ 
	ys = cubic_spline2(n, x, y, left, right, xs);
	
	output(n, left, right, xs, ys);
	
	free(x);
	free(y);
	
	return 0;
}

//����������ֵ 
double cubic_spline(int n, double *x, double *y, double left, double right, double xs){
/*
�ó������ڼ�������������ֵ 
*/

/*
h ��ʾÿ��С���䳤��h[j]=x[j]-x[j-1]
a ��ʾ�´ζԽ���
b ��ʾ�Խ���
c ��ʾ�ϴζԽ���
d ��ʾ�Ҷ�
M ��ʾ���������е�ϵ�� 
*/	
	double *h, *a, *b, *c, *d;
	double *M;
	double **A;
	double ys;
	int i,j;
	
	//��ʼ�� 
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
	
	//����ÿ��С����ĳ��� 
	for(i=1; i<=n; i++){
		h[i] = x[i] - x[i-1];
	}
	
	//����a b c 
	compute_diags(n, x, h, a, b, c);
	//����d 
	right_side(n, y, h, d);
	
	//���ݱ߽���������������Ҷ˵ĵ�һ�к����һ�� 
	//fixed_boundary(n, h, y, a, b, c, d);
	natural_boundary(n, a, b, c, d);
	//cycle_boundary(n, h, y, a, b, c, d);
	//natural_boundary(n, a, b, c, d);
	
	//����a,b,c���ɾ���A
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


	//���AM=d, �����M
	ColElim(n, A, d, M);
	
	//���������ϵ��M����xs������������ֵys 
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

//���ڱ߽�����������������ֵ 
double cubic_spline2(int n, double *x, double *y, double left, double right, double xs){
/*
�ó���������ڱ߽����������ڼ�������������ֵ  
*/

/*
h ��ʾÿ��С���䳤��h[j]=x[j]-x[j-1]
a ��ʾ�´ζԽ���
b ��ʾ�Խ���
c ��ʾ�ϴζԽ���
d ��ʾ�Ҷ�
M ��ʾ���������е�ϵ�� 
*/	
	double *h, *a, *b, *c, *d;
	double *M;
	double **A;
	double ys;
	int i,j;
	
	//�����ڴ沢��ʼ�� 
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
	
	//����ÿ��С����ĳ��� 
	for(i=1; i<=n; i++){
		h[i] = x[i] - x[i-1];
	}
	h[0]=h[1];
	
	//���� a b c 
	compute_diags3(n, x, h, a, b, c);
	//���� d 
	right_side3(n, y, h, d);
	
	//����a,b,c���ɾ��� A
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

	//���AM=d, �����M
	ColElim(n, A, d, M);

	
	//���������ϵ��M����xs������������ֵys 
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

//a b c ��ֵ 
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

//�����Ⱦ����� 
void uniform_partition(int n, double left, double right, double *xs) {
	/*
	�ó������ڶ�������о����ʷ�, �����ʷֵ���������˵�

	n ��ʾ�ʷ�С������
	left ��ʾ������˵�
	right ��ʾ�����Ҷ˵�
	xs ��ʾ�ʷֵ�����
	*/
	int i;
	double h;

	//h��ʾ�ʷֲ���
	h = (right-left) / n;

	//ÿ���ʷֵ����꼴Ϊ��˵�+����*����
	for(i=0; i<=n; i++) {
		xs[i] = left + h*i;
	}

	return;
}

//��ȡ����ֵ���� 
void assign_function(int n, double *x, double *y){
	int i;
	for(i=0;i<=n;i++){
		y[i]=cos(PI*x[i]);
	} 
} 

//�������� d 
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

//M0,M1=0�ı߽����� 
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

//��Ȼ�߽����� 
void natural_boundary(int n, double *a, double *b, double *c, double *d){
	b[0] = 1;
	c[0] = 0;
	d[0] = 0;
	a[n] = 0;
	b[n] = 1;
	d[n] = 0;
}

//���� ����A 
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

//����Ԫ�ĸ�˹��ȥ�� 
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
       	}   // ѡ����Ԫ
       	if(r!=k){
          	for(i=k; i<n; i++){
             	ki = k*n+i;
             	ri = r*n+i;             
             	t = A[k][i];
             	A[k][i] = A[r][i];
             	A[r][i] = t;
          	} // �������� A �� r,k ����Ԫ��
          	t = b[k];
          	b[k] = b[r];
          	b[r] = t;  // ���� b �� r,k ����Ԫ�� 
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

//��ȡ�����Ĳ�ֵ 
double compute_value(int n, double *x, double *y, double *h, double *M, double left, double xs){
	int j;
	double ys;
	//����xs���ڵ�С����[x_j-1,x_j] 
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
	//����xs���ڵ�С����[x_j-1,x_j] 
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
	//����xs���ڵ�С����[x_j-1,x_j] 
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
�ó�������������

xs ��ʾ�����
ys ��ʾ������ֵ����ֵ

n ��ʾ������ֵ���

h ��ʾ���� 
*/
	double h;
	
	//���㲽�� 
	h = (right-left) / n;
	
	//������ 
	printf("�ֶ���Ϊ: %d\n",n);
	printf("����Ϊ: %f\n", h);
	printf("�����Ϊ: %f\n", xs);
	printf("���������������ֵΪ: %f\n", ys); 
	printf("���Ϊ: %e\n", fabs(ys-cos(PI*xs)) );
	
	return;
}

