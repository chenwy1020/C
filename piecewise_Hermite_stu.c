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
�ó������ڼ���ֶ�����Hermite��ֵ 
*/ 

/*
x_node ��ʾ����������ֵ�����
y_node ��ʾ����������ֵ㴦����ֵ
dy_node ��ʾ����������ֵ㴦����ֵ

left ��ʾ������˵�
right ��ʾ�����Ҷ˵�

xs ��ʾ�����
ys ��ʾ������ֵ����ֵ

n ��ʾ�������� 
*/ 
	double *x_node, *y_node, *dy_node;
	double left = -1.0, right = 1.0;
	double xs = 1.0/3.0;
	double ys=0;
	
	int n = 80;
	int i;

//Ϊx_node y_node dy_node xs ys�����ڴ�ռ� 
	x_node = (double *)malloc(sizeof(double)*(n+1));
	y_node = (double *)malloc(sizeof(double)*(n+1));
	dy_node = (double *)malloc(sizeof(double)*(n+1));
	for(i=0;i<n+1;i++){
		x_node[i]=0;
		y_node[i]=0;
		dy_node[i]=0;
	}

//x_nodeΪ������[left, right]�ֳ�n_cells�ݺ��������ֵ� 
	uniform_partition(n, left, right, x_node);
	
//���ݾ�ȷ������������ֵ㴦�ĺ���ֵ�͵���ֵ 
	assign_function(n, x_node, y_node, dy_node);

//���зֶ�����Hermite��ֵ, ��������xs���ĺ�����ֵys 
	ys = piecewise_Hermite3(n, x_node, y_node, dy_node, xs, left, right);
	
//��������� 
	output(n, left, right, xs, ys);

//�ͷ��ڴ�ռ� 
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
	
	//��ȡ[x(j),x(j+1)]�����Ϣ 
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
�ó�������������

xs ��ʾ�����
ys ��ʾ������ֵ����ֵ

n ��ʾ�������� 

h ��ʾ���� 
*/
	double h;
	
	//���㲽�� 
	h = (right-left) / n;
	
	//������ 
	printf("�����Ϊ: %f\n", xs);
	printf("�����ֶ�Hermite��ֵΪ: %lf\n", ys);
	printf("�ֶ���Ϊ: %d\n",n);
	printf("����Ϊ: %lf\n", h);
	printf("���Ϊ: %e\n", fabs(ys-cos(xs)) );
	printf("����Ϊ: %e\n", pow(h,4)/384.0);
	
	return;
}



