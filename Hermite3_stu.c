#include<stdio.h>
#include<math.h>

double Hermite3(double *x, double *y, double *dy, double xs);

int main(){
/*
�ó������ڼ�������Hermite��ֵ 
*/
	
/*
x ��ʾ��ֵ������
y ��ʾ��ֵ�㺯��ֵ
dy ��ʾ��ֵ�㵼��ֵ 
xs ��ʾ���������
ys ��ʾ�����Hermite��ֵ 
*/ 
	double x[2] = {0.0, 1.0};
	double y[2] = {0.0, 1.0};
	double dy[2] = {0.0, 3.0};
	double xs = 0.5;
	double ys = 0.0;
	
//����Hermite��ֵys 
	ys = Hermite3(x, y, dy, xs);
	
//��������� 
	printf("�����Ϊ: %lf\n", xs);
	printf("���������Hermite��ֵΪ: %lf\n", ys);
	
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































