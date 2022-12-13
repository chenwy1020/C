#include <stdio.h>
//ʹ��malloc������Ҫinclude��׼���ļ�stdlib.h 
#include <stdlib.h>

void Lagrange(int n, double *x, double *y, int m, double *xs, double *ys);
void output(int n, double *x, double *y, int m, double *xs, double *ys);
void uniform_partition(int m, double left, double right, double *xs); 

int main(){
	/*
	�ó������ڼ���Lagrange��ֵ, �㷨�μ��α�1.1���Լ���2�ογ�PPT 
	
	����˵��:
	n ��ʾ��ֵ�����
	x ��ʾ��ֵ��
	y ��ʾ��ֵ�㺯��ֵ
	m ��ʾ�Ⱦ��ʷ�С������, m+1Ϊ�������� 
	xs ��ʾ�����
	ys ��ʾ������Lagrange��ֵ
	*/
	int n = 3;
	double x[] = {1.0, 1.3, 1.5};
	double y[] = {1.0, 1.69, 2.25};
	int m = 5;
	double *xs;
	double *ys;

	//Ϊxs, ys�����ڴ�ռ�
	xs = (double *)malloc(sizeof(double)*(m+1));
	ys = (double *)malloc(sizeof(double)*(m+1));
	 
	//���������Ϊ��[1,2]����Ⱦ�ֳ�5�����õ�6�����
	uniform_partition(m, 1.0, 2.0, xs); 

    //����Lagrange��ֵ 
	Lagrange(n, x, y, m+1, xs, ys);

	//������, ��x y xs ys�������Ļ
	output(n, x, y, m+1, xs, ys);
	
	//�ͷ��ڴ�ռ�
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
	�ú����������������
	��ֵ���x, ��ֵ��㺯��ֵy, �����xs�ʹ�����Lagrange��ֵys�������Ļ
	
	�������:
	n ��ʾ��ֵ�����
	x ��ʾ��ֵ��
	y ��ʾ��ֵ�㺯��ֵ
	m ��ʾ��������
	xs ��ʾ�����
	ys ��ʾ������Lagrange��ֵ
	*/
	int i;
	
	//�����ֵ��� 
	printf("��ֵ���Ϊ: "); 
	for(i=0;i<n;i++){	
		printf(" %f", x[i]);
	}
	printf("\n");
	
	//�����ֵ��㺯��ֵ 
	printf("��ֵ��㺯��ֵΪ: "); 
	for(i=0;i<n;i++){
		printf(" %f", y[i]);
	}
	printf("\n");

	//�������� 
	printf("�����Ϊ: "); 
	for(i=0;i<m;i++){	
		printf(" %f", xs[i]);
	}
	printf("\n");
	
	//��������Lagrange��ֵ 
	printf("�����Lagrange��ֵΪ: "); 
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

 
