#include <stdio.h>

// �궨��MAX_PT��ʾ��������20����ֵ��� 
#define MAX_PT 20

double Newton(int n, double *x, double *y, double dt[MAX_PT][MAX_PT], double xs, double ys);//����Newton��ֵ 
void difference_table(int row, int col, double *x, double *y, double dt[MAX_PT][MAX_PT]);//���̱���� 
void print_table(int n, double dt[MAX_PT][MAX_PT]);//�����ֵ 

int main(){
	/*
	�ó������ڼ���Newton��ֵ, �㷨�μ��α�1.2���Լ��γ̵�4��PPT 
	
	����˵��:
	x ��ʾ��ֵ��
	y ��ʾ��ֵ�㺯��ֵ
	xs ��ʾ�����
	ys ��ʾ������Newton��ֵ
	n ��ʾ��ֵ�����
	dt ��ʾ���̱� 
	*/
	int n;
	double x[MAX_PT];
	double y[MAX_PT];
	double xs;
	double ys = 0.0;
	double dt[MAX_PT][MAX_PT];

	//��dt�ĳ�ֵ��Ϊ0 
	int i,j;
	for(i=0; i<MAX_PT; i++){
		for(j=0; j<MAX_PT; j++){
			dt[i][j] = 0.0;
		}
	}
	
	//��������xs 
	printf("����������: ");
	scanf("%lf", &xs); 
	
	//�Խ����nѭ��, ÿ������һ����ֵ���, �����ֵ���� 
	for(n=1; n<MAX_PT; n++){
		
		//���������Ĳ�ֵ������� 
		printf("\n�������%d����ֵ���: ", n);
		scanf("%lf", &x[n-1]);
		
		//���������Ĳ�ֵ��㺯��ֵ 
		printf("�������%d����ĺ���ֵ: ", n);
		scanf("%lf", &y[n-1]);
		
		/*
		//���������Ĳ�ֱ�����²��̱�
		for(i=0;i<n;i++){
			difference_table(n-1,i,x,y,dt);
		}
		*/
		 
		//���������Ĳ�ֵ������¼������㺯��ֵ 
		ys = Newton(n, x, y, dt, xs, ys);
		
		//������̱� 
		print_table(n, dt);
		//���ys 
		printf("�����Newton��ֵΪ: %f\n", ys);
		
	}
	
	return 0;
	
}

double Newton(int n, double *x, double *y, double dt[MAX_PT][MAX_PT], double xs, double ys){
	int i;
	double w;
	
	//���������Ĳ�ֱ�����²��̱�
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
	�ó�������������̱� 
	
	�������:
	n ��ʾ��ֵ������ 
	dt ��ʾ���̱� 
	*/ 
	int i,j;
	
	printf("\n��%d�����̱�Ϊ:\n",n);
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			printf("%f\t",dt[i][j]);
		}
		printf("\n");
	} 
	printf("\n");
	
	return;
}

