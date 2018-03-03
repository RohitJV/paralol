#include<stdio.h>

int n;

void mat_mult(int a[][10], int b[][10], int c[][10]) {
    int i,j,k;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            c[i][j]=0;
            for(k=0;k<n;k++) {
                c[i][j]+=(a[i][k]*b[k][j]);
            }
        }
    }
}

int main() {
    int a[10][10], b[10][10], c[10][10];
    int i,j;
    printf("Enter the size of the square matrix : ");
    scanf("%d",&n);
    printf("Enter the elements of the first square matrix (size nxn) row-wise : \n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            scanf("%d",&a[i][j]);
        }
    }
    printf("Enter the elements of the second square matrix (size nxn) row-wise : \n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            scanf("%d",&b[i][j]);
        }
    }

    mat_mult(a,b,c);

    printf("The elements of the resultant vector 'y' (size n): \n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            printf("%d ",c[i][j]);
        }
        printf("\n");
    }
}
