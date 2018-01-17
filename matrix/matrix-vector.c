#include<stdio.h>

int n;

void mat_vect(int a[][10], int x[], int y[]) {
    int i,j;
    for(i=0;i<n;i++) {
        y[i]=0;
        for(j=0;j<n;j++) {
            y[i]=y[i]+a[i][j]*x[j];
        }
    }
}

int main() {
    int a[10][10], x[10], y[10];
    int i,j;
    printf("Enter the size of the square matrix : ");
    scanf("%d",&n);
    printf("Enter the elements of the square matrix (size nxn) row-wise : \n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            scanf("%d",&a[i][j]);
        }
    }
    printf("Enter the elements of the vector (size n): \n");
    for(i=0;i<n;i++) {
        scanf("%d",&x[i]);
    }

    mat_vect(a,x,y);

    printf("The elements of the resultant vector 'y' (size n): \n");
    for(i=0;i<n;i++) {
        printf("%d ",y[i]);
    }
}
