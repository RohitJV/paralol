#include<stdio.h>

float a[10][10], b[10], y[10];

void gaussian_elimination() {
    int i,j,k;
    for(i=0;i<n;i++) {
        // set the diagonal element 1, and change the elements to the right accordingly
        for(j=i+1;j<n;j++) {
            // assumption that a[i][i]!=0
            a[i][j]=a[i][j]/a[i][i];
        }
        y[i]=b[i]/a[i][i];
        a[i][i]=1;
    }
}

int main() {
    int i,j;
    printf("Enter the size of the square matrix : ");
    scanf("%d",&n);
    printf("Enter the elements of the square matrix 'A'(size nxn) in 'Ax=y', row-wise : \n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            scanf("%f",&a[i][j]);
        }
    }
    printf("Enter the elements of the vector 'y' in 'Ax=y' (size n): \n");
    for(i=0;i<n;i++) {
        scanf("%f",&b[i]);
    }

    gaussian_elimination();


}
