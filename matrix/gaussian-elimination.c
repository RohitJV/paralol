#include<stdio.h>

float a[10][10], b[10], y[10];
int n;


// implementation according to the textbook
/*
    Ax=y -> Ux=b, where U is a upper triangular matrix
    In the below code, we replace elements in A with elements of U
*/
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

        /*
        The below invariant establishes that the matrix is indeed a Upper Triangular matrix:
        For each row below, apply operations such that the element a[.][i]=0
        */
        for(j=i+1;j<n;j++) {
            for(k=i+1;k<n;k++) {
                a[j][k] = a[j][k] - a[i][k]*a[j][i];
            }
            a[j][i]=0;
            b[j]=b[j]-a[j][i]*y[i];
        }
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

    printf("Elements of Upper Triangular matrix 'U' in 'Ux=y' : \n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++)
            printf("%f ",a[i][j]);
        printf("\n");
    }
    printf("Elements of 'y' in 'Ux=y' : \n");
    for(i=0;i<n;i++) {
        printf("%f ",y[i]);
    }
}
