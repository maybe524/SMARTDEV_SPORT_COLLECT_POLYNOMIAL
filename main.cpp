/*
* 本实验根据数组x[], y[]列出的一组数据，用最小二乘法求它的拟合曲线。
* 近似解析表达式为y = a0 + a1 * x + a2 * x^2 + a3 * x^3;
*/
#include <stdio.h>
#include <math.h>

#define DEBUG
#define MAXN    12
#define RANK    3

int main()
{
    int i, j, k;
    double x[MAXN] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
    double y[MAXN] = {0, 1.27, 2.16, 2.86, 3.44, 3.87, 4.15, 4.37, 4.51, 4.58, 4.02, 4.64};
    double atemp[2 * (RANK + 1)] = {0}, b[RANK + 1] = {0}, a[RANK + 1][RANK + 1];

    for (i = 0; i < MAXN; i++) {  //
        atemp[1] += x[i];
        atemp[2] += pow(x[i], 2);
        atemp[3] += pow(x[i], 3);
        atemp[4] += pow(x[i], 4);
        atemp[5] += pow(x[i], 5);
        atemp[6] += pow(x[i], 6);
        b[0] += y[i];
        b[1] += x[i] * y[i];
        b[2] += pow(x[i], 2) * y[i];
        b[3] += pow(x[i], 3) * y[i];
    }

    atemp[0] = MAXN;

#ifdef DEBUG
    for(i = 0; i <= 2 * RANK; i++)  printf("atemp[%d] = %f\n", i, atemp[i]);
    printf("\n");
    for(i = 0; i <= RANK; i++)  printf("b[%d] = %f\n", i, b[i]);
    printf("\n");
#endif

    for (i = 0; i < RANK + 1; i++) {  // 构建线性方程组系数矩阵，b[]不变
        k = i;
        for (j = 0; j < RANK + 1; j++) a[i][j] = atemp[k++];
    }

#ifdef DEBUG
    for(i = 0; i < RANK + 1; i++){
        for(j = 0; j < RANK + 1; j++)  printf("a[%d][%d] = %-17f  ", i, j, a[i][j]);
        printf("\n");
    }
    printf("\n");
#endif

    // 以下为高斯列主元消去法解线性方程组
    for (k = 0; k < RANK + 1 - 1; k++) {  // n - 1列
        int column = k;
        double mainelement = a[k][k];

        for (i = k; i < RANK + 1; i++)  // 找主元素
            if (fabs(a[i][k]) > mainelement) {
                mainelement = fabs(a[i][k]);
                column = i;
            }
        for (j = k; j < RANK + 1; j++) {  // 交换两行
            double atemp = a[k][j];
            a[k][j] = a[column][j];
            a[column][j] = atemp;
        }
        double btemp = b[k];
        b[k] = b[column];
        b[column] = btemp;

        for (i = k + 1; i < RANK + 1; i++) {  // 消元过程
            double Mik = a[i][k] / a[k][k];
            for (j = k; j < RANK + 1; j++) a[i][j] -= Mik * a[k][j];
            b[i] -= Mik * b[k];
        }
    }

#ifdef DEBUG
    for(i = 0; i < RANK + 1; i++) {  // 经列主元高斯消去法得到的上三角阵(最后一列为常系数)
        for(j = 0; j < RANK + 1; j++)  printf("%20f", a[i][j]);
        printf("%20f\n", b[i]);
    }
    printf("\n");
#endif

    b[RANK + 1 - 1] /= a[RANK + 1 - 1][RANK + 1 - 1];  // 回代过程
    for (i = RANK + 1 - 2; i >= 0; i--) {
        double sum = 0;
        for (j = i + 1; j < RANK + 1; j++) sum += a[i][j] * b[j];
        b[i] = (b[i] - sum) / a[i][i];
    }
    // 高斯列主元消去法结束

    printf("P(x) = %f%+fx%+fx^2%+fx^3\n\n", b[0], b[1], b[2], b[3]);

#ifdef DEBUG
    for(i = 0; i < MAXN; i++) {  // 误差比较
        double temp = b[0] + b[1] * x[i] + b[2] * x[i] * x[i] + b[3] * x[i] * x[i] * x[i];
        printf("%f    %f    error: %f\n", y[i], temp, temp - y[i]);
    }
#endif

    return 0;
}