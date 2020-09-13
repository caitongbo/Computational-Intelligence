/*
 * 测试PSO在CEC2014上函数F1、F3、F8、F12、F18、F20、F25、F30上的算法性能
 * 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <boost/random.hpp>

#define c1 2 //加速度因子一般是根据大量实验所得
#define c2 2
#define maxgen 500  // 迭代次数
#define sizepop 50 // 种群规模
#define popmax 100 // 个体最大取值
#define popmin -100 // 个体最小取值
#define Vmax 1 // 速度最大值
#define Vmin -1 //速度最小值
#define dim 10 // 粒子的维数
#define PI 3.1415926 //圆周率

double pop[sizepop][dim]; // 定义种群数组
double V[sizepop][dim]; // 定义种群速度数组
double fitness[sizepop]; // 定义种群的适应度数组
double result[maxgen];  //定义存放每次迭代种群最优值的数组
double pbest[sizepop][dim];  // 个体极值的位置
double gbest[dim]; //群体极值的位置
double fitnesspbest[sizepop]; //个体极值适应度的值
double fitnessgbest; // 群体极值适应度值
double genbest[maxgen][dim]; //每一代最优值取值粒子
double* f; //返回cec计算结果

double* OShift, * M, * y, * z, * x_bound;   //CEC2014外部定义
int ini_flag = 0, n_flag, func_flag, * SS;  //CEC2014外部定义

int func_num; //CEC2014功能函数

void cec14_test_func(double*, double*, int, int, int);  //CEC2014函数测试集


//适应度函数
double func(double* arr)
{
    double *x = arr; //x数组

    cec14_test_func(x, f, dim, sizepop, func_num); //调用CEC2014测试集函数

    double fitness = f[0];

    return fitness;

}

//(0,1)随机数产生函数  
double randval_x()
{
    double val;

    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_x(0, 1);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

    val = random_real_num_x();

    return(val);
}

//种群随机数产生函数  
double randval_pop()
{
    double val;

    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_x(-100, 100);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

    val = random_real_num_x();

    return(val);
}

//速度随机数产生函数  
double randval_v()
{
    double val;

    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_x(-1, 1);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_v(generator, uniform_real_generate_x);

    val = random_real_num_v();

    return(val);
}
// 种群初始化
void pop_init(void)
{
    f = (double*)malloc(sizeof(double) * sizepop);

    for (int i = 0;i < sizepop;i++)	//sizepop 种群规模
    {
        for (int j = 0;j < dim;j++)//dim 粒子的维数
        {
            pop[i][j] = randval_pop(); //-100到100之间的随机数   //pop[i][j]  种群数组
            V[i][j] = randval_v(); //-1到1之间             // V[i][j] 种群速度数组
        }
        fitness[i] = func(pop[i]); //计算适应度函数值

        //printf("/n/nfitness = %f\n",fitness[i]);

        //exit(1);s
    }
}
// max()函数定义
double* max(double* fit, int size)//寻找数组内的最大值和最大值所处的位置
{
    int index = 0; // 初始化序号
    double max = *fit; // 初始化最大值为数组第一个元素
    static double best_fit_index[2];
    for (int i = 1;i < size;i++)
    {
        if (*(fit + i) > max)
        {
            max = *(fit + i);
            index = i;
        }
    }
    best_fit_index[0] = index;
    best_fit_index[1] = max;
    return best_fit_index;

}
// 迭代寻优
void PSO_func(void)
{
    pop_init();
    double* best_fit_index; // 用于存放群体极值和其位置(序号)
    best_fit_index = max(fitness, sizepop); //求群体极值
    int index = (int)(*best_fit_index);
    // 群体极值位置
    for (int i = 0;i < dim;i++)
    {
        gbest[i] = pop[index][i];
    }
    // 个体极值位置
    for (int i = 0;i < sizepop;i++)
    {
        for (int j = 0;j < dim;j++)
        {
            pbest[i][j] = pop[i][j];
        }
    }
    // 个体极值适应度值
    for (int i = 0;i < sizepop;i++)
    {
        fitnesspbest[i] = fitness[i];
    }
    //群体极值适应度值
    double bestfitness = *(best_fit_index + 1);
    fitnessgbest = bestfitness;

    //迭代寻优
    for (int i = 0;i < maxgen;i++)
    {
        for (int j = 0;j < sizepop;j++)
        {
            //速度更新及粒子更新
            for (int k = 0;k < dim;k++)
            {
                // 速度更新
                double rand1 = randval_x(); //0到1之间的随机数
                double rand2 = randval_x(); //0到1之间的随机数
                V[j][k] = V[j][k] + c1 * rand1 * (pbest[j][k] - pop[j][k]) + c2 * rand2 * (gbest[k] - pop[j][k]);
                if (V[j][k] > Vmax)
                    V[j][k] = Vmax;
                if (V[j][k] < Vmin)
                    V[j][k] = Vmin;
                // 粒子更新
                pop[j][k] = pop[j][k] + V[j][k];
                if (pop[j][k] > popmax)
                    pop[j][k] = popmax;
                if (pop[j][k] < popmin)
                    pop[j][k] = popmin;
            }
            fitness[j] = func(pop[j]); //新粒子的适应度值
        }
        for (int j = 0;j < sizepop;j++)
        {
            // 个体极值更新
            if (fitness[j] > fitnesspbest[j])
            {
                for (int k = 0;k < dim;k++)
                {
                    pbest[j][k] = pop[j][k];
                }
                fitnesspbest[j] = fitness[j];
            }
            // 群体极值更新
            if (fitness[j] > fitnessgbest)
            {
                for (int k = 0;k < dim;k++)
                    gbest[k] = pop[j][k];
                fitnessgbest = fitness[j];
            }
        }
        for (int k = 0;k < dim;k++)
        {
            genbest[i][k] = gbest[k]; // 每一代最优值取值粒子位置记录
        }
        result[i] = fitnessgbest; // 每代的最优值记录到数组
    }
}

//功能函数
void start_Function(int n) {

    func_num = n;
    
    printf("\nF%d开始计算，请稍后...\n",func_num);

    clock_t start, finish; //程序开始和结束时间
    start = clock(); //开始计时
    srand((unsigned)time(NULL)); // 初始化随机数种子
    PSO_func();

    double* best_arr;
    best_arr = max(result, maxgen);   //result  存放每次迭代种群最优值的数组   //maxgen 迭代次数
    int best_gen_number = *best_arr; // 最优值所处的代数
    double best = *(best_arr + 1); //最优值

    printf("迭代了%d次，在第%d次取到最优值，最优值为:%lf.\n", maxgen, best_gen_number + 1, best);
    printf("取到最优值的位置为(%lf,%lf).\n", genbest[best_gen_number][0], genbest[best_gen_number][1]);
    finish = clock(); //结束时间
    double duration = (double)(finish - start) / CLOCKS_PER_SEC; // 程序运行时间
    printf("程序运行耗时:%lf\n", duration);
    printf("F%d计算完成。\n", func_num);


}
// 主函数
int main(void)
{
    //开启测试函数，如不进行某测试函数，通过"//"注释掉即可
    start_Function(1); //Rotated High Conditioned Elliptic Funcktion

    start_Function(3); //Rotated Discus Function

    start_Function(8); //Shifted Rastrigin's Function

    start_Function(12); //Shiftedand Rotated Katsuura Function

    start_Function(18); //Hybrid Function 2 (N=3) 

    start_Function(20); //Hybrid Function 4 (N=4) 

    start_Function(25); //Composition Function 3 (N=3) 

    start_Function(30); //Composition Function 8 (N=3) 

}