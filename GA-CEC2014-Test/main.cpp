/*
    测试GA在CEC2014上函数F1、F3、F8、F12、F18、F20、F25、F30上的算法性能
*/

#include <windows.h>
#include <time.h>
#include <stdio.h>  
#include <stdlib.h>  
#include <math.h> 
#include <malloc.h>
#include <boost/random.hpp>
#include <fstream>
#include <string>
#include<iostream>

using namespace std;

//请根据你的需要来修改以下参数   
#define MAXGENS 20000 //最大迭代次数，每2000代记录一次 [0,20000] 
const int NVARS = 10; //问题变量个数 [2,10,20,30,50,100]

#define TRUE 1 
#define FALSE 0    
#define PI 3.1415926

struct genotype //种群的一个基因的结构体类型
{
    double gene[NVARS]; //变量
    double fitness; //基因的适应度
    double upper[NVARS]; //基因变量的上界
    double lower[NVARS]; //基因变量的下界
    double rfitness; //比较适应度
    double cfitness; //积累适应度
};

struct genotype* population; //种群
struct genotype* newpopulation; //新种群取代旧的基因

//以下是一些函数声明 
void initialize(void);  //初始化随机变量
double randval(double, double);    //随机数产生函数
void evaluate(void);             //评价函数，可以由用户自定义，该函数取得每个基因的适应度
void keep_the_best(void);       //保存每次遗传后的最佳基因
void elitist(void);  //搜寻杰出个体函数：找出最好和最坏的个体。如果某代的最好个体比前一代的最好个体要坏，那么后者将会取代当前种群的最坏个体 
void select(void);  //选择函数：用于最大化合并杰出模型的标准比例选择，保证最优秀的个体得以生存
void crossover(void);  //杂交函数：选择两个个体来杂交，这里用单点杂交
void Xover(int, int);   //交叉 
void swap(double*, double*);  //交换
void mutate(void);    //变异函数：被该函数选中后会使得某一变量被一个随机的值所取代 
void report(void);   //报告模拟进展情况
void start_Function(int);   //开启功能函数
void test_popsize(void); //测试种群大小对算法的性能影响 [20,200]
void test_pxover(void);  //测试交叉大小对算法的性能影响 [0.9,0.1]
void test_pmutation(void);   //测试变异大小对算法的性能影响 [0.01,0.1]
void cec14_test_func(double*, double*, int, int, int);  //CEC2014函数测试集

double* OShift, * M, * y, * z, * x_bound;   //CEC2014外部定义
double* f, * x; //f是适应度大小，x是未知量的值
double PXOVER; //交叉概率 [0.9,0.1] 步长-0.1
double PMUTATION; //变异概率 [0.01,0.1] 步长+0.01
int POPSIZE; //种群大小 [20,200] 步长+20
int ini_flag = 0, n_flag, func_flag, * SS;  //CEC2014外部定义
int func_num; //标记功能函数，1、3、8、12、18、20、25、30
int generation; //当前基因个数
int cur_best; //最优个体
FILE* txtFile;  //结果输出文件

double lbound = -100, ubound = 100; //定义随机数上下界 (-100,100)

//全局随机数生成器，用于变异、交叉等操作
boost::mt19937 generator(time(0)* rand());
boost::uniform_real<> uniform_real_generate_r(0, 1);
boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r(generator, uniform_real_generate_r);


//随机数初始化
void initialize(void)
{
    int i, j;

    x = (double*)malloc(NVARS * POPSIZE * sizeof(double));
    f = (double*)malloc(sizeof(double) * POPSIZE);

    population = new genotype[POPSIZE + 1]; //种群
    newpopulation = new genotype[POPSIZE + 1]; //新种群,取代旧的基因

    for (i = 0; i < NVARS; i++)
    {
        for (j = 0; j < POPSIZE; j++)
        {
            population[j].fitness = 0;     //基因的适应度
            population[j].rfitness = 0;    //比较适应度
            population[j].cfitness = 0;    //积累适应度
            population[j].lower[i] = lbound;     //基因变量的上界
            population[j].upper[i] = ubound;      //基因变量的下界
            population[j].gene[i] = randval(population[j].lower[i], population[j].upper[i]);  //变量 
        }
    }
}

//随机数产生函数  
double randval(double low, double high)
{
    double val;

    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_x(-100, 100);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

    val = random_real_num_x();

    return(val);
}


//评价函数，可以由用户自定义，该函数取得每个基因的适应度    
void evaluate()
{
    int i, j;

    for (j = 0; j < POPSIZE; j++)
    {
        for (i = 0; i < NVARS; i++)
            x[i] = population[j].gene[i];

        for (i = 0; i < NVARS; i++)
        {
            x[1 * NVARS + i] = 0.0;
            //printf("%d,%Lf\n", 1 * NVARS + j,x[1* NVARS +j]); //输出x变量的取值情况
        }

        cec14_test_func(x, f, NVARS, POPSIZE, func_num); //调用CEC2014测试集函数

        //返回x1、x2、x3....xn在CEC2014计算后的适应度大小
        //printf("\nG = %d, N = %d, 适应度 = %lf\n", generation, k++, f[0]);
        population[j].fitness = f[0];  //得到每个种群的适应度

    }

}


//保存每次遗传后的最佳基因 
void keep_the_best()
{
    int mem;
    int i;
    cur_best = 0;

    //保存最佳个体的索引  
    for (mem = 0; mem < POPSIZE; mem++)
    {
        if (population[mem].fitness > population[POPSIZE].fitness)
        {
            cur_best = mem;
            population[POPSIZE].fitness = population[mem].fitness;
        }
    }

    //一旦找到种群的最佳个体，就拷贝他的基因   
    for (i = 0; i < NVARS; i++)
        population[POPSIZE].gene[i] = population[cur_best].gene[i];
}


//搜寻杰出个体函数：找出最好和最坏的个体。  
//如果某代的最好个体比前一代的最好个体要坏，那么后者将会取代当前种群的最坏个体 
void elitist()
{
    int i;
    double best, worst; //最好和最坏个体的适应度值
    int best_mem, worst_mem; //最好和最坏个体的索引
    best = population[0].fitness;
    worst = population[0].fitness;

    for (i = 0; i < POPSIZE - 1; ++i)
    {
        if (population[i].fitness > population[i + 1].fitness)
        {
            if (population[i].fitness >= best)
            {
                best = population[i].fitness;
                best_mem = i;
            }
            if (population[i + 1].fitness <= worst)
            {
                worst = population[i + 1].fitness;
                worst_mem = i + 1;
            }
        }
        else
        {
            if (population[i].fitness <= worst)
            {
                worst = population[i].fitness;
                worst_mem = i;
            }
            if (population[i + 1].fitness >= best)
            {
                best = population[i + 1].fitness;
                best_mem = i + 1;
            }
        }
    }

    //如果新种群中的最好个体比前一代的最好个体要强的话，
    //那么就把新种群的最好个体拷贝出来。  
    //否则就用前一代的最好个体取代这次的最坏个体  
    if (best >= population[POPSIZE].fitness)
    {
        for (i = 0; i < NVARS; i++)
            population[POPSIZE].gene[i] = population[best_mem].gene[i];
        population[POPSIZE].fitness = population[best_mem].fitness;
    }
    else
    {
        for (i = 0; i < NVARS; i++)
            population[worst_mem].gene[i] = population[POPSIZE].gene[i];
        population[worst_mem].fitness = population[POPSIZE].fitness;
    }
}


//选择函数：用于最大化合并杰出模型的标准比例选择，保证最优秀的个体得以生存 
void select(void)
{
    int mem, j, i;
    double sum = 0;
    double p;

    //找出种群的适应度之和   
    for (mem = 0; mem < POPSIZE; mem++)
    {
        sum += population[mem].fitness;
    }

    //计算相对适应度    
    for (mem = 0; mem < POPSIZE; mem++)
    {
        population[mem].rfitness = population[mem].fitness / sum;
    }
    population[0].cfitness = population[0].rfitness;

    //计算累加适应度   
    for (mem = 1; mem < POPSIZE; mem++)
    {
        population[mem].cfitness = population[mem - 1].cfitness + population[mem].rfitness;
    }

    //用累加适应度作出选择  
    for (i = 0; i < POPSIZE; i++)
    {
        p = random_real_num_r();

        if (p < population[0].cfitness)
            newpopulation[i] = population[0];
        else
        {
            for (j = 0; j < POPSIZE; j++)
                if (p >= population[j].cfitness && p < population[j + 1].cfitness)
                    newpopulation[i] = population[j + 1];
        }
    }

    //当一个新种群建立的时候，将其拷贝回去  
    for (i = 0; i < POPSIZE; i++)
        population[i] = newpopulation[i];
}


//杂交函数：选择两个个体来杂交，这里用单点杂交 
void crossover(void)
{
    int mem, one;
    int first = 0;
    double x;

    for (mem = 0; mem < POPSIZE; ++mem)
    {
        x = random_real_num_r();

        if (x < PXOVER)
        {
            ++first;
            if (first % 2 == 0)
                Xover(one, mem);
            else
                one = mem;
        }
    }
}


//交叉   
void Xover(int one, int two)
{
    int i;
    int point;

    if (NVARS > 1)
    {
        if (NVARS == 2)
            point = 1;
        else {
            point = (rand() % (NVARS - 1)) + 1;
        }

        for (i = 0; i < point; i++)
            swap(&population[one].gene[i], &population[two].gene[i]);
    }
}

//交换操作
void swap(double* x, double* y)
{
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

//变异函数：被该函数选中后会使得某一变量被一个随机的值所取代 
void mutate(void)
{
    int i, j;
    double x;

    for (i = 0; i < POPSIZE; i++)
        for (j = 0; j < NVARS; j++)
        {
            x = random_real_num_r();

            if (x < PMUTATION)
            {
                population[i].gene[j] = randval(lbound, ubound);
            }
        }
}


//报告模拟进展情况，输出文件中的数据用空格隔开
void report(void)
{
    int i;
    double best_val;  //最佳种群适应度
    double avg;        //平均种群适应度
    double stddev;     //种群适应度偏差 
    double sum_square;  //各个个体平方之和
    double square_sum;  //平均值的平方乘个数
    double sum;        //所有种群适应度之和
    sum = 0.0;
    sum_square = 0.0;

    for (i = 0; i < POPSIZE; i++)
    {
        sum += population[i].fitness;
        sum_square += population[i].fitness * population[i].fitness;
    }       avg = sum / (double)POPSIZE;
    square_sum = avg * avg * POPSIZE;
    stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
    best_val = population[POPSIZE].fitness;

    //以间隔每2000代进行输出
    if (generation == 2000 || generation == 4000 || generation == 6000 || generation == 8000 || generation == 10000 ||
        generation == 12000 || generation == 14000 || generation == 16000 || generation == 18000 || generation == 20000)
    {

        printf("\n F%d\t第%d代\t%E\n", func_num, generation, best_val);
        fprintf(txtFile, "\n %5d  %E  %E  %E \n", generation, best_val, avg, stddev);

    }
}


//功能函数
void start_Function(int n) {

    char result_file[100];
    int i;
    func_num = n;
    generation = 0;

    //结果txt文件输出操作
    sprintf(result_file, "Result_F%d_P%d_PX%lf_PM%lf.txt", func_num, POPSIZE, PXOVER, PMUTATION);
    if ((txtFile = fopen(result_file, "w")) == NULL)
    {
        printf("文件打开失败！");
        exit(1);
    }

    printf("\n F%d 正在计算，请稍后...\n", func_num);

    fprintf(txtFile, "\n generation     best          average           standard \n");
    fprintf(txtFile, " number         value          fitness            deviation \n");

    //initialize();

    evaluate();    //评价函数，可以由用户自定义，该函数取得每个基因的适应度
    keep_the_best();    //保存每次遗传后的最佳基因

    while (generation < MAXGENS)
    {
        generation++;
        select();     //选择函数：用于最大化合并杰出模型的标准比例选择，保证最优秀的个体得以生存
        crossover();  //杂交函数：选择两个个体来杂交，这里用单点杂交 
        mutate();     //变异函数：被该函数选中后会使得某一变量被一个随机的值所取代 
        report();     //报告模拟进展情况
        evaluate();   //评价函数，可以由用户自定义，该函数取得每个基因的适应度
        elitist();    //搜寻杰出个体函数：找出最好和最坏的个体。如果某代的最好个体比前一代的最好个体要坏，那么后者将会取代当前种群的最坏个体 
    }

    fprintf(txtFile, "\n计算完成\n");
    fprintf(txtFile, "\n变量: \n");
    for (i = 0; i < NVARS; i++)
    {
        fprintf(txtFile, "\n x%d = %LF", i + 1, population[POPSIZE].gene[i]);
    }
    fprintf(txtFile, "\n\n最佳适应度 = %LF", population[POPSIZE].fitness);
    fclose(txtFile);

    population[POPSIZE].fitness = 0;  //修正最终适应度结果

    printf("\n F%d 计算完成，结果已生成：%s\n", func_num, result_file);

}


//测试种群大小对算法的性能影响
//POPSIZE [20,200]
//PXOVER = 0.7 , PMUTATION = 0.06
void test_popsize(double pxover, double pmutation)
{
    int popsize;

    PXOVER = pxover;
    PMUTATION = pmutation;

    for (popsize = 20; popsize <= 200; popsize += 20)
    {
        POPSIZE = popsize;
        printf("\n\n POPSIZE = %d, PMUTATION = %lf, PXOVER = %lf\n", POPSIZE, PMUTATION, PXOVER);

        initialize();   //随机数随种群大小改变

        //测试随机数的生成值
       /* for (int j = 0; j < POPSIZE; j++)
        {
            for (int i = 0; i < NVARS; i++)
            {
                printf("P = %d, N = %d, x%d = %lf\n", j+1,i+1, i+1,population[j].gene[i]);
            }
        }*/

        //开启测试函数，如不进行某测试函数，通过"//"注释掉即可
        start_Function(1); //Rotated High Conditioned Elliptic Funcktion

        start_Function(3); //Rotated Discus Function

        start_Function(8); //Shifted Rastrigin's Function

        start_Function(12); //Shiftedand Rotated Katsuura Function

        start_Function(18); //Hybrid Function 2 (N=3) 

        start_Function(20); //Hybrid Function 4 (N=4) 

        start_Function(25); //Composition Function 3 (N=3) 

        start_Function(30); //Composition Function 8 (N=3) 

        //释放内存空间
        free(x);
        free(f);
        delete[]population;
        delete[]newpopulation;
    }
}


//测试交叉大小对算法的性能影响 
//PXOVER [0.9.0.1]
//POPSIZE = 50, PMUTATION = 0.06
void test_pxover(int popsize, double pmutation)
{
    double pxover;

    POPSIZE = popsize;
    PMUTATION = pmutation;

    initialize();   //保持随机数不变

    for (pxover = 0.9; pxover >= 0.1; pxover -= 0.1)
    {
        PXOVER = pxover;
        printf("\n\n POPSIZE = %d, PMUTATION = %lf, PXOVER = %lf\n", POPSIZE, PMUTATION, PXOVER);

        //测试随机数的生成值
        /* for (int j = 0; j < POPSIZE; j++)
         {
             for (int i = 0; i < NVARS; i++)
             {
                 printf("P = %d, N = %d, x%d = %lf\n", j + 1, i + 1, i + 1, population[j].gene[i]);
             }
         }*/

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

    //释放内存空间
    free(x);
    free(f);
    delete[]population;
    delete[]newpopulation;
}


//测试变异大小对算法的性能影响
//PMUTATION [0.01.0.1]
//POPSIZE = 50, PXOVER = 0.7
void test_pmutation(int popsize, double pxover)
{
    double pmutation;

    POPSIZE = popsize;
    PXOVER = pxover;

    initialize();   //保持随机数不变

    for (pmutation = 0.01; pmutation <= 0.1; pmutation += 0.01)
    {
        PMUTATION = pmutation;
        printf("\n\n POPSIZE = %d, PMUTATION = %lf, PXOVER = %lf\n", POPSIZE, PMUTATION, PXOVER);

        //测试随机数的生成值
        /* for (int j = 0; j < POPSIZE; j++)
         {
             for (int i = 0; i < NVARS; i++)
             {
                 printf("P = %d, N = %d, x%d = %lf\n", j + 1, i + 1, i + 1, population[j].gene[i]);
             }
         }*/

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

    //释放内存空间
    free(x);
    free(f);
    delete[]population;
    delete[]newpopulation;
}


//主函数
void main(void) {

    int popsize;
    double pxover, pmutation;

    popsize = 50;   //种群大小 [20,200]
    pxover = 0.7;   //交叉率  [0.9,0.1]
    pmutation = 0.06;   //变异率 [0.01,0.1]

    test_popsize(pxover, pmutation); //测试种群大小对算法的性能影响

    test_pxover(popsize, pmutation); //测试交叉大小对算法的性能影响

    test_pmutation(popsize, pxover); //测试变异大小对算法的性能影响

}