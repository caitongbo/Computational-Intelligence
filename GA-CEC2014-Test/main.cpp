/*
    ����GA��CEC2014�Ϻ���F1��F3��F8��F12��F18��F20��F25��F30�ϵ��㷨����
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

//����������Ҫ���޸����²���   
#define MAXGENS 20000 //������������ÿ2000����¼һ�� [0,20000] 
const int NVARS = 10; //����������� [2,10,20,30,50,100]

#define TRUE 1 
#define FALSE 0    
#define PI 3.1415926

struct genotype //��Ⱥ��һ������Ľṹ������
{
    double gene[NVARS]; //����
    double fitness; //�������Ӧ��
    double upper[NVARS]; //����������Ͻ�
    double lower[NVARS]; //����������½�
    double rfitness; //�Ƚ���Ӧ��
    double cfitness; //������Ӧ��
};

struct genotype* population; //��Ⱥ
struct genotype* newpopulation; //����Ⱥȡ���ɵĻ���

//������һЩ�������� 
void initialize(void);  //��ʼ���������
double randval(double, double);    //�������������
void evaluate(void);             //���ۺ������������û��Զ��壬�ú���ȡ��ÿ���������Ӧ��
void keep_the_best(void);       //����ÿ���Ŵ������ѻ���
void elitist(void);  //��Ѱ�ܳ����庯�����ҳ���ú���ĸ��塣���ĳ������ø����ǰһ������ø���Ҫ������ô���߽���ȡ����ǰ��Ⱥ������� 
void select(void);  //ѡ������������󻯺ϲ��ܳ�ģ�͵ı�׼����ѡ�񣬱�֤������ĸ����������
void crossover(void);  //�ӽ�������ѡ�������������ӽ��������õ����ӽ�
void Xover(int, int);   //���� 
void swap(double*, double*);  //����
void mutate(void);    //���캯�������ú���ѡ�к��ʹ��ĳһ������һ�������ֵ��ȡ�� 
void report(void);   //����ģ���չ���
void start_Function(int);   //�������ܺ���
void test_popsize(void); //������Ⱥ��С���㷨������Ӱ�� [20,200]
void test_pxover(void);  //���Խ����С���㷨������Ӱ�� [0.9,0.1]
void test_pmutation(void);   //���Ա����С���㷨������Ӱ�� [0.01,0.1]
void cec14_test_func(double*, double*, int, int, int);  //CEC2014�������Լ�

double* OShift, * M, * y, * z, * x_bound;   //CEC2014�ⲿ����
double* f, * x; //f����Ӧ�ȴ�С��x��δ֪����ֵ
double PXOVER; //������� [0.9,0.1] ����-0.1
double PMUTATION; //������� [0.01,0.1] ����+0.01
int POPSIZE; //��Ⱥ��С [20,200] ����+20
int ini_flag = 0, n_flag, func_flag, * SS;  //CEC2014�ⲿ����
int func_num; //��ǹ��ܺ�����1��3��8��12��18��20��25��30
int generation; //��ǰ�������
int cur_best; //���Ÿ���
FILE* txtFile;  //�������ļ�

double lbound = -100, ubound = 100; //������������½� (-100,100)

//ȫ������������������ڱ��졢����Ȳ���
boost::mt19937 generator(time(0)* rand());
boost::uniform_real<> uniform_real_generate_r(0, 1);
boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r(generator, uniform_real_generate_r);


//�������ʼ��
void initialize(void)
{
    int i, j;

    x = (double*)malloc(NVARS * POPSIZE * sizeof(double));
    f = (double*)malloc(sizeof(double) * POPSIZE);

    population = new genotype[POPSIZE + 1]; //��Ⱥ
    newpopulation = new genotype[POPSIZE + 1]; //����Ⱥ,ȡ���ɵĻ���

    for (i = 0; i < NVARS; i++)
    {
        for (j = 0; j < POPSIZE; j++)
        {
            population[j].fitness = 0;     //�������Ӧ��
            population[j].rfitness = 0;    //�Ƚ���Ӧ��
            population[j].cfitness = 0;    //������Ӧ��
            population[j].lower[i] = lbound;     //����������Ͻ�
            population[j].upper[i] = ubound;      //����������½�
            population[j].gene[i] = randval(population[j].lower[i], population[j].upper[i]);  //���� 
        }
    }
}

//�������������  
double randval(double low, double high)
{
    double val;

    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_x(-100, 100);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

    val = random_real_num_x();

    return(val);
}


//���ۺ������������û��Զ��壬�ú���ȡ��ÿ���������Ӧ��    
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
            //printf("%d,%Lf\n", 1 * NVARS + j,x[1* NVARS +j]); //���x������ȡֵ���
        }

        cec14_test_func(x, f, NVARS, POPSIZE, func_num); //����CEC2014���Լ�����

        //����x1��x2��x3....xn��CEC2014��������Ӧ�ȴ�С
        //printf("\nG = %d, N = %d, ��Ӧ�� = %lf\n", generation, k++, f[0]);
        population[j].fitness = f[0];  //�õ�ÿ����Ⱥ����Ӧ��

    }

}


//����ÿ���Ŵ������ѻ��� 
void keep_the_best()
{
    int mem;
    int i;
    cur_best = 0;

    //������Ѹ��������  
    for (mem = 0; mem < POPSIZE; mem++)
    {
        if (population[mem].fitness > population[POPSIZE].fitness)
        {
            cur_best = mem;
            population[POPSIZE].fitness = population[mem].fitness;
        }
    }

    //һ���ҵ���Ⱥ����Ѹ��壬�Ϳ������Ļ���   
    for (i = 0; i < NVARS; i++)
        population[POPSIZE].gene[i] = population[cur_best].gene[i];
}


//��Ѱ�ܳ����庯�����ҳ���ú���ĸ��塣  
//���ĳ������ø����ǰһ������ø���Ҫ������ô���߽���ȡ����ǰ��Ⱥ������� 
void elitist()
{
    int i;
    double best, worst; //��ú���������Ӧ��ֵ
    int best_mem, worst_mem; //��ú�����������
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

    //�������Ⱥ�е���ø����ǰһ������ø���Ҫǿ�Ļ���
    //��ô�Ͱ�����Ⱥ����ø��忽��������  
    //�������ǰһ������ø���ȡ����ε������  
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


//ѡ������������󻯺ϲ��ܳ�ģ�͵ı�׼����ѡ�񣬱�֤������ĸ���������� 
void select(void)
{
    int mem, j, i;
    double sum = 0;
    double p;

    //�ҳ���Ⱥ����Ӧ��֮��   
    for (mem = 0; mem < POPSIZE; mem++)
    {
        sum += population[mem].fitness;
    }

    //���������Ӧ��    
    for (mem = 0; mem < POPSIZE; mem++)
    {
        population[mem].rfitness = population[mem].fitness / sum;
    }
    population[0].cfitness = population[0].rfitness;

    //�����ۼ���Ӧ��   
    for (mem = 1; mem < POPSIZE; mem++)
    {
        population[mem].cfitness = population[mem - 1].cfitness + population[mem].rfitness;
    }

    //���ۼ���Ӧ������ѡ��  
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

    //��һ������Ⱥ������ʱ�򣬽��俽����ȥ  
    for (i = 0; i < POPSIZE; i++)
        population[i] = newpopulation[i];
}


//�ӽ�������ѡ�������������ӽ��������õ����ӽ� 
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


//����   
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

//��������
void swap(double* x, double* y)
{
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

//���캯�������ú���ѡ�к��ʹ��ĳһ������һ�������ֵ��ȡ�� 
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


//����ģ���չ���������ļ��е������ÿո����
void report(void)
{
    int i;
    double best_val;  //�����Ⱥ��Ӧ��
    double avg;        //ƽ����Ⱥ��Ӧ��
    double stddev;     //��Ⱥ��Ӧ��ƫ�� 
    double sum_square;  //��������ƽ��֮��
    double square_sum;  //ƽ��ֵ��ƽ���˸���
    double sum;        //������Ⱥ��Ӧ��֮��
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

    //�Լ��ÿ2000���������
    if (generation == 2000 || generation == 4000 || generation == 6000 || generation == 8000 || generation == 10000 ||
        generation == 12000 || generation == 14000 || generation == 16000 || generation == 18000 || generation == 20000)
    {

        printf("\n F%d\t��%d��\t%E\n", func_num, generation, best_val);
        fprintf(txtFile, "\n %5d  %E  %E  %E \n", generation, best_val, avg, stddev);

    }
}


//���ܺ���
void start_Function(int n) {

    char result_file[100];
    int i;
    func_num = n;
    generation = 0;

    //���txt�ļ��������
    sprintf(result_file, "Result_F%d_P%d_PX%lf_PM%lf.txt", func_num, POPSIZE, PXOVER, PMUTATION);
    if ((txtFile = fopen(result_file, "w")) == NULL)
    {
        printf("�ļ���ʧ�ܣ�");
        exit(1);
    }

    printf("\n F%d ���ڼ��㣬���Ժ�...\n", func_num);

    fprintf(txtFile, "\n generation     best          average           standard \n");
    fprintf(txtFile, " number         value          fitness            deviation \n");

    //initialize();

    evaluate();    //���ۺ������������û��Զ��壬�ú���ȡ��ÿ���������Ӧ��
    keep_the_best();    //����ÿ���Ŵ������ѻ���

    while (generation < MAXGENS)
    {
        generation++;
        select();     //ѡ������������󻯺ϲ��ܳ�ģ�͵ı�׼����ѡ�񣬱�֤������ĸ����������
        crossover();  //�ӽ�������ѡ�������������ӽ��������õ����ӽ� 
        mutate();     //���캯�������ú���ѡ�к��ʹ��ĳһ������һ�������ֵ��ȡ�� 
        report();     //����ģ���չ���
        evaluate();   //���ۺ������������û��Զ��壬�ú���ȡ��ÿ���������Ӧ��
        elitist();    //��Ѱ�ܳ����庯�����ҳ���ú���ĸ��塣���ĳ������ø����ǰһ������ø���Ҫ������ô���߽���ȡ����ǰ��Ⱥ������� 
    }

    fprintf(txtFile, "\n�������\n");
    fprintf(txtFile, "\n����: \n");
    for (i = 0; i < NVARS; i++)
    {
        fprintf(txtFile, "\n x%d = %LF", i + 1, population[POPSIZE].gene[i]);
    }
    fprintf(txtFile, "\n\n�����Ӧ�� = %LF", population[POPSIZE].fitness);
    fclose(txtFile);

    population[POPSIZE].fitness = 0;  //����������Ӧ�Ƚ��

    printf("\n F%d ������ɣ���������ɣ�%s\n", func_num, result_file);

}


//������Ⱥ��С���㷨������Ӱ��
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

        initialize();   //���������Ⱥ��С�ı�

        //���������������ֵ
       /* for (int j = 0; j < POPSIZE; j++)
        {
            for (int i = 0; i < NVARS; i++)
            {
                printf("P = %d, N = %d, x%d = %lf\n", j+1,i+1, i+1,population[j].gene[i]);
            }
        }*/

        //�������Ժ������粻����ĳ���Ժ�����ͨ��"//"ע�͵�����
        start_Function(1); //Rotated High Conditioned Elliptic Funcktion

        start_Function(3); //Rotated Discus Function

        start_Function(8); //Shifted Rastrigin's Function

        start_Function(12); //Shiftedand Rotated Katsuura Function

        start_Function(18); //Hybrid Function 2 (N=3) 

        start_Function(20); //Hybrid Function 4 (N=4) 

        start_Function(25); //Composition Function 3 (N=3) 

        start_Function(30); //Composition Function 8 (N=3) 

        //�ͷ��ڴ�ռ�
        free(x);
        free(f);
        delete[]population;
        delete[]newpopulation;
    }
}


//���Խ����С���㷨������Ӱ�� 
//PXOVER [0.9.0.1]
//POPSIZE = 50, PMUTATION = 0.06
void test_pxover(int popsize, double pmutation)
{
    double pxover;

    POPSIZE = popsize;
    PMUTATION = pmutation;

    initialize();   //�������������

    for (pxover = 0.9; pxover >= 0.1; pxover -= 0.1)
    {
        PXOVER = pxover;
        printf("\n\n POPSIZE = %d, PMUTATION = %lf, PXOVER = %lf\n", POPSIZE, PMUTATION, PXOVER);

        //���������������ֵ
        /* for (int j = 0; j < POPSIZE; j++)
         {
             for (int i = 0; i < NVARS; i++)
             {
                 printf("P = %d, N = %d, x%d = %lf\n", j + 1, i + 1, i + 1, population[j].gene[i]);
             }
         }*/

         //�������Ժ������粻����ĳ���Ժ�����ͨ��"//"ע�͵�����
        start_Function(1); //Rotated High Conditioned Elliptic Funcktion

        start_Function(3); //Rotated Discus Function

        start_Function(8); //Shifted Rastrigin's Function

        start_Function(12); //Shiftedand Rotated Katsuura Function

        start_Function(18); //Hybrid Function 2 (N=3) 

        start_Function(20); //Hybrid Function 4 (N=4) 

        start_Function(25); //Composition Function 3 (N=3) 

        start_Function(30); //Composition Function 8 (N=3) 

    }

    //�ͷ��ڴ�ռ�
    free(x);
    free(f);
    delete[]population;
    delete[]newpopulation;
}


//���Ա����С���㷨������Ӱ��
//PMUTATION [0.01.0.1]
//POPSIZE = 50, PXOVER = 0.7
void test_pmutation(int popsize, double pxover)
{
    double pmutation;

    POPSIZE = popsize;
    PXOVER = pxover;

    initialize();   //�������������

    for (pmutation = 0.01; pmutation <= 0.1; pmutation += 0.01)
    {
        PMUTATION = pmutation;
        printf("\n\n POPSIZE = %d, PMUTATION = %lf, PXOVER = %lf\n", POPSIZE, PMUTATION, PXOVER);

        //���������������ֵ
        /* for (int j = 0; j < POPSIZE; j++)
         {
             for (int i = 0; i < NVARS; i++)
             {
                 printf("P = %d, N = %d, x%d = %lf\n", j + 1, i + 1, i + 1, population[j].gene[i]);
             }
         }*/

         //�������Ժ������粻����ĳ���Ժ�����ͨ��"//"ע�͵�����
        start_Function(1); //Rotated High Conditioned Elliptic Funcktion

        start_Function(3); //Rotated Discus Function

        start_Function(8); //Shifted Rastrigin's Function

        start_Function(12); //Shiftedand Rotated Katsuura Function

        start_Function(18); //Hybrid Function 2 (N=3) 

        start_Function(20); //Hybrid Function 4 (N=4) 

        start_Function(25); //Composition Function 3 (N=3) 

        start_Function(30); //Composition Function 8 (N=3) 

    }

    //�ͷ��ڴ�ռ�
    free(x);
    free(f);
    delete[]population;
    delete[]newpopulation;
}


//������
void main(void) {

    int popsize;
    double pxover, pmutation;

    popsize = 50;   //��Ⱥ��С [20,200]
    pxover = 0.7;   //������  [0.9,0.1]
    pmutation = 0.06;   //������ [0.01,0.1]

    test_popsize(pxover, pmutation); //������Ⱥ��С���㷨������Ӱ��

    test_pxover(popsize, pmutation); //���Խ����С���㷨������Ӱ��

    test_pmutation(popsize, pxover); //���Ա����С���㷨������Ӱ��

}