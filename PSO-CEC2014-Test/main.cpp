/*
 * ����PSO��CEC2014�Ϻ���F1��F3��F8��F12��F18��F20��F25��F30�ϵ��㷨����
 * 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <boost/random.hpp>

#define c1 2 //���ٶ�����
#define c2 2
#define MAXGENS 500  // ��������
#define POPSIZE 50 // ��Ⱥ��ģ
#define ubound 100 // �������ȡֵ
#define lbound -100 // ������Сȡֵ
#define Vmax 1 // �ٶ����ֵ
#define Vmin -1 //�ٶ���Сֵ
#define NVARS 10 // ���ӵ�ά��
#define PI 3.1415926 //Բ����

double pop[POPSIZE][NVARS]; // ������Ⱥ����
double V[POPSIZE][NVARS]; // ������Ⱥ�ٶ�����
double fitness[POPSIZE]; // ������Ⱥ����Ӧ������
double result[MAXGENS];  //������ÿ�ε�����Ⱥ����ֵ������
double pbest[POPSIZE][NVARS];  // ���弫ֵ��λ��
double gbest[NVARS]; //Ⱥ�弫ֵ��λ��
double fitnesspbest[POPSIZE]; //���弫ֵ��Ӧ�ȵ�ֵ
double fitnessgbest; // Ⱥ�弫ֵ��Ӧ��ֵ
double genbest[MAXGENS][NVARS]; //ÿһ������ֵȡֵ����
double* f; //����cec������
FILE* txtFile;

double* OShift, * M, * y, * z, * x_bound;   //CEC2014�ⲿ����
int ini_flag = 0, n_flag, func_flag, * SS;  //CEC2014�ⲿ����

int func_num; //CEC2014���ܺ���

void cec14_test_func(double*, double*, int, int, int);  //CEC2014�������Լ�


//��Ӧ�Ⱥ���
double func(double* arr)
{
    double *x = arr; //x����

    cec14_test_func(x, f, NVARS, POPSIZE, func_num); //����CEC2014���Լ�����

    double fitness = f[0];

    return fitness;

}

//(0,1)�������������  
double randval_x()
{
    double val;

    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_x(0, 1);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

    val = random_real_num_x();

    return(val);
}

//��Ⱥ�������������  
double randval_pop()
{
    double val;

    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_x(-100, 100);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

    val = random_real_num_x();

    return(val);
}

//�ٶ��������������  
double randval_v()
{
    double val;

    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_x(-1, 1);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

    val = random_real_num_x();

    return(val);
}
// ��Ⱥ��ʼ��
void pop_init(void)
{
    f = (double*)malloc(sizeof(double) * POPSIZE);

    for (int i = 0;i < POPSIZE;i++)	//POPSIZE ��Ⱥ��ģ
    {
        for (int j = 0;j < NVARS;j++)//NVARS ���ӵ�ά��
        {
            pop[i][j] = randval_pop(); //-100��100֮��������   //pop[i][j]  ��Ⱥ����
            V[i][j] = randval_v(); //-1��1֮��             // V[i][j] ��Ⱥ�ٶ�����
        }
        fitness[i] = func(pop[i]); //������Ӧ�Ⱥ���ֵ

    }
}
// max()��������
double* max(double* fit, int size)//Ѱ�������ڵ����ֵ�����ֵ������λ��
{
    int index = 0; // ��ʼ�����
    double max = *fit; // ��ʼ�����ֵΪ�����һ��Ԫ��
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
// ����Ѱ��
void PSO_func(void)
{
    pop_init();
    double* best_fit_index; // ���ڴ��Ⱥ�弫ֵ����λ��(���)
    best_fit_index = max(fitness, POPSIZE); //��Ⱥ�弫ֵ
    int index = (int)(*best_fit_index);

    char result_file[100];

    //���txt�ļ��������
    sprintf(result_file, "result/Result_F%d_P%d.txt", func_num, POPSIZE);
    if ((txtFile = fopen(result_file, "w")) == NULL)
    {
        printf("�ļ���ʧ�ܣ�");
        exit(1);
    }

    // Ⱥ�弫ֵλ��
    for (int i = 0;i < NVARS;i++)
    {
        gbest[i] = pop[index][i];
    }
    // ���弫ֵλ��
    for (int i = 0;i < POPSIZE;i++)
    {
        for (int j = 0;j < NVARS;j++)
        {
            pbest[i][j] = pop[i][j];
        }
    }
    // ���弫ֵ��Ӧ��ֵ
    for (int i = 0;i < POPSIZE;i++)
    {
        fitnesspbest[i] = fitness[i];
    }
    //Ⱥ�弫ֵ��Ӧ��ֵ
    double bestfitness = *(best_fit_index + 1);
    fitnessgbest = bestfitness;

    //����Ѱ��
    for (int i = 0;i < MAXGENS;i++)
    {

        for (int j = 0;j < POPSIZE;j++)
        {
            //�ٶȸ��¼����Ӹ���
            for (int k = 0;k < NVARS;k++)
            {
                // �ٶȸ���
                double rand1 = randval_x(); //0��1֮��������
                double rand2 = randval_x(); //0��1֮��������
                V[j][k] = V[j][k] + c1 * rand1 * (pbest[j][k] - pop[j][k]) + c2 * rand2 * (gbest[k] - pop[j][k]);
                if (V[j][k] > Vmax)
                    V[j][k] = Vmax;
                if (V[j][k] < Vmin)
                    V[j][k] = Vmin;
                // ���Ӹ���
                pop[j][k] = pop[j][k] + V[j][k];
                if (pop[j][k] > ubound)
                    pop[j][k] = ubound;
                if (pop[j][k] < lbound)
                    pop[j][k] = lbound;
            }
            fitness[j] = func(pop[j]); //�����ӵ���Ӧ��ֵ
        }
        for (int j = 0;j < POPSIZE;j++)
        {
            // ���弫ֵ����
            if (fitness[j] > fitnesspbest[j])
            {
                for (int k = 0;k < NVARS;k++)
                {
                    pbest[j][k] = pop[j][k];
                }
                fitnesspbest[j] = fitness[j];
            }
            // Ⱥ�弫ֵ����
            if (fitness[j] > fitnessgbest)
            {
                for (int k = 0;k < NVARS;k++)
                    gbest[k] = pop[j][k];
                fitnessgbest = fitness[j];
            }
        }
        for (int k = 0;k < NVARS;k++)
        {
            genbest[i][k] = gbest[k]; // ÿһ������ֵȡֵ����λ�ü�¼
        }
        result[i] = fitnessgbest; // ÿ��������ֵ��¼������

        printf("\n F%d\t��%d��\t%E", func_num, i + 1, result[i]);
        fprintf(txtFile,"\n F%d\t��%d��\t%E", func_num,i + 1, result[i]);
    }
}

//���ܺ���
void start_Function(int n) 
{
    char result_file[100];

    func_num = n;
    
    printf("\nF%d��ʼ���㣬���Ժ�...\n",func_num);

    clock_t start, finish; //����ʼ�ͽ���ʱ��
    start = clock(); //��ʼ��ʱ

    PSO_func();

    double* best_arr;
    best_arr = max(result, MAXGENS);   //result  ���ÿ�ε�����Ⱥ����ֵ������   //MAXGENS ��������
    int best_gen_number = *best_arr; // ����ֵ�����Ĵ���
    double best = *(best_arr + 1); //����ֵ

    printf("\n\n������%d�Σ��ڵ�%d��ȡ������ֵ������ֵΪ:%lf.\n", MAXGENS, best_gen_number + 1, best);
    fprintf(txtFile,"\n\n������%d�Σ��ڵ�%d��ȡ������ֵ������ֵΪ:%lf.\n", MAXGENS, best_gen_number + 1, best);

    printf("ȡ������ֵ��λ��Ϊ(%lf,%lf).\n", genbest[best_gen_number][0], genbest[best_gen_number][1]);
    fprintf(txtFile,"ȡ������ֵ��λ��Ϊ(%lf,%lf).\n", genbest[best_gen_number][0], genbest[best_gen_number][1]);

    finish = clock(); //����ʱ��
    double duration = (double)(finish - start) / CLOCKS_PER_SEC; // ��������ʱ��
    printf("�������к�ʱ:%lf\n", duration);
    fprintf(txtFile,"�������к�ʱ:%lf\n", duration);

    fclose(txtFile);
    printf("F%d������ɡ�\n", func_num);


}
// ������
int main(void)
{
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