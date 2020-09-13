/*
	CEC14 Test Function Suite for Differential Evolution
	Tongbo Cai(email: 461892420@qq.com; 461892420@njtech.edu.cn) 
	Sep. 13th 2020
*/

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <boost/random.hpp>

#define NVARS 10	//问题变量个数
#define MAXGENS 500	//最大迭代次数
#define NP 50	//种群规模
#define F 0.5	//变异的控制参数
#define CR 0.3	//交叉的控制参数


struct individual	//个体
{
	double position[NVARS];
	double fitness;
};

struct genotype	//基因
{
	struct individual ptcle[NP];
	double f;
	double cr;
	double(*function)(double*); //通过指针存储适应度
	double global_fitness;
	double global_solution[NVARS];
	double lbound;
	double ubound;
};

struct individual mutation[NP];
struct individual crossover[NP];
struct genotype population;

double* f; //返回cec计算结果
double* OShift, * M, * y, * z, * x_bound;   //声明CEC2014外部定义
int ini_flag = 0, n_flag, func_flag, * SS;  //声明CEC2014外部定义
int func_num; //声明CEC2014功能函数
FILE* txtFile;  //结果输出文件

double lboud = -100; //最小值
double ubound = 100; //最大值

void cec14_test_func(double*, double*, int, int, int);  //声明CEC2014函数测试集


//适应度计算函数
double function_fitness(double* var)
{
	double fitness;
	double* x = var; //x数组

	cec14_test_func(x, f, NVARS, NP, func_num); //调用CEC2014测试集函数
	fitness = f[0];

	return fitness;
}

//(0,1)随机数产生函数  
double RandX()
{
	double val;

	boost::mt19937 generator(time(0) * rand());
	boost::uniform_real<> uniform_real_generate_x(0, 1);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

	val = random_real_num_x();

	return(val);
}

//(0,NVARS)随机数产生函数  
int RandNVARS()
{
	int val;

	boost::mt19937 generator(time(0) * rand());
	boost::uniform_real<> uniform_real_generate_x(0, NVARS);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

	val = (int)random_real_num_x();

	return(val);
}

//(0,NP)随机数产生函数  
int RandNP()
{
	int val;

	boost::mt19937 generator(time(0) * rand());
	boost::uniform_real<> uniform_real_generate_x(0, NP);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

	val = (int)random_real_num_x();

	return(val);
}

//初始化
void Initialize()
{
	f = (double*)malloc(sizeof(double) * NP);	//适应度存放

	population.f = F;
	population.cr = CR;
	population.lbound = lboud;
	population.ubound = ubound;
	population.function = function_fitness;

	srand((unsigned)time(NULL));

	for (int i = 0; i < NP; i++)
	{
		for (int j = 0; j < NVARS; j++)
		{
			double randx = RandX();
			population.ptcle[i].position[j] = population.lbound + (population.ubound - population.lbound) * randx;
		}
		population.ptcle[i].fitness = population.function(population.ptcle[i].position);
	}

	memcpy(&population.global_solution, &population.ptcle[0].position, sizeof(double) * NVARS);
	population.global_fitness = population.ptcle[0].fitness;
}

//变异
void Mutation() {
	for (int i = 0; i < NP; i++)
	{
		int rand1 = RandNP();
		int rand2 = RandNP();
		int rand3 = RandNP();

		if (rand1 == i || rand2 == i || rand3 == i || rand1 == rand2 || rand2 == rand3)
		{

			rand1 = RandNP();
			rand2 = RandNP();
			rand3 = RandNP();
		}
		for (int j = 0; j < NVARS; j++)
		{
			mutation[i].position[j] = population.ptcle[rand1].position[j] + population.f * (population.ptcle[rand2].position[j] - population.ptcle[rand3].position[j]);
			if (mutation[i].position[j] > population.ubound)
			{
				mutation[i].position[j] = population.ubound;
			}
			if (mutation[i].position[j] < population.lbound)
			{
				mutation[i].position[j] = population.lbound;
			}
		}
	}

}

//交叉
void CrossOver() {
	for (int i = 0; i < NP; i++)
	{
		int randc = RandNVARS();
		for (int j = 0; j < NVARS; j++)
		{
			double rand_cr = RandX();
			if ((j == randc) || (rand_cr <= population.cr))
			{
				crossover[i].position[j] = mutation[i].position[j];
			}
			else
			{
				crossover[i].position[j] = population.ptcle[i].position[j];
			}
			if (crossover[i].position[j] > population.ubound)
			{
				crossover[i].position[j] = population.ubound;
			}
			if (crossover[i].position[j] < population.lbound)
			{
				crossover[i].position[j] = population.lbound;
			}
		}
		crossover[i].fitness = population.function(crossover[i].position);
	}
}

//选择
void Selection() {
	for (int i = 0; i < NP; i++)
	{
		if (crossover[i].fitness < population.ptcle[i].fitness)
		{
			population.ptcle[i].fitness = crossover[i].fitness;
			memcpy(&population.ptcle[i].position, &crossover[i].position, sizeof(double) * NVARS);

			if (population.ptcle[i].fitness < population.global_fitness)
			{
				population.global_fitness = population.ptcle[i].fitness;
				memcpy(&population.global_solution, &population.ptcle[i].position, sizeof(double) * NVARS);
			}
		}

	}

}

//功能函数
void start_Function(int n)
	
{	
	char result_file[100];
	func_num = n;
	clock_t pro_start, pro_finished;


	//结果txt文件输出操作
	sprintf(result_file, "result/Result_F%d_NP%d_F%lf_CR%lf.txt", func_num, NP, F, CR);
	if ((txtFile = fopen(result_file, "w")) == NULL)
	{
		printf("文件打开失败！");
		exit(1);
	}

	printf("\n *************** NP = %d, P = %lf, CR = %lf ***************\n ", NP, F, CR);
	fprintf(txtFile, "\n*************** NP = %d, P = %lf, CR = %lf ***************\n", NP, F, CR);

	printf("\n F%d开始计算，请稍后...\n", func_num);
	pro_start = clock();

	Initialize();

	for (int generation = 1; generation <= MAXGENS; generation++)
	{
		Mutation();
		CrossOver();
		Selection();
		// 一次迭代结束
		printf("\n F%d\t第%d代\t%E", func_num,generation, population.global_fitness);
		fprintf(txtFile, "\nF%d\t第%d代\t%E", func_num, generation, population.global_fitness);

	}

	fprintf(txtFile, "\n\nF%d 计算完成\n",func_num);
	printf("\n\n 最佳适应度：%LF\n", population.global_fitness);
	fprintf(txtFile, "\n\n最佳适应度：%LF\n", population.global_fitness);
	printf("\n 变量: \n");
	fprintf(txtFile, "\n变量: \n");

	for (int i = 0; i < NVARS; i++)
	{
		printf("\n x%d = %LF", i + 1, population.global_solution[i]);
		fprintf(txtFile, "\nx%d = %LF", i + 1, population.global_solution[i]);

	}
	printf("\n\n F%d 计算完成，结果已保存至：%s\n", func_num, result_file);
	fclose(txtFile);

	//计算耗时
	pro_finished = clock();
	double Total_time = (double)(pro_finished - pro_start) / CLOCKS_PER_SEC;
	printf("\n 运行耗时: %3f second\n\n", Total_time);
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