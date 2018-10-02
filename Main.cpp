#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
#include "string.h"

#define POPSIZE	1000     
#define MAXNVARS	50			
#define MAXGENS	1000
#define RECORD_FRE 100 
#define total_task_num	10
#define RUNS	30
#define rmp		0.3
#define GLB	0
#define GUB 1
#define SBX_rate 0.9
#define mu 1.0
int job;

const double PI = acos(-1.0);
double MAX_GENE[MAXNVARS], MIN_GENE[MAXNVARS];
double fbest_value[total_task_num][MAXGENS / RECORD_FRE + 1][RUNS];

double	fbest[total_task_num];
int		generation;

struct gene
{
	double	x[MAXNVARS];
	double	f[total_task_num];
	int		tao;
    int		scale_f;
    int		rank[total_task_num];
};

struct gene population[2*POPSIZE];


int task_num;
int NVARS;
int cnt;


#include"benchmark.h"


double randval(double a, double b)
{
	return a + (b-a) * rand()/(double)RAND_MAX;
}

struct VI
{
	double v;
	int I;
};
struct VI vis[2 * POPSIZE];

int cmp_vi( const void *pa , const void *pb )
{
	struct VI * p = (struct VI *) pa;
	struct VI * q = (struct VI *) pb;
	return p->v > q->v ? 1 : -1;
}

double objective(double x[], int g)
{
	double f = 0;
	if (g == 0){
		f = Sphere(x, g);
	}
	else if (g == 1){
		f = Sphere2(x, g);
	}
	else if (g == 2){
		f = Sphere3(x, g);
	}
	else if (g == 3){
		f = Weierstrass25D(x, g);
	}
	else if (g == 4){
		f = Rosenbrock(x, g);
	}
	else if (g == 5){
		f = Ackley(x, g);
	}
	else if (g == 6){
		f = Weierstrass50D(x, g);
	}
	else if (g == 7){
		f = Schwefel(x, g);
	}
	else if (g == 8) {
		f = Griewank(x, g);
	}
	else if (g == 9){
		f = Rastrgin(x, g);
	}
	return f;
}

static int phase = 0;
double gaussian()
{
	static double V1, V2, S;
	double X;
	if (phase == 0) {
		do {
			double U1 = (double)rand() / (double)RAND_MAX;
			double U2 = (double)rand() / (double)RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);
	phase = 1 - phase;

	return X;
}

double gauss(double a, double b)
{
	return a + gaussian() * b;
}

double max_gene(int pos)
{
	double d = 0;
	for (int i = 0; i < POPSIZE; i++)
		if (population[i].x[pos] > d) d = population[i].x[pos];
	return d;
}

double min_gene(int pos)
{
	double d = 1;
	for (int i = 0; i < POPSIZE; i++)
		if (population[i].x[pos] < d) d = population[i].x[pos];
	return d;
}

void update_fbest()
{
	for (int i = 0; i < POPSIZE; i++)
		if (population[i].f[population[i].tao] < fbest[population[i].tao]) fbest[population[i].tao] = population[i].f[population[i].tao];
}

void initialize()
{
	int i, j;
	for (i = 0; i < task_num; i++) fbest[i] = 1e20;

	for(i = 0 ; i < POPSIZE * 2; i++){
		for (j = 0; j < MAXNVARS; j++){
			population[i].x[j] = randval(GLB, GUB);		
		}
		if (i < POPSIZE){
			for (j = 0; j < task_num; j++) {
				population[i].f[j] = objective(population[i].x, j);
			}
		}
	}
	
	
	//rank all objectives .....................
    for(i = 0; i < task_num; i++){
        for(j = 0; j < POPSIZE; j++){
            vis[j].I = j;
            vis[j].v = population[j].f[i];
        }
        qsort(vis, POPSIZE, sizeof(vis[0]), cmp_vi);	//将一个task中的个体从小到大排
        for(j = 0; j < POPSIZE; j++){
            population[vis[j].I].rank[i] = j;			//记录每个个体在task_i上的排名
        }
    }

    //calcualte skill-factor.
    for(i = 0; i < POPSIZE; i++){
		population[i].tao = 0;
		population[i].scale_f = population[i].rank[0];
		for (j = 1; j < task_num; j++){
            if(population[i].rank[j] < population[i].rank[population[i].tao]){
                population[i].tao = j;
                population[i].scale_f = population[i].rank[j];
            }
        }
    }

}	

void SBX(int a, int b)
{
	int i;
	double rand;
	double y1, y2, yl, yu;
	double c1, c2;
	double alpha, beta, betaq;
	if (randval(0, 1) < SBX_rate)
	{
		for (i = 0; i < MAXNVARS; i++)
		{
			if (randval(0, 1) < 0.9)	//交叉率
			{
				if (fabs(population[a].x[i] - population[b].x[i]) > 0)	//如果两个不同才交叉
				{
					if (population[a].x[i] < population[b].x[i])
					{
						y1 = population[a].x[i];
						y2 = population[b].x[i];
					}
					else
					{
						y1 = population[b].x[i];
						y2 = population[a].x[i];
					}
					yl = MIN_GENE[i];
					yu = MAX_GENE[i];
					rand = randval(0, 1);
					beta = 1.0 + (2.0*(y1 - yl) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(mu + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (mu + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (mu + 1.0)));
					}
					c1 = 0.5*((y1 + y2) - betaq*(y2 - y1));
					beta = 1.0 + (2.0*(yu - y2) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(mu + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (mu + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (mu + 1.0)));
					}
					c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));
					if (c1<yl)
						c1 = yl;
					if (c2<yl)
						c2 = yl;
					if (c1>yu)
						c1 = yu;
					if (c2>yu)
						c2 = yu;
					if (randval(0, 1) <= 0.5)
					{
						population[cnt].x[i] = c2;
						population[cnt + 1].x[i] = c1;
					}
					else
					{
						population[cnt].x[i] = c1;
						population[cnt + 1].x[i] = c2;
					}
				}
				else
				{
					population[cnt].x[i] = population[a].x[i];
					population[cnt + 1].x[i] = population[b].x[i];
				}
			}
			else
			{
				population[cnt].x[i] = population[a].x[i];
				population[cnt + 1].x[i] = population[b].x[i];
			}
		}
	}
	else
	{
		for (i = 0; i < MAXNVARS; i++)
		{
			population[cnt].x[i] = population[a].x[i];
			population[cnt + 1].x[i] = population[b].x[i];
		}
	}

	//---------vertical cultral transmission--------------
	if (randval(0, 1) < 0.5){
		population[cnt].tao = population[a].tao;
		population[cnt + 1].tao = population[b].tao;
	}
	else{
		population[cnt].tao = population[b].tao;
		population[cnt + 1].tao = population[a].tao;
	}

	population[cnt].f[population[cnt].tao] = objective(population[cnt].x, population[cnt].tao);
	for (i = 0; i < task_num; i++)
		if (i != population[cnt].tao) population[cnt].f[i] = 1e20;
	cnt++;
	population[cnt].f[population[cnt].tao] = objective(population[cnt].x, population[cnt].tao);
	for (i = 0; i < task_num; i++)
		if (i != population[cnt].tao) population[cnt].f[i] = 1e20;
	cnt++;
}

void gaussian_mutation(int a)
{
	double p = 1. / MAXNVARS;
	for (int i = 0; i < MAXNVARS; i++){
		if (randval(0, 1) < p){
			double t = population[a].x[i] + gauss(0, 0.1 * MAX_GENE[i]);
			if (t > GUB) t = randval(population[a].x[i], GUB);
			if (t < GLB) t = randval(GLB, population[a].x[i]);
			population[cnt].x[i] = t;
		}
	}
	population[cnt].tao = population[a].tao;
	population[cnt].f[population[cnt].tao] = objective(population[cnt].x, population[cnt].tao);
	for (int i = 0; i < task_num; i++)
		if (i != population[cnt].tao) population[cnt].f[i] = 1e20;
	cnt++;
}

void production()
{
	int i, j;	
	
	//Generate POPSIZE offsprings and evaluate them according to their τ
	cnt = POPSIZE;
	for (i = 0; i < POPSIZE / 2; i++) {
		int a = rand() % POPSIZE;
		int b = rand() % POPSIZE;
		if (randval(0, 1) < rmp || population[a].tao == population[b].tao){	//Crossover p_a and p_b
			SBX(a, b);

		} 
		else{	//Mutation to generate two offsprings
			gaussian_mutation(a);
			gaussian_mutation(b);
		}
	}

	//rank all objective .....................
    for(i = 0; i < task_num; i++){
        for(j = 0; j < cnt; j++){
            vis[j].I = j;
            vis[j].v = population[j].f[i];
        }
        qsort(vis, cnt, sizeof(vis[0]), cmp_vi);
        for(j = 0; j < cnt; j++){
            population[vis[j].I].rank[i] = j;
        }
    }

	for (i = 0; i < cnt; i++)
		population[i].scale_f = population[i].rank[population[i].tao];

	//selection
    for(i = 0; i < cnt - 1; i++){
		for (j = i + 1; j < cnt; j++){
			if (population[i].scale_f > population[j].scale_f){
				gene p = population[i];
				population[i] = population[j];
				population[j] = p;
			}
		}
    }
}

void MFEA()
{
	int i, j;
	initialize();
	update_fbest();
	generation = 0;
	//------------record the initial results---------------------
	printf("generation:%d\n", generation);
	for (i = 0; i < task_num; i++) {
		fbest_value[i][0][job] = fbest[i];
		printf("fbest%d=%lf\n", i + 1, fbest_value[i][0][job]);
	}
	
	while(generation < MAXGENS){
		for (int i = 0; i < MAXNVARS; i++){
			MAX_GENE[i] = max_gene(i);
			MIN_GENE[i] = min_gene(i);
		}
		production();

		if ((generation + 1) % RECORD_FRE == 0){
			update_fbest();
			int k = (generation + 1) / RECORD_FRE;
			for (i = 0; i < task_num; i++) fbest_value[i][k][job] = fbest[i];
			printf("generation:%d\n", generation + 1);
			for (i = 0; i < task_num; i++) 
				printf("fbest%d=%lf\n", i + 1, fbest_value[i][k][job]);
			if (generation == 1) system("pause");
		}
		generation++;
	}
}

void final_report()
{
	FILE *pf;
	char name[30];
	double sum[total_task_num][MAXGENS / RECORD_FRE + 1];
	int i, j, l;

	for (i = 0; i < task_num; i++){
		sprintf(name, "Result\\30 RUNS results for task %d.txt", i + 1);
		fopen_s(&pf, name, "w");
		for (j = 0; j < RUNS; j++){
			if (j < RUNS - 1) fprintf(pf, "%lf,", fbest_value[i][MAXGENS / RECORD_FRE][j]);
			else  fprintf(pf, "%lf", fbest_value[i][MAXGENS / RECORD_FRE][j]);
		}
	}
	fclose(pf);

	for (i = 0; i < task_num; i++)
		for (j = 0; j <= MAXGENS / RECORD_FRE; j++)
			sum[i][j] = 0;


	for (i = 0; i <= MAXGENS / RECORD_FRE; i++)
	{
		for (j = 0; j < RUNS; j++){
			for (l = 0; l < task_num; l++)	{
				sum[l][i] += fbest_value[l][i][j];
			}
		}
	}

	fopen_s(&pf, "Result\\final results.txt", "w");
	for (i = 0; i < task_num; i++){
		fprintf(pf, "task%d:%lf\n", i + 1, sum[i][MAXGENS / RECORD_FRE] / RUNS);
	}
	fclose(pf);

	sprintf(name, "Result\\Convergence Curve.txt");
	fopen_s(&pf, name, "w");
	for (i = 0; i <= MAXGENS / RECORD_FRE; i++){
		for (j = 0; j < task_num; j++)
			fprintf(pf, "%lf\t", sum[j][i] / RUNS);
		fprintf(pf, "\n");

	}
	fclose(pf);

}

int main()
{
	task_num = 10;
	for (job = 0; job < RUNS; job++) {
			srand(job); 
			MFEA();
	}
	final_report();
}