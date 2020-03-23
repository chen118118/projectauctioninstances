
#include <iostream>
#include <fstream>

using namespace std;

#define NN  10		//max number of task nodes in a project graph 
#define NA	5		//max number of pre/post nodes of each task node in a project graph
#define NP  20		//max number of task paths in a project graph 
#define ND	100		//max duration of each task
#define NB	10		//max number of bidders of each task

struct task
{
	int index;

	int npre; 
	int npost;
	int pre[NA];
	int post[NA];

	int dmin;	//minimum quotable duration 
	int dmax;	 //maximum quotable duration

	double wmin;	//minimum workload
	double wmax;	//maximum workload
	double wmean;	//mean work load
} ;


struct bidder
{
	int index; //index

	int task;  //index of the task that can be bid for

	double k;  //type

	int c[ND]; //execution cost of the task for each quotable duration
}; 


std::ofstream df;   //data file


void main()
{
	void instancegenerator(int projectnetwork, int workloadscenario, int instance);

	int projectnetwork, workloadscenario, instance;

	df.open("projectauction-instancedata.txt", ios::app);

	for (projectnetwork = 1; projectnetwork <= 5; projectnetwork++)
	for (workloadscenario = 1; workloadscenario <= 3; workloadscenario++)
	for (instance = 1; instance <= 10; instance++)
	{

		instancegenerator(projectnetwork, workloadscenario, instance);

	}
	
	df.close();
}


void instancegenerator(int projectnetwork, int workloadscenario, int instance)
{

	double uniform01();
	double taskexecutioncost(double ci, double cd, double k, double wmin, double wmax, int d);

	int nn;  //number of tasks (nodes)
	int nb[NN]; //number of bidders of each task

	double cu;  //user cost of each project

	double rotask[NN];  //hidden user cost of each task 
	double citask[NN], cdtask[NN];  //incentive and disincentive of each task

	struct task task[NN];
	struct bidder bidder[NN][NB];   // bidder[i][m] denotes the m-th bidder for task i

	int imin[NN]; //best bidder for each task
	double kmin;

	int imax[NN]; //worst bidder for each task
	double kmax;

	double weight[NN]; //weight of each task in the project
	
	int par; 

	int v;
	double vd;

	int i, j, m;

	df << "Project Network " << projectnetwork << ", " << "Workload Scenario " << workloadscenario << ", " << "Instance " << instance << ": " << endl;

	//input the network structure for each of the 5 instance project sets

	//project instance (a)

	if (projectnetwork == 1)
	{
		nn = 3;

		task[0].npre = 0; task[0].npost = 0;
		task[1].npre = 0; task[1].npost = 0;
		task[2].npre = 0; task[2].npost = 0;

		weight[0] = 1.0 / 3.0; weight[1] = 1.0 / 3.0; weight[2] = 1.0 / 3.0;


	}

	//project instance (b)
	
	if (projectnetwork == 2)
	{

		nn = 4;

		task[0].npre = 0; task[0].npost = 1;
		task[0].post[0] = 3;

		task[1].npre = 0; task[1].npost = 1;
		task[1].post[0] = 3;

		task[2].npre = 0; task[2].npost = 0;

		task[3].npre = 2; task[3].npost = 0;
		task[3].pre[0] = 0; task[3].pre[1] = 1;

		weight[0] = 1.0 / 3.0; weight[1] = 1.0 / 3.0; weight[2] = 1.0 / 3.0; weight[3] = 2.0 / 3.0;
	}


	//project instance (c)

	if (projectnetwork == 3)
	{

		nn = 5;

		task[0].npre = 0; task[0].npost = 2;
		task[0].post[0] = 2; task[0].post[1] = 3;

		task[1].npre = 0; task[1].npost = 1;
		task[1].post[0] = 4;

		task[2].npre = 1; task[2].npost = 1;
		task[2].pre[0] = 0;
		task[2].post[0] = 4;

		task[3].npre = 1; task[3].npost = 0;
		task[3].pre[0] = 0;

		task[4].npre = 2; task[4].npost = 0;
		task[4].pre[0] = 1; task[4].pre[1] = 2;

		weight[0] = 2.0 / 3.0; weight[1] = 1.0 / 3.0; weight[2] = 1.0 / 3.0; weight[3] = 1.0 / 3.0; weight[4] = 2.0 / 3.0;
	}

	//project instance (d)

	if (projectnetwork == 4)
	{
		nn = 6;

		task[0].npre = 0; task[0].npost = 2;
		task[0].post[0] = 2; task[0].post[1] = 3;

		task[1].npre = 0; task[1].npost = 1;
		task[1].post[0] = 4;

		task[2].npre = 1; task[2].npost = 1;
		task[2].pre[0] = 0;
		task[2].post[0] = 5;

		task[3].npre = 1; task[3].npost = 0;
		task[3].pre[0] = 0;

		task[4].npre = 1; task[4].npost = 0;
		task[4].pre[0] = 1;

		task[5].npre = 1; task[5].npost = 0;
		task[5].pre[0] = 2;

		weight[0] = 2.0 / 3.0; weight[1] = 1.0 / 3.0; weight[2] = 1.0 / 3.0; weight[3] = 1.0 / 3.0; weight[4] = 1.0 / 3.0; weight[5] = 1.0 / 3.0;
	}

	//project instance (e)
	
	if (projectnetwork == 5)
	{

		nn = 7;

		task[0].npre = 0; task[0].npost = 1;
		task[0].post[0] = 2;

		task[1].npre = 0; task[1].npost = 2;
		task[1].post[0] = 3; task[1].post[1] = 4;

		task[2].npre = 1; task[2].npost = 1;
		task[2].pre[0] = 0;
		task[2].post[0] = 5;

		task[3].npre = 1; task[3].npost = 1;
		task[3].pre[0] = 1;
		task[3].post[0] = 5;

		task[4].npre = 1; task[4].npost = 1;
		task[4].pre[0] = 1;
		task[4].post[0] = 6;

		task[5].npre = 2; task[5].npost = 0;
		task[5].pre[0] = 2; task[5].pre[1] = 3;


		task[6].npre = 1; task[6].npost = 0;
		task[6].pre[0] = 4;

		weight[0] = 1.0 / 3.0; weight[1] = 2.0 / 3.0; weight[2] = 1.0 / 3.0; weight[3] = 1.0 / 3.0; weight[4] = 1.0 / 3.0; weight[5] = 2.0 / 3.0; weight[6] = 1.0 / 3.0;
	}


	df << "  number of tasks: " << nn << endl;

	df << "  weight of each task: ";
	for (i = 0; i < nn; i++) df << weight[i] << ", ";
	df << endl;

	//set the number of bidders for each task

	for (i = 0; i < nn; i++) nb[i] = 3;

	df << "  number of bidders of each task: ";
	for (i = 0; i < nn; i++) df << nb[i] << ", ";
	df << endl;

	//generate the type of each bidder for each task
	for (i = 0; i < nn; i++)
	for (m = 0; m < nb[i]; m++)
	{
		bidder[i][m].index = m;
		bidder[i][m].task = i;
		bidder[i][m].k = uniform01() * 1 + 1;
	}

	df << "  types of the bidders for each task: " << endl;

	for (i = 0; i < nn; i++)
	{
		df << "    bidders for task " << i+1 << ": ";
		for (m = 0; m < nb[i]; m++)
		{
			df << bidder[i][m].k << ", ";
		}
		df << endl;
	}

	//identify the best bidder of each task

	for (i = 0; i < nn; i++)
	{
		kmin = bidder[i][0].k;
		imin[i] = 0;

		for (m = 1; m < nb[i]; m++)
		{
			if (bidder[i][m].k < kmin)
			{
				kmin = bidder[i][m].k;
				imin[i] = m;
			}

		}
	}


	//identify the worst bidder of each task

	for (i = 0; i < nn; i++)
	{

		kmax = bidder[i][0].k;
		imax[i] = 0;

		for (m = 1; m < nb[i]; m++)
		{

			if (bidder[i][m].k > kmax)
			{

				kmax = bidder[i][m].k;
				imax[i] = m;
			}

		}
	}


	//generate the mean workload, minimum workload and maximum workload of each task

	for (i = 0; i < nn; i++)
	{

		task[i].wmean = uniform01() * 10 + 5;

		if (workloadscenario == 1)
		{
			task[i].wmin = task[i].wmean;
			task[i].wmax = task[i].wmean;
		}

		if (workloadscenario == 2)
		{
			task[i].wmin = task[i].wmean*0.9;
			task[i].wmax = task[i].wmean*1.1;
		}

		if (workloadscenario == 3)
		{
			task[i].wmin = task[i].wmean*0.8;
			task[i].wmax = task[i].wmean*1.2;
		}

	}

	df << "  mean, minimum, and maximum work load of each task: " << endl;
	for (i = 0; i < nn; i++)
	{
		df << "    task " << i + 1 << ": " << task[i].wmean << ", " << task[i].wmin << ", " << task[i].wmax << endl;
	}

	//set cu for the project and ro, ci, cd for each task

	cu = 10;

	df << "  user cost of the project: " << cu << endl;


	//set the minimum due date and maximum due date of each task

	for (i = 0; i < nn; i++)
	{
		task[i].dmin = 0;

		//task[i].dmax = int(2 * sqrt(bidder[i][imin[i]].k / (weight[i] * cu))*task[i].wmax) + 1;
		task[i].dmax = int(2 * sqrt(bidder[i][imax[i]].k / (weight[i] * cu))*task[i].wmax) + 1;
	}


	df << "  minimum and maximum duration of each task: " << endl;
	for (i = 0; i < nn; i++)
	{
		df << "    task " << i + 1 << ": " << task[i].dmin << ", " << task[i].dmax << endl;
	}

	df << endl;


	for (par = 1; par <= 9; par++)
	{
		if (par == 1)
		{

			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = weight[i] * cu;
				cdtask[i] = weight[i] * cu;
			}


		}

		if (par == 2)
		{
			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = 0.5 * weight[i] * cu;
				cdtask[i] = weight[i] * cu;
			}

		}

		if (par == 3)
		{
			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = 0;
				cdtask[i] = weight[i] * cu;
			}

		}


		if (par == 4)
		{

			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = weight[i] * cu;
				cdtask[i] = 0.5 *(weight[i] + 1) * cu;
			}


		}

		if (par == 5)
		{

			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = 0.5 * weight[i] * cu;
				cdtask[i] = 0.5 *(weight[i] + 1) * cu;
			}

		}

		if (par == 6)
		{

			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = 0;
				cdtask[i] = 0.5 *(weight[i] + 1) * cu;
			}

		}


		if (par == 7)
		{

			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = weight[i] * cu;
				cdtask[i] = cu;
			}


		}

		if (par == 8)
		{

			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = 0.5 * weight[i] * cu;
				cdtask[i] = cu;
			}

		}

		if (par == 9)
		{

			for (i = 0; i < nn; i++)
			{

				rotask[i] = weight[i] * cu;
				citask[i] = 0;
				cdtask[i] = cu;
			}

		}


		//calculate the cost function for each task

		for (i = 0; i < nn; i++)
		for (m = 0; m < nb[i]; m++)
		{

			for (j = task[i].dmin; j <= task[i].dmax; j++)
			{
				vd = taskexecutioncost(citask[i], cdtask[i], bidder[i][m].k, task[i].wmin, task[i].wmax, j);
				v = int(vd);
				if (vd - v < 0.5) bidder[i][m].c[j - task[i].dmin] = v;
				else bidder[i][m].c[j - task[i].dmin] = v + 1;

			}

		}

	}


}



#define max(x,y) (((x)>(y))?(x):(y))
#define min(x,y) (((x)<(y))?(x):(y))


double uniform01()
{

	return rand() / 32768.0;

}



//calcualte the integral of f(w) from a to x, where f(w) is the density function of a PERT distributed random variable

double F(double a, double c, double x)
{

	double v;

	if (fabs(c - a) < 0.000000001)
	{
		if (a <= x && x <= c) return 1;
		else
		{
			cout << "error in function F()!" << endl;
			exit(1);
		}

	}


	v = 10 * pow(x - a, 3)*pow(c - x, 2) + 5 * pow(x - a, 4)*(c - x) + pow(x - a, 5);

	v = v / (pow(c - a, 5));

	return v;

}

//calcualte the integral of wf(w) from a to x, where f(w) is the density function of a PERT distributed random variable
double E(double a, double c, double x)
{

	double F(double a, double c, double x);

	double v;


	if (fabs(c - a) < 0.000000001)
	{
		if (a <= x && x <= c) return a;
		else
		{
			cout << "error in function E()!" << endl;
			exit(1);
		}

	}

	v = 7.5 * pow(x - a, 4)*pow(c - x, 2) + 3 * pow(x - a, 5)*(c - x) + 0.5*pow(x - a, 6);

	v = v / (pow(c - a, 5));

	v = a * F(a, c, x) + v;

	return v;
}


//calcualte the integral of w*w*f(w) from a to x, where f(w) is the density function of a PERT distributed random variable
double M(double a, double c, double x)
{
	;
	double F(double a, double c, double x);
	double E(double a, double c, double x);

	double v;


	if (fabs(c - a) < 0.000000001)
	{
		if (a <= x && x <= c) return a * a;
		else
		{
			cout << "error in function M()!" << endl;
			exit(1);
		}

	}

	v = 6 * pow(x - a, 5)*pow(c - x, 2) + 2 * pow(x - a, 6)*(c - x) + (2.0 / 7)*pow(x - a, 7);

	v = v / (pow(c - a, 5));

	v = -a * a*F(a, c, x) + 2 * a*E(a, c, x) + v;

	return v;
}


//calculate the expected cost of a company of type k for executing a task
//with a quoted duration d and the workload PERT symmetrically distributed on [wmin, wmax], including incentive/desincentive costs ci and cd
double taskexecutioncost(double ci, double cd, double k, double wmin, double wmax, int d)
{
	double F(double a, double c, double x);
	double E(double a, double c, double x);
	double M(double a, double c, double x);

	double delta, b;
	double w1, w2;
	double v;

	if (fabs(wmax - wmin) < 0.000000001)
	{
		if (sqrt(k / cd)*wmin <= d && d <= sqrt(k / ci)*wmin) v = k * wmin*wmin / d;
		if (d < sqrt(k / cd)*wmin) v = 2 * sqrt(k*cd)*wmin - cd * d;
		if (d > sqrt(k / ci)*wmin) v = 2 * sqrt(k*ci)*wmin - ci * d;

		return v;

	}


	delta = wmax - wmin;
	b = (wmin + wmax) / 2;

	if (d == 0) v = 2 * sqrt(k*cd)*b;  //???

	if (d > 0)
	{
		w1 = max(wmin, min(wmax, sqrt(ci / k)*d));
		w2 = min(wmax, max(wmin, sqrt(cd / k)*d));


		v = -ci * d*F(wmin, wmax, w1) - cd * d*(1 - F(wmin, wmax, w2));
		v = v + 2 * sqrt(k*ci)*E(wmin, wmax, w1) + 2 * sqrt(k*cd)*(b - E(wmin, wmax, w2));
		v = v + (k / d) * (M(wmin, wmax, w2) - M(wmin, wmax, w1));
	}

	return v;

}
