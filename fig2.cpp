#include<iostream>
#include<fstream>
#include<iomanip>
#include<time.h>
#include<cmath>
#include<cstdlib>
#include"big.c"
using namespace std;

class Attribute
{
public:
	int degree_AA,degree_DA,state,neigh_I;
};

const int N=100000,seeds=1000,min_=0,max_=9,count=100;
int main()//int argc,char**argv
{
	time_t begin,end;
    begin=clock();
	int ave_degree_k=4,number;
	double beta=0.12,mu0=0.5,alpha=1,omega=1,r=0.01,nn;
	cout<<"enter na :";
	cin>>nn;
	rnd_seed((unsigned long)time(0));
	static Attribute *indiv;
	indiv=new Attribute[N];
	static int check,**neighbor,*ss,*num_I,kmax;
	neighbor=new int *[N];
	ss=new int [N];
	num_I=new int [N];
	static double *R,*mu,*betaI_temp;
	R=new double [N];
	mu=new double [N];
	double u,v;
	int m,n;
	char fname[100];
	sprintf(fname,"fig12-0.01seed_alpha_%2.3f_omega_%2.3f_r_%2.3f_n_%2.3f.dat",alpha,omega,r,nn);
	ofstream outfile(fname);
	if(!outfile)
	{
		cerr<<"error!";
		exit(1);
	}
	char fname1[100];
	sprintf(fname1,"DAER2_r_%2.3f_n_%2.3f.net",r,nn);
	ifstream infile(fname1);
	if(!infile)
	{
		cerr<<"error!";
		exit(1);
	}
	infile>>check;
	cout<<"check="<<check<<" "<<"N="<<N<<endl;
	infile>>kmax;//max-degree
	betaI_temp=new double [kmax+1];
	for(int i=0;i<=kmax;i++)
		betaI_temp[i]=1-pow((1-beta),i);
	for(int i=0;i<N;i++)
	{
		int temp;
		infile>>temp;
		indiv[i].degree_AA=temp;
		infile>>temp;
		indiv[i].degree_DA=temp;
		int degree_temp=indiv[i].degree_AA+indiv[i].degree_DA;
		neighbor[i]=new int [degree_temp];
		for(int j=0;j<degree_temp;j++)
			infile>>neighbor[i][j];
	}
	infile.close();

	double temp_R=mu0/(1+1/alpha*exp(-omega/r));
	double result[count]={0};
	int temp_count=0;
	while(temp_count<count)
	{cout<<temp_count<<endl;
		int A=(1-r)*N;
		int S=N,I=0;
		for(int i=0;i<N;i++)
		{
			indiv[i].state=1;
			ss[i]=i;
			indiv[i].neigh_I=0;
		}
		for(int i=0;i<seeds;i++)
		{
			m=rnd()%S;
			int temp=ss[m];
			indiv[temp].state=2;
			S--;
			I++;
			ss[m]=ss[S];
			if(temp<A)
			{
				for(int j=0;j<indiv[temp].degree_AA+indiv[temp].degree_DA;j++)
				{
					int neighbor_temp=neighbor[temp][j];
					indiv[neighbor_temp].neigh_I++;
				}
			}
		}
		int t=0,t_sum=900,num_sum=0;
		double re_temp=0.0;
		while(I&&t<1000)
		{
			for(int i=0;i<N;i++)
			{
				num_I[i]=indiv[i].neigh_I;
				R[i]=0;
			}
			for(int i=A;i<N;i++)
			{
				if(indiv[i].state==2)	mu[i]=temp_R;
				else
				{
					for(int j=0;j<indiv[i].degree_DA;j++)
						if(indiv[neighbor[i][j]].state==2) R[neighbor[i][j]]+=1.0/(r*num_I[i]);
				}
			}
			for(int i=0;i<A;i++)
				mu[i]=mu0/(1+1/alpha*exp(-omega*R[i]));
			for(int i=0;i<N;i++)
			{
				if(indiv[i].state==1)
				{
					double beta_I=betaI_temp[num_I[i]];
					v=drnd();
					if(v<beta_I)
					{
						indiv[i].state=2;
						I++;
						S--;
						if(i<A)
						{
							for(int j=0;j<indiv[i].degree_AA+indiv[i].degree_DA;j++)
							{
								int neighbor_i=neighbor[i][j];
								indiv[neighbor_i].neigh_I++;
							}
						}
					}
				}
				else if(indiv[i].state==2)
				{
					v=drnd();
					if(v<mu[i])
					{
						indiv[i].state=1;
						S++;
						I--;
						if(i<A)
						{
							for(int j=0;j<indiv[i].degree_AA+indiv[i].degree_DA;j++)
							{
								int neighbor_i=neighbor[i][j];
								indiv[neighbor_i].neigh_I--;
							}
						}
					}
				}
			}
			t++;
			if(t>t_sum)
			{
				re_temp+=double(I)/N;
				num_sum++;
			}
		}
		if(num_sum!=0) result[temp_count]=re_temp/num_sum;
		else result[temp_count]=0;
		temp_count++;
	}

    double result_ave=0,result_chi,temp1=0;
    for(int i=0;i<count;i++)
    	result_ave+=result[i]/count;
    for(int i=0;i<count;i++)
    	temp1 += (result[i]-result_ave)*(result[i]-result_ave);
    result_chi=sqrt(temp1/(count-1));
    outfile<<nn<<" "<<result_ave<<" "<<result_chi<<endl;
    
	outfile.close();
	end=clock();
    cout<<"runtime: "<<double(end-begin)/CLOCKS_PER_SEC<<endl;
	return 0;
}