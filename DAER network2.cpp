#include<iostream>
#include<fstream>
#include<iomanip>
#include<time.h>
#include<cmath>
#include<cstdlib>
#include"big.c"
using namespace std;
//const int N=10000000,ave_degree_k=3,number=9;
class Link
{
public:
	int num;
	Link *next;
};
class Attribute
{
public:
	int degree_AA,degree_DA;
};
int main()
{
	int N=100000,ave_degree_k=4,number;
	double r,nn;
	cout<<"enter r and nn : ";
	cin>>r>>nn;
	char fname[100];
	sprintf(fname,"DAER2_r_%2.3f_n_%2.3f.net",r,nn);
	ofstream outfile(fname);
	if(!outfile)
	{
		cerr<<"error!";
		exit(1);
	}
	static Attribute *indiv;
	indiv=new Attribute[N];
	rnd_seed((unsigned long)time(0));
	Link **head,**tail,**move,**add;
	head=new Link *[N];
	tail=new Link *[N];
	move=new Link *[N];
	add=new Link *[N];
	int A=(1-r)*N;
	int D=N-A;
	int linkmax=ave_degree_k*A/2,m,n;
	int linkmax2=nn*A;
	for(int i=0;i<N;i++)
	{
		indiv[i].degree_AA=indiv[i].degree_DA=0;
		add[i]=new Link;
		add[i]->next=NULL;
		add[i]->num=i;
		head[i]=move[i]=tail[i]=add[i];
	}
	do
	{
		bool flag;
		do
		{
			flag=false;
			m=rnd()%A;
			n=rnd()%A;
			if(m==n)
				flag=true;
			else 
			{
				while(move[m])
				{
					if(move[m]->num==n)
					{
						flag=true;
						break;
					}
					else 
						move[m]=move[m]->next;
				}
				move[m]=head[m];
			}
		}while(flag);
		add[m]=new Link;
		add[m]->next=NULL;
		add[m]->num=n;
		tail[m]->next=add[m];
		tail[m]=add[m];
		indiv[m].degree_AA++;
		add[n]=new Link;
		add[n]->next=NULL;
		add[n]->num=m;
		tail[n]->next=add[n];
		tail[n]=add[n];
		indiv[n].degree_AA++;
		linkmax--;
	}while(linkmax);

	while(linkmax2)
	{
		bool flag;
		do
		{
			flag=false;
			m=A+rnd()%D;
			n=rnd()%A;
			while(move[m])
			{
				if(move[m]->num==n)
				{
					flag=true;
					break;
				}
				else 
					move[m]=move[m]->next;
			}
			move[m]=head[m];
		}while(flag);
		add[m]=new Link;
		add[m]->next=NULL;
		add[m]->num=n;
		tail[m]->next=add[m];
		tail[m]=add[m];
		indiv[m].degree_DA++;
		add[n]=new Link;
		add[n]->next=NULL;
		add[n]->num=m;
		tail[n]->next=add[n];
		tail[n]=add[n];
		indiv[n].degree_DA++;
		linkmax2--;
	}

	outfile<<N<<" ";
	int kmax=0;
	for(int i=0;i<N;i++)
		if((indiv[i].degree_AA+indiv[i].degree_DA)>kmax)
			kmax=indiv[i].degree_AA+indiv[i].degree_DA;
	outfile<<kmax<<endl;
	cout<<kmax<<endl;
	for(int i=0;i<N;i++)
	{
		outfile<<indiv[i].degree_AA<<" "<<indiv[i].degree_DA<<" ";
		move[i]=move[i]->next;
		while(move[i])
		{
			outfile<<move[i]->num<<" ";
			move[i]=move[i]->next;
		}
		outfile<<endl;
	}
	outfile.close();
	return 0;
}